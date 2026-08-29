#include "Leiden.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

constexpr double kTolerance = 1.0e-10;

class UnsupportedSparseCPM final : public QualityFunction {
public:
    explicit UnsupportedSparseCPM(double gamma) : delegate_(gamma) {}
    double quality(const Graph& G, const LeidenGraphStats& stats,
                   const LeidenPartition& partition) const override {
        return delegate_.quality(G, stats, partition);
    }
    double deltaMove(const Graph& G, const LeidenGraphStats& stats,
                     const LeidenPartition& partition, Vertex v,
                     Community target) const override {
        return delegate_.deltaMove(G, stats, partition, v, target);
    }
    double deltaMoveFromWeights(const LeidenGraphStats& stats,
                                const LeidenPartition& partition, Vertex v,
                                Community target, double source_weight,
                                double target_weight) const override {
        return delegate_.deltaMoveFromWeights(stats, partition, v, target,
                                              source_weight, target_weight);
    }
    double refinementNodeMass(const LeidenGraphStats& stats,
                              Vertex v) const override {
        return delegate_.refinementNodeMass(stats, v);
    }
    double refinementResolution(const LeidenGraphStats& stats) const override {
        return delegate_.refinementResolution(stats);
    }
private:
    CPMQualityFunction delegate_;
};

void Fail(const std::string& test_name,
          const std::string& expected,
          const std::string& actual,
          const std::string& difference)
{
    std::cerr << test_name << " failed\n"
              << "expected value: " << expected << "\n"
              << "actual value:   " << actual << "\n"
              << "difference:     " << difference << "\n";
    std::exit(EXIT_FAILURE);
}

void CheckTrue(const std::string& test_name, bool condition)
{
    if (!condition) {
        Fail(test_name, "true", "false", "condition is false");
    }
}

void CheckNear(const std::string& test_name, double expected, double actual)
{
    const double diff = std::abs(expected - actual);
    if (diff > kTolerance) {
        Fail(test_name,
             std::to_string(expected),
             std::to_string(actual),
             std::to_string(diff));
    }
}

void CheckThrowsInvalidArgument(const std::string& test_name,
                                const Graph& G,
                                const LeidenGraphStats& stats,
                                const LeidenPartition& partition,
                                const QualityFunction& qf,
                                double theta)
{
    std::mt19937 rng(1);
    try {
        (void)RefinePartition(G, stats, partition, qf, theta, rng);
    } catch (const std::invalid_argument&) {
        return;
    }
    Fail(test_name, "std::invalid_argument", "no exception", "theta accepted");
}

void CheckValidPartition(const std::string& test_name,
                         const Graph& G,
                         const LeidenPartition& partition)
{
    CheckTrue(test_name + " assignment size",
              static_cast<int>(partition.community_of.size()) == num_vertices(G));
    const std::size_t nc = partition.community_size.size();
    CheckTrue(test_name + " statistic sizes",
              partition.community_strength.size() == nc &&
              partition.internal_edge_weight.size() == nc);

    for (Vertex v = 0; v < num_vertices(G); ++v) {
        const Community c = partition.community_of[v];
        CheckTrue(test_name + " assigned", c >= 0);
        CheckTrue(test_name + " in range", static_cast<std::size_t>(c) < nc);
    }
}

void CheckStatsAgainstRecompute(const std::string& test_name,
                                const Graph& G,
                                const LeidenGraphStats& stats,
                                const LeidenPartition& partition)
{
    const LeidenPartition recomputed =
        MakePartition(G, stats, partition.community_of);
    const std::size_t nc = std::max(partition.community_size.size(),
                                    recomputed.community_size.size());

    for (std::size_t c = 0; c < nc; ++c) {
        const double expected_size =
            (c < recomputed.community_size.size()) ? recomputed.community_size[c] : 0.0;
        const double actual_size =
            (c < partition.community_size.size()) ? partition.community_size[c] : 0.0;
        const double expected_strength =
            (c < recomputed.community_strength.size()) ? recomputed.community_strength[c] : 0.0;
        const double actual_strength =
            (c < partition.community_strength.size()) ? partition.community_strength[c] : 0.0;
        const double expected_internal =
            (c < recomputed.internal_edge_weight.size()) ? recomputed.internal_edge_weight[c] : 0.0;
        const double actual_internal =
            (c < partition.internal_edge_weight.size()) ? partition.internal_edge_weight[c] : 0.0;

        CheckNear(test_name + " community_size", expected_size, actual_size);
        CheckNear(test_name + " community_strength", expected_strength, actual_strength);
        CheckNear(test_name + " internal_edge_weight", expected_internal, actual_internal);
    }
}

void MarkSubsetForTest(const Graph& G,
                       const std::vector<Vertex>& subset,
                       std::vector<std::size_t>& subset_mark,
                       std::size_t& subset_generation)
{
    if (subset_generation == std::numeric_limits<std::size_t>::max()) {
        std::fill(subset_mark.begin(), subset_mark.end(), 0);
        subset_generation = 1;
    } else {
        ++subset_generation;
    }
    for (Vertex v : subset) {
        ValidateVertex(v, G);
        subset_mark[v] = subset_generation;
    }
}

bool IsInSubsetForTest(const std::vector<std::size_t>& subset_mark,
                       std::size_t subset_generation,
                       Vertex v)
{
    return subset_mark[v] == subset_generation;
}

RefinementCommunityStats DirectRefinementCommunityStats(
    const Graph& G,
    const LeidenGraphStats& stats,
    const LeidenPartition& refined,
    const QualityFunction& qf,
    const std::vector<Vertex>& subset,
    const std::vector<std::size_t>& subset_mark,
    std::size_t subset_generation)
{
    RefinementCommunityStats direct;
    direct.entries.reserve(subset.size());

    for (Vertex v : subset) {
        const Community community = refined.community_of[v];
        RefinementCommunityEntry& entry = direct.entries[community];
        ++entry.member_count;
        entry.mass += qf.refinementNodeMass(stats, v);
    }

    for (const auto& item : direct.entries) {
        const Community community = item.first;
        for (Vertex v : subset) {
            if (refined.community_of[v] != community) {
                continue;
            }
            for (const Edge& e : G.adj[v]) {
                if (e.to != v &&
                    subset_mark[e.to] == subset_generation &&
                    refined.community_of[e.to] != community) {
                    direct.entries[community].external_weight += e.weight;
                }
            }
        }
    }

    direct.active_communities.reserve(direct.entries.size());
    for (const auto& item : direct.entries) {
        if (item.second.member_count > 0) {
            direct.active_communities.push_back(item.first);
        }
    }
    std::sort(direct.active_communities.begin(),
              direct.active_communities.end());

    return direct;
}

RefinementCommunityEntry GetEntry(const RefinementCommunityStats& stats,
                                  Community community)
{
    const auto it = stats.entries.find(community);
    if (it == stats.entries.end()) {
        return {};
    }
    return it->second;
}

void CheckEntryNear(const std::string& test_name,
                    const RefinementCommunityStats& stats,
                    Community community,
                    int member_count,
                    double mass,
                    double external_weight)
{
    const RefinementCommunityEntry entry = GetEntry(stats, community);
    CheckTrue(test_name + " member_count",
              entry.member_count == member_count);
    CheckNear(test_name + " mass", mass, entry.mass);
    CheckNear(test_name + " external_weight",
              external_weight,
              entry.external_weight);
}

void CheckRefinementCommunityStatsEqual(
    const std::string& test_name,
    const RefinementCommunityStats& expected,
    const RefinementCommunityStats& actual)
{
    for (const auto& item : expected.entries) {
        const RefinementCommunityEntry actual_entry =
            GetEntry(actual, item.first);
        CheckTrue(test_name + " member_count",
                  item.second.member_count == actual_entry.member_count);
        CheckNear(test_name + " mass",
                  item.second.mass,
                  actual_entry.mass);
        CheckNear(test_name + " external_weight",
                  item.second.external_weight,
                  actual_entry.external_weight);
    }
    for (const auto& item : actual.entries) {
        if (expected.entries.find(item.first) == expected.entries.end() &&
            item.second.member_count == 0) {
            continue;
        }
        const RefinementCommunityEntry expected_entry =
            GetEntry(expected, item.first);
        CheckTrue(test_name + " reverse member_count",
                  expected_entry.member_count == item.second.member_count);
        CheckNear(test_name + " reverse mass",
                  expected_entry.mass,
                  item.second.mass);
        CheckNear(test_name + " reverse external_weight",
                  expected_entry.external_weight,
                  item.second.external_weight);
    }
    std::vector<Community> expected_active = expected.active_communities;
    std::vector<Community> actual_active = actual.active_communities;
    std::sort(expected_active.begin(), expected_active.end());
    std::sort(actual_active.begin(), actual_active.end());
    CheckTrue(test_name + " active_communities",
              expected_active == actual_active);
}

void CheckRefinementOfPartition(const std::string& test_name,
                                const LeidenPartition& partition,
                                const LeidenPartition& refined)
{
    for (std::size_t u = 0; u < refined.community_of.size(); ++u) {
        for (std::size_t v = u + 1; v < refined.community_of.size(); ++v) {
            if (refined.community_of[u] == refined.community_of[v]) {
                CheckTrue(test_name + " no cross-P merge",
                          partition.community_of[u] == partition.community_of[v]);
            }
        }
    }
}

void CheckRefinementRun(const std::string& test_name,
                        const Graph& G,
                        const LeidenGraphStats& stats,
                        const LeidenPartition& partition,
                        const QualityFunction& qf,
                        double theta,
                        unsigned int seed)
{
    const LeidenPartition singleton = MakeSingletonPartition(G, stats);
    const double quality_before = qf.quality(G, stats, singleton);
    std::mt19937 rng(seed);
    const LeidenPartition refined =
        RefinePartition(G, stats, partition, qf, theta, rng);
    const double quality_after = qf.quality(G, stats, refined);

    CheckValidPartition(test_name, G, refined);
    CheckRefinementOfPartition(test_name, partition, refined);
    CheckStatsAgainstRecompute(test_name, G, stats, refined);
    CheckTrue(test_name + " quality from singleton nondecreasing",
              quality_after + kTolerance >= quality_before);
}

void CheckSamePartition(const std::string& test_name,
                        const LeidenPartition& expected,
                        const LeidenPartition& actual)
{
    CheckTrue(test_name + " community_of",
              expected.community_of == actual.community_of);
    CheckTrue(test_name + " community_size",
              expected.community_size == actual.community_size);
    CheckTrue(test_name + " community_strength",
              expected.community_strength == actual.community_strength);
    CheckTrue(test_name + " internal_edge_weight",
              expected.internal_edge_weight == actual.internal_edge_weight);
}

Graph MakeUnweightedGraph()
{
    Graph G = MakeGraph(6);
    add_undirected_edge(G, 0, 1);
    add_undirected_edge(G, 1, 2);
    add_undirected_edge(G, 2, 0);
    add_undirected_edge(G, 2, 3);
    add_undirected_edge(G, 3, 4);
    add_undirected_edge(G, 4, 5);
    add_undirected_edge(G, 5, 3);
    return G;
}

Graph MakeWeightedGraph()
{
    Graph G = MakeGraph(5);
    add_undirected_edge(G, 0, 1, 2.5);
    add_undirected_edge(G, 1, 2, 0.75);
    add_undirected_edge(G, 2, 3, 3.0);
    add_undirected_edge(G, 3, 4, 1.25);
    add_undirected_edge(G, 0, 4, 4.5);
    add_undirected_edge(G, 1, 4, 0.5);
    return G;
}

Graph MakeSelfLoopGraph()
{
    Graph G = MakeGraph(4);
    add_undirected_edge(G, 0, 0, 2.0);
    add_undirected_edge(G, 0, 1, 1.5);
    add_undirected_edge(G, 1, 2, 0.5);
    add_undirected_edge(G, 2, 2, 3.0);
    add_undirected_edge(G, 2, 3, 2.5);
    return G;
}

void TestRefinementOfMoveNodesFastPartition()
{
    const Graph G = MakeUnweightedGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const CPMQualityFunction cpm(0.4);
    std::mt19937 move_rng(123);
    const MoveNodesFastResult moved =
        MoveNodesFast(G, stats, MakeSingletonPartition(G, stats), cpm, move_rng);

    CheckRefinementRun("Test A refinement of P",
                       G,
                       stats,
                       moved.partition,
                       cpm,
                       0.01,
                       456);
}

void TestNoIllegalCrossMerge()
{
    const Graph G = MakeUnweightedGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition partition =
        MakePartition(G, stats, {0, 0, 0, 1, 1, 1});
    const CPMQualityFunction cpm(0.2);
    const ModularityQualityFunction modularity(1.0);

    CheckRefinementRun("Test B CPM no cross merge", G, stats, partition, cpm, 0.05, 10);
    CheckRefinementRun("Test B modularity no cross merge",
                       G,
                       stats,
                       partition,
                       modularity,
                       0.05,
                       10);
}

void TestWeighted()
{
    const Graph G = MakeWeightedGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition partition =
        MakePartition(G, stats, {0, 0, 1, 1, 1});
    const CPMQualityFunction cpm(0.35);
    const ModularityQualityFunction modularity(0.8);

    CheckRefinementRun("Test C CPM weighted", G, stats, partition, cpm, 0.1, 7);
    CheckRefinementRun("Test C modularity weighted",
                       G,
                       stats,
                       partition,
                       modularity,
                       0.1,
                       7);
}

void TestSelfLoop()
{
    const Graph G = MakeSelfLoopGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition partition =
        MakePartition(G, stats, {0, 0, 1, 1});
    const CPMQualityFunction cpm(0.2);
    const ModularityQualityFunction modularity(1.2);
    std::vector<std::size_t> subset_mark(num_vertices(G), 0);
    std::size_t subset_generation = 0;
    MarkSubsetForTest(G, {0, 1}, subset_mark, subset_generation);

    CheckNear("Test D self-loop excluded from E(v,S-v)",
              1.5,
              EdgeWeightFromNodeToSubset(G,
                                         0,
                                         subset_mark,
                                         subset_generation));
    CheckRefinementRun("Test D CPM self-loop", G, stats, partition, cpm, 0.1, 42);
    CheckRefinementRun("Test D modularity self-loop",
                       G,
                       stats,
                       partition,
                       modularity,
                       0.1,
                       42);
}

void TestThetaReproducibility()
{
    const Graph G = MakeWeightedGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition partition =
        MakePartition(G, stats, {0, 0, 1, 1, 1});
    const CPMQualityFunction cpm(0.35);
    std::mt19937 rng1(2026);
    std::mt19937 rng2(2026);

    const LeidenPartition a =
        RefinePartition(G, stats, partition, cpm, 0.05, rng1);
    const LeidenPartition b =
        RefinePartition(G, stats, partition, cpm, 0.05, rng2);

    CheckTrue("Test E theta reproducibility", a.community_of == b.community_of);
}

void TestSubsetSeedDeterminism()
{
    const std::uint64_t seed_a = MakeSubsetSeed(1234, 5, 9);
    const std::uint64_t seed_b = MakeSubsetSeed(1234, 5, 9);
    const std::uint64_t different_parent = MakeSubsetSeed(1234, 5, 10);
    const std::uint64_t different_level = MakeSubsetSeed(1234, 6, 9);

    CheckTrue("Test E seed deterministic", seed_a == seed_b);
    CheckTrue("Test E seed parent differs", seed_a != different_parent);
    CheckTrue("Test E seed level differs", seed_a != different_level);
}

void TestParallelRefinementDeterminismBaseline()
{
    const Graph G = MakeWeightedGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition partition =
        MakePartition(G, stats, {0, 0, 1, 1, 1});
    const CPMQualityFunction cpm(0.35);

    const LeidenPartition a =
        RefinePartition(G, stats, partition, cpm, 0.05, 2026, 3);
    const LeidenPartition b =
        RefinePartition(G, stats, partition, cpm, 0.05, 2026, 3);

    CheckSamePartition("Test E parallel deterministic baseline", a, b);
    CheckValidPartition("Test E parallel deterministic baseline", G, a);
    CheckRefinementOfPartition("Test E parallel deterministic baseline",
                               partition,
                               a);
}

void TestDifferentSeeds()
{
    const Graph G = MakeWeightedGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition partition =
        MakePartition(G, stats, {0, 0, 1, 1, 1});
    const ModularityQualityFunction modularity(1.0);

    CheckRefinementRun("Test F seed 1", G, stats, partition, modularity, 0.1, 1);
    CheckRefinementRun("Test F seed 2", G, stats, partition, modularity, 0.1, 2);
}

void TestAggregateNodeSizeCPM()
{
    Graph G = MakeGraph(3);
    add_undirected_edge(G, 0, 1, 1.0);
    add_undirected_edge(G, 1, 2, 0.5);
    const std::vector<double> node_size = {10.0, 1.0, 1.0};
    const LeidenGraphStats stats = BuildLeidenGraphStats(G, node_size);
    const LeidenPartition partition = MakePartition(G, stats, {0, 0, 1});
    const CPMQualityFunction cpm(0.2);
    std::vector<std::size_t> subset_mark(num_vertices(G), 0);
    std::size_t subset_generation = 0;
    MarkSubsetForTest(G, {0, 1}, subset_mark, subset_generation);

    CheckTrue("Test G CPM node_size threshold false",
              !IsNodeWellConnectedToSubset(G,
                                           stats,
                                           cpm,
                                           0,
                                           subset_mark,
                                           subset_generation,
                                           11.0));

    Graph H = MakeGraph(3);
    add_undirected_edge(H, 0, 1, 3.0);
    const LeidenGraphStats h_stats = BuildLeidenGraphStats(H, node_size);
    const LeidenPartition h_partition = MakePartition(H, h_stats, {0, 0, 1});
    std::mt19937 rng(77);
    const LeidenPartition refined =
        RefinePartition(H, h_stats, h_partition, cpm, 1.0e-6, rng);

    CheckValidPartition("Test G CPM aggregate node size", H, refined);
    CheckRefinementOfPartition("Test G CPM aggregate node size",
                               h_partition,
                               refined);
    CheckStatsAgainstRecompute("Test G CPM aggregate node size",
                               H,
                               h_stats,
                               refined);
    CheckTrue("Test G aggregate singleton can merge",
              refined.community_of[0] == refined.community_of[1]);

    (void)partition;
}

void TestModularityMass()
{
    Graph G = MakeGraph(3);
    add_undirected_edge(G, 0, 1, 1.0);
    add_undirected_edge(G, 1, 2, 9.0);
    const std::vector<double> huge_node_size = {1000.0, 1000.0, 1.0};
    const LeidenGraphStats default_stats = BuildLeidenGraphStats(G);
    const LeidenGraphStats huge_stats =
        BuildLeidenGraphStats(G, huge_node_size);
    const LeidenPartition default_partition =
        MakePartition(G, default_stats, {0, 0, 1});
    const LeidenPartition huge_partition =
        MakePartition(G, huge_stats, {0, 0, 1});
    const ModularityQualityFunction modularity(1.0);
    std::vector<std::size_t> subset_mark(num_vertices(G), 0);
    std::size_t subset_generation = 0;
    MarkSubsetForTest(G, {0, 1}, subset_mark, subset_generation);

    CheckTrue("Test H modularity strength threshold",
              IsNodeWellConnectedToSubset(G,
                                          huge_stats,
                                          modularity,
                                          0,
                                          subset_mark,
                                          subset_generation,
                                          11.0));

    std::mt19937 rng1(5);
    std::mt19937 rng2(5);
    const LeidenPartition refined_default =
        RefinePartition(G, default_stats, default_partition, modularity, 0.1, rng1);
    const LeidenPartition refined_huge =
        RefinePartition(G, huge_stats, huge_partition, modularity, 0.1, rng2);

    CheckTrue("Test H modularity ignores node_size for mass",
              refined_default.community_of == refined_huge.community_of);
}

void TestThetaValidation()
{
    const Graph G = MakeUnweightedGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition partition =
        MakePartition(G, stats, {0, 0, 0, 1, 1, 1});
    const CPMQualityFunction cpm(0.2);

    CheckThrowsInvalidArgument("Test I theta zero", G, stats, partition, cpm, 0.0);
    CheckThrowsInvalidArgument("Test I theta negative", G, stats, partition, cpm, -1.0);
    CheckRefinementRun("Test I tiny positive theta",
                       G,
                       stats,
                       partition,
                       cpm,
                       1.0e-12,
                       33);
}

Graph MakeStatsUpdateGraph()
{
    Graph G = MakeGraph(5);
    add_undirected_edge(G, 0, 0, 10.0);
    add_undirected_edge(G, 0, 1, 1.0);
    add_undirected_edge(G, 0, 2, 2.0);
    add_undirected_edge(G, 0, 3, 3.0);
    add_undirected_edge(G, 0, 4, 4.0);
    add_undirected_edge(G, 1, 2, 0.25);
    add_undirected_edge(G, 2, 3, 0.5);
    return G;
}

void CheckCachedStatsAgainstDirect(const std::string& test_name,
                                   const Graph& G,
                                   const LeidenGraphStats& stats,
                                   const LeidenPartition& refined,
                                   const QualityFunction& qf,
                                   const std::vector<Vertex>& subset,
                                   const std::vector<std::size_t>& subset_mark,
                                   std::size_t subset_generation,
                                   const RefinementCommunityStats& cached)
{
    const RefinementCommunityStats direct =
        DirectRefinementCommunityStats(G,
                                       stats,
                                       refined,
                                       qf,
                                       subset,
                                       subset_mark,
                                       subset_generation);
    CheckRefinementCommunityStatsEqual(test_name, direct, cached);
}

void TestRefinementCommunityStatsCache()
{
    const Graph G = MakeStatsUpdateGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    LeidenPartition refined = MakePartition(G, stats, {0, 0, 1, 2, 3});
    const CPMQualityFunction cpm(0.1);
    const std::vector<Vertex> subset = {0, 1, 2, 3};
    std::vector<std::size_t> subset_mark(num_vertices(G), 0);
    std::size_t subset_generation = 0;
    MarkSubsetForTest(G, subset, subset_mark, subset_generation);

    RefinementCommunityStats cached =
        BuildRefinementCommunityStats(G,
                                      stats,
                                      refined,
                                      cpm,
                                      subset,
                                      subset_mark,
                                      subset_generation);
    CheckCachedStatsAgainstDirect("Test J initial cached stats",
                                  G,
                                  stats,
                                  refined,
                                  cpm,
                                  subset,
                                  subset_mark,
                                  subset_generation,
                                  cached);

    UpdateRefinementCommunityStatsForMove(G,
                                          stats,
                                          refined,
                                          cpm,
                                          subset_mark,
                                          subset_generation,
                                          0,
                                          1,
                                          cached);
    MoveNodeToCommunity(G, stats, refined, 0, 1);
    CheckCachedStatsAgainstDirect("Test J single move cached stats",
                                  G,
                                  stats,
                                  refined,
                                  cpm,
                                  subset,
                                  subset_mark,
                                  subset_generation,
                                  cached);

    UpdateRefinementCommunityStatsForMove(G,
                                          stats,
                                          refined,
                                          cpm,
                                          subset_mark,
                                          subset_generation,
                                          1,
                                          1,
                                          cached);
    MoveNodeToCommunity(G, stats, refined, 1, 1);
    CheckCachedStatsAgainstDirect("Test J source empty cached stats",
                                  G,
                                  stats,
                                  refined,
                                  cpm,
                                  subset,
                                  subset_mark,
                                  subset_generation,
                                  cached);

    UpdateRefinementCommunityStatsForMove(G,
                                          stats,
                                          refined,
                                          cpm,
                                          subset_mark,
                                          subset_generation,
                                          3,
                                          1,
                                          cached);
    MoveNodeToCommunity(G, stats, refined, 3, 1);
    CheckCachedStatsAgainstDirect("Test J multiple moves cached stats",
                                  G,
                                  stats,
                                  refined,
                                  cpm,
                                  subset,
                                  subset_mark,
                                  subset_generation,
                                  cached);
}

void TestLocalStatsConstructionValues()
{
    const Graph G = MakeStatsUpdateGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition refined = MakePartition(G, stats, {0, 0, 1, 2, 3});
    const CPMQualityFunction cpm(0.1);
    const std::vector<Vertex> subset = {0, 1, 2, 3};
    std::vector<std::size_t> subset_mark(num_vertices(G), 0);
    std::size_t subset_generation = 0;
    MarkSubsetForTest(G, subset, subset_mark, subset_generation);

    const RefinementCommunityStats community_stats =
        BuildRefinementCommunityStats(G,
                                      stats,
                                      refined,
                                      cpm,
                                      subset,
                                      subset_mark,
                                      subset_generation);

    CheckTrue("Test K local stats entry count",
              community_stats.entries.size() == 3U);
    CheckTrue("Test K local stats active order",
              community_stats.active_communities ==
                  std::vector<Community>({0, 1, 2}));
    CheckEntryNear("Test K community 0", community_stats, 0, 2, 2.0, 5.25);
    CheckEntryNear("Test K community 1", community_stats, 1, 1, 1.0, 2.75);
    CheckEntryNear("Test K community 2", community_stats, 2, 1, 1.0, 3.5);
}

void TestLocalStatsMoveUpdateValues()
{
    const Graph G = MakeStatsUpdateGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    LeidenPartition refined = MakePartition(G, stats, {0, 0, 1, 2, 3});
    const CPMQualityFunction cpm(0.1);
    const std::vector<Vertex> subset = {0, 1, 2, 3};
    std::vector<std::size_t> subset_mark(num_vertices(G), 0);
    std::size_t subset_generation = 0;
    MarkSubsetForTest(G, subset, subset_mark, subset_generation);

    RefinementCommunityStats community_stats =
        BuildRefinementCommunityStats(G,
                                      stats,
                                      refined,
                                      cpm,
                                      subset,
                                      subset_mark,
                                      subset_generation);

    UpdateRefinementCommunityStatsForMove(G,
                                          stats,
                                          refined,
                                          cpm,
                                          subset_mark,
                                          subset_generation,
                                          0,
                                          1,
                                          community_stats);
    MoveNodeToCommunity(G, stats, refined, 0, 1);

    CheckCachedStatsAgainstDirect("Test L local stats after move",
                                  G,
                                  stats,
                                  refined,
                                  cpm,
                                  subset,
                                  subset_mark,
                                  subset_generation,
                                  community_stats);
    CheckTrue("Test L local stats entry count",
              community_stats.entries.size() == 3U);
    CheckTrue("Test L local stats active order",
              community_stats.active_communities ==
                  std::vector<Community>({0, 1, 2}));
    CheckEntryNear("Test L community 0", community_stats, 0, 1, 1.0, 1.25);
    CheckEntryNear("Test L community 1", community_stats, 1, 2, 2.0, 4.75);
    CheckEntryNear("Test L community 2", community_stats, 2, 1, 1.0, 3.5);
}

void TestLocalStatsSparseCommunityIds()
{
    Graph G = MakeGraph(2);
    add_undirected_edge(G, 0, 1, 7.0);
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    LeidenPartition refined;
    refined.community_of = {3, 1000000};
    const CPMQualityFunction cpm(0.1);
    const std::vector<Vertex> subset = {0, 1};
    std::vector<std::size_t> subset_mark(num_vertices(G), 0);
    std::size_t subset_generation = 0;
    MarkSubsetForTest(G, subset, subset_mark, subset_generation);

    const RefinementCommunityStats community_stats =
        BuildRefinementCommunityStats(G,
                                      stats,
                                      refined,
                                      cpm,
                                      subset,
                                      subset_mark,
                                      subset_generation);

    CheckTrue("Test M sparse local entry count",
              community_stats.entries.size() == 2U);
    CheckTrue("Test M sparse active order",
              community_stats.active_communities ==
                  std::vector<Community>({3, 1000000}));
    CheckEntryNear("Test M community 3", community_stats, 3, 1, 1.0, 7.0);
    CheckEntryNear("Test M community 1000000",
                   community_stats,
                   1000000,
                   1,
                   1.0,
                   7.0);
}

void TestSubsetMarkerMembershipEquivalence()
{
    const Graph G = MakeGraph(10);
    const std::vector<Vertex> subset = {1, 4, 7};
    std::vector<std::size_t> subset_mark(num_vertices(G), 0);
    std::size_t subset_generation = 0;
    MarkSubsetForTest(G, subset, subset_mark, subset_generation);

    for (Vertex v = 0; v < num_vertices(G); ++v) {
        const bool expected = (v == 1 || v == 4 || v == 7);
        CheckTrue("Test N marker membership",
                  IsInSubsetForTest(subset_mark,
                                    subset_generation,
                                    v) == expected);
    }
}

void TestSubsetMarkerConsecutiveSubsets()
{
    const Graph G = MakeGraph(8);
    std::vector<std::size_t> subset_mark(num_vertices(G), 0);
    std::size_t subset_generation = 0;

    MarkSubsetForTest(G, {1, 4}, subset_mark, subset_generation);
    CheckTrue("Test O generation 1 member 1",
              IsInSubsetForTest(subset_mark, subset_generation, 1));
    CheckTrue("Test O generation 1 member 4",
              IsInSubsetForTest(subset_mark, subset_generation, 4));

    MarkSubsetForTest(G, {2, 5}, subset_mark, subset_generation);
    CheckTrue("Test O generation 2 member 2",
              IsInSubsetForTest(subset_mark, subset_generation, 2));
    CheckTrue("Test O generation 2 member 5",
              IsInSubsetForTest(subset_mark, subset_generation, 5));
    CheckTrue("Test O generation 2 stale 1",
              !IsInSubsetForTest(subset_mark, subset_generation, 1));
    CheckTrue("Test O generation 2 stale 4",
              !IsInSubsetForTest(subset_mark, subset_generation, 4));
}

void TestExactSparseTargetCorners()
{
    Graph G = MakeGraph(4);
    add_undirected_edge(G, 1, 2, 1.25);
    add_undirected_edge(G, 1, 2, 2.75); // parallel weighted edge
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition singleton = MakeSingletonPartition(G, stats);
    NeighborCommunityScratch scratch;
    scratch.weights.assign(singleton.community_size.size(), 0.0);
    scratch.marks.assign(singleton.community_size.size(), 0);
    BuildNeighborCommunityWeights(G, singleton, 1, scratch);

    RefinementCommunityStats community_stats;
    for (Community c = 0; c < 4; ++c) {
        RefinementCommunityEntry entry;
        entry.member_count = 1;
        entry.mass = stats.node_strength[static_cast<std::size_t>(c)];
        community_stats.entries.emplace(c, entry);
        community_stats.active_communities.push_back(c);
        if (entry.mass <= 0.0) {
            community_stats.nonpositive_mass_active_communities.push_back(c);
        }
    }

    const ModularityQualityFunction modularity(1.0);
    std::size_t exceptions = 0;
    const std::vector<Community> targets =
        BuildExactSparseRefinementTargets(stats,
                                          modularity,
                                          community_stats,
                                          scratch,
                                          1,
                                          1,
                                          &exceptions);
    CheckTrue("Test P modularity zero-strength exception",
              targets == std::vector<Community>({0, 1, 2, 3}));
    CheckTrue("Test P modularity exception count", exceptions == 2);
    CheckNear("Test P modularity nonneighbor zero delta",
              0.0,
              modularity.deltaMoveFromWeights(
                  stats, singleton, 1, 3, 0.0, 0.0));

    NeighborCommunityScratch isolated_scratch;
    isolated_scratch.weights.assign(singleton.community_size.size(), 0.0);
    isolated_scratch.marks.assign(singleton.community_size.size(), 0);
    BuildNeighborCommunityWeights(G, singleton, 0, isolated_scratch);
    const std::vector<Community> isolated_targets =
        BuildExactSparseRefinementTargets(stats,
                                          modularity,
                                          community_stats,
                                          isolated_scratch,
                                          0,
                                          0);
    CheckTrue("Test P isolated source full fallback",
              isolated_targets == community_stats.active_communities);

    const LeidenGraphStats zero_size_stats =
        BuildLeidenGraphStats(G, {1.0, 1.0, 1.0, 0.0});
    RefinementCommunityStats cpm_stats = community_stats;
    cpm_stats.nonpositive_mass_active_communities = {3};
    for (Community c = 0; c < 4; ++c) {
        cpm_stats.entries[c].mass =
            zero_size_stats.node_size[static_cast<std::size_t>(c)];
    }
    const CPMQualityFunction cpm(0.5);
    exceptions = 0;
    const std::vector<Community> cpm_targets =
        BuildExactSparseRefinementTargets(zero_size_stats,
                                          cpm,
                                          cpm_stats,
                                          scratch,
                                          1,
                                          1,
                                          &exceptions);
    CheckTrue("Test P CPM zero-size exception",
              cpm_targets == std::vector<Community>({1, 2, 3}));
    CheckTrue("Test P CPM exception count", exceptions == 1);

    const UnsupportedSparseCPM unsupported(0.5);
    const std::vector<Community> fallback_targets =
        BuildExactSparseRefinementTargets(zero_size_stats,
                                          unsupported,
                                          cpm_stats,
                                          scratch,
                                          1,
                                          1);
    CheckTrue("Test P unknown quality full fallback",
              fallback_targets == cpm_stats.active_communities);
}

} // namespace

int main()
{
    TestRefinementOfMoveNodesFastPartition();
    TestNoIllegalCrossMerge();
    TestWeighted();
    TestSelfLoop();
    TestThetaReproducibility();
    TestSubsetSeedDeterminism();
    TestParallelRefinementDeterminismBaseline();
    TestDifferentSeeds();
    TestAggregateNodeSizeCPM();
    TestModularityMass();
    TestThetaValidation();
    TestRefinementCommunityStatsCache();
    TestLocalStatsConstructionValues();
    TestLocalStatsMoveUpdateValues();
    TestLocalStatsSparseCommunityIds();
    TestSubsetMarkerMembershipEquivalence();
    TestSubsetMarkerConsecutiveSubsets();
    TestExactSparseTargetCorners();

    std::cout << "All refinement tests passed.\n";
    return EXIT_SUCCESS;
}
