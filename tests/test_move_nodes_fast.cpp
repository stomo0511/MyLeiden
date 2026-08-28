#include "Leiden.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

constexpr double kTolerance = 1.0e-10;
constexpr double kLocalOptimalityTolerance = 1.0e-12;

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

void CheckTrue(const std::string& test_name, bool condition)
{
    if (!condition) {
        Fail(test_name, "true", "false", "condition is false");
    }
}

double LookupWeight(const std::unordered_map<Community, double>& weights,
                    Community community)
{
    const auto it = weights.find(community);
    return (it == weights.end()) ? 0.0 : it->second;
}

Community EmptyCommunityForTest(const LeidenPartition& partition)
{
    return EmptyCommunityForMove(partition);
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
              partition.internal_edge_weight.size() == nc &&
              partition.community_is_empty.size() == nc);
    CheckTrue(test_name + " smallest empty range",
              partition.smallest_empty_community >= 0 &&
              static_cast<std::size_t>(partition.smallest_empty_community) <= nc);

    for (Vertex v = 0; v < num_vertices(G); ++v) {
        const Community c = partition.community_of[v];
        CheckTrue(test_name + " non-negative community", c >= 0);
        CheckTrue(test_name + " community in statistic range",
                  static_cast<std::size_t>(c) < nc);
    }

    for (Community c = 0; c < static_cast<Community>(nc); ++c) {
        CheckTrue(test_name + " empty-community index",
                  (partition.community_size[c] == 0.0) ==
                      (partition.community_is_empty[c] != 0));
    }
}

void CheckPartitionStatsAgainstRecompute(const std::string& test_name,
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

void CheckLocalOptimality(const std::string& test_name,
                          const Graph& G,
                          const LeidenGraphStats& stats,
                          const LeidenPartition& partition,
                          const QualityFunction& qf)
{
    // A tiny tolerance absorbs roundoff near zero; MoveNodesFast itself uses
    // the paper's strict best_delta > 0.0 rule.
    for (Vertex v = 0; v < num_vertices(G); ++v) {
        const Community source = partition.community_of[v];
        const auto weights = BuildNeighborCommunityWeights(G, partition, v);
        std::vector<Community> candidates;
        candidates.reserve(weights.size() + 1);
        for (const auto& item : weights) {
            candidates.push_back(item.first);
        }
        candidates.push_back(EmptyCommunityForTest(partition));
        std::sort(candidates.begin(), candidates.end());
        candidates.erase(std::unique(candidates.begin(), candidates.end()),
                         candidates.end());

        const double weight_to_source = LookupWeight(weights, source);
        for (Community target : candidates) {
            const double delta =
                qf.deltaMoveFromWeights(stats,
                                        partition,
                                        v,
                                        target,
                                        weight_to_source,
                                        LookupWeight(weights, target));
            if (delta > kLocalOptimalityTolerance) {
                Fail(test_name + " local optimality",
                     "delta <= " + std::to_string(kLocalOptimalityTolerance),
                     std::to_string(delta),
                     "v=" + std::to_string(v) +
                         " target=" + std::to_string(target));
            }
        }
    }
}

MoveNodesFastResult RunAndCheck(const std::string& test_name,
                                const Graph& G,
                                const QualityFunction& qf,
                                unsigned int seed,
                                bool require_move)
{
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition initial = MakeSingletonPartition(G, stats);
    const double quality_before = qf.quality(G, stats, initial);
    std::mt19937 rng(seed);

    MoveNodesFastResult result =
        MoveNodesFast(G, stats, initial, qf, rng);
    const double quality_after = qf.quality(G, stats, result.partition);

    CheckValidPartition(test_name, G, result.partition);
    CheckTrue(test_name + " quality nondecreasing",
              quality_after + kTolerance >= quality_before);
    if (require_move) {
        CheckTrue(test_name + " at least one move", result.num_moves > 0);
    }
    CheckPartitionStatsAgainstRecompute(test_name, G, stats, result.partition);
    CheckLocalOptimality(test_name, G, stats, result.partition, qf);
    return result;
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

void TestCPMUnweighted()
{
    const CPMQualityFunction cpm(0.4);
    RunAndCheck("Test A CPM unweighted", MakeUnweightedGraph(), cpm, 1234, true);
}

void TestModularityUnweighted()
{
    const ModularityQualityFunction modularity(1.0);
    RunAndCheck("Test B modularity unweighted",
                MakeUnweightedGraph(),
                modularity,
                1234,
                true);
}

void TestWeighted()
{
    const Graph G = MakeWeightedGraph();
    const CPMQualityFunction cpm(0.35);
    const ModularityQualityFunction modularity(0.8);
    RunAndCheck("Test C CPM weighted", G, cpm, 777, true);
    RunAndCheck("Test C modularity weighted", G, modularity, 777, true);
}

void TestSelfLoop()
{
    const Graph G = MakeSelfLoopGraph();
    const CPMQualityFunction cpm(0.2);
    const ModularityQualityFunction modularity(1.2);
    RunAndCheck("Test D CPM self-loop", G, cpm, 99, false);
    RunAndCheck("Test D modularity self-loop", G, modularity, 99, false);
}

void TestReproducibility()
{
    const Graph G = MakeWeightedGraph();
    const ModularityQualityFunction modularity(1.0);
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition initial = MakeSingletonPartition(G, stats);

    std::mt19937 rng1(2026);
    std::mt19937 rng2(2026);
    const MoveNodesFastResult a =
        MoveNodesFast(G, stats, initial, modularity, rng1);
    const MoveNodesFastResult b =
        MoveNodesFast(G, stats, initial, modularity, rng2);

    CheckTrue("Test E reproducibility assignment",
              a.partition.community_of == b.partition.community_of);
    CheckTrue("Test E reproducibility moves", a.num_moves == b.num_moves);
    CheckNear("Test E reproducibility quality",
              modularity.quality(G, stats, a.partition),
              modularity.quality(G, stats, b.partition));
}

void TestDifferentSeeds()
{
    const Graph G = MakeWeightedGraph();
    const CPMQualityFunction cpm(0.35);
    RunAndCheck("Test F CPM seed 1", G, cpm, 1, false);
    RunAndCheck("Test F CPM seed 2", G, cpm, 2, false);

    const ModularityQualityFunction modularity(1.0);
    RunAndCheck("Test F modularity seed 1", G, modularity, 1, false);
    RunAndCheck("Test F modularity seed 2", G, modularity, 2, false);
}

void TestEmptyCommunityManagement()
{
    const Graph G = MakeGraph(8);
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);

    const LeidenPartition singleton = MakeSingletonPartition(G, stats);
    CheckTrue("Test G singleton flags",
              std::all_of(singleton.community_is_empty.begin(),
                          singleton.community_is_empty.end(),
                          [](unsigned char flag) { return flag == 0; }));
    CheckTrue("Test G singleton hint",
              singleton.smallest_empty_community == 8);
    CheckTrue("Test G singleton next community",
              EmptyCommunityForMove(singleton) == 8);

    LeidenPartition partition =
        MakePartition(G, stats, {0, 1, 2, 3, 4, 5, 6, 7});
    MoveNodeToCommunity(G, stats, partition, 3, 0);
    CheckTrue("Test G community becomes empty",
              partition.community_is_empty[3] != 0);
    CheckTrue("Test G newly empty updates hint",
              partition.smallest_empty_community <= 3);
    CheckTrue("Test G newly empty is smallest",
              EmptyCommunityForMove(partition) == 3);

    const Graph sparse_graph = MakeGraph(12);
    const LeidenGraphStats sparse_stats = BuildLeidenGraphStats(sparse_graph);
    LeidenPartition sparse = MakePartition(
        sparse_graph,
        sparse_stats,
        {0, 0, 1, 3, 4, 5, 6, 8, 9, 11, 11, 11});
    CheckTrue("Test G sparse empty communities",
              sparse.community_is_empty[2] != 0 &&
              sparse.community_is_empty[7] != 0 &&
              sparse.community_is_empty[10] != 0);
    CheckTrue("Test G sparse smallest hint",
              sparse.smallest_empty_community == 2);
    CheckTrue("Test G smallest sparse empty",
              EmptyCommunityForMove(sparse) == 2);
    MoveNodeToCommunity(sparse_graph, sparse_stats, sparse, 0, 2);
    CheckTrue("Test G reused community removed",
              sparse.community_is_empty[2] == 0);
    CheckTrue("Test G next sparse empty",
              EmptyCommunityForMove(sparse) == 7);

    Graph expansion_graph = MakeGraph(5);
    const LeidenGraphStats expansion_stats =
        BuildLeidenGraphStats(expansion_graph);
    LeidenPartition expanded =
        MakeSingletonPartition(expansion_graph, expansion_stats);
    MoveNodeToCommunity(expansion_graph, expansion_stats, expanded, 0, 7);
    CheckTrue("Test G expansion gaps empty",
              expanded.community_is_empty[0] != 0 &&
              expanded.community_is_empty[5] != 0 &&
              expanded.community_is_empty[6] != 0);
    CheckTrue("Test G expansion target active",
              expanded.community_is_empty[7] == 0);
    CheckTrue("Test G expansion smallest",
              EmptyCommunityForMove(expanded) == 0);

    const LeidenPartition copied = expanded;
    CheckTrue("Test G copy preserves empty communities",
              copied.community_is_empty == expanded.community_is_empty &&
              copied.smallest_empty_community ==
                  expanded.smallest_empty_community);
}

Community ReferenceEmptyCommunityForMove(const std::set<Community>& empty,
                                         std::size_t community_count)
{
    if (!empty.empty()) {
        return *empty.begin();
    }
    return static_cast<Community>(community_count);
}

void ReferenceEnsureCommunity(std::set<Community>& empty,
                              std::size_t& community_count,
                              Community community)
{
    const std::size_t need = static_cast<std::size_t>(community) + 1;
    for (std::size_t c = community_count; c < need; ++c) {
        empty.insert(static_cast<Community>(c));
    }
    community_count = std::max(community_count, need);
}

void TestEmptyCommunityReferenceEquivalence()
{
    const Graph G = MakeGraph(6);
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    LeidenPartition partition =
        MakePartition(G, stats, {0, 0, 1, 3, 5, 5});

    std::set<Community> reference_empty = {2, 4};
    std::size_t reference_count = 6;
    CheckTrue("Test H initial reference empty",
              EmptyCommunityForMove(partition) ==
                  ReferenceEmptyCommunityForMove(reference_empty,
                                                 reference_count));

    const std::vector<std::pair<Vertex, Community>> moves = {
        {0, 2},
        {2, 4},
        {3, 0},
        {1, 3},
        {4, 1},
        {5, 7},
        {0, 6}
    };

    std::vector<Community> reference_assignment = partition.community_of;
    std::vector<double> reference_size = partition.community_size;
    for (std::size_t i = 0; i < moves.size(); ++i) {
        const Vertex v = moves[i].first;
        const Community target = moves[i].second;
        const Community source = reference_assignment[v];

        ReferenceEnsureCommunity(reference_empty, reference_count, target);
        reference_size.resize(reference_count, 0.0);
        reference_size[source] -= stats.node_size[v];
        reference_assignment[v] = -1;
        if (reference_size[source] == 0.0) {
            reference_empty.insert(source);
        }
        reference_assignment[v] = target;
        reference_size[target] += stats.node_size[v];
        reference_empty.erase(target);

        MoveNodeToCommunity(G, stats, partition, v, target);

        const std::string label =
            "Test H step " + std::to_string(i);
        CheckTrue(label,
                  EmptyCommunityForMove(partition) ==
                      ReferenceEmptyCommunityForMove(reference_empty,
                                                     reference_count));
    }
}

} // namespace

int main()
{
    TestCPMUnweighted();
    TestModularityUnweighted();
    TestWeighted();
    TestSelfLoop();
    TestReproducibility();
    TestDifferentSeeds();
    TestEmptyCommunityManagement();
    TestEmptyCommunityReferenceEquivalence();

    std::cout << "All MoveNodesFast tests passed.\n";
    return EXIT_SUCCESS;
}
