#include "Leiden.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
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
    for (Community c = 0;
         c < static_cast<Community>(partition.community_size.size());
         ++c) {
        if (partition.community_size[c] == 0.0) {
            return c;
        }
    }
    return static_cast<Community>(partition.community_size.size());
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
        CheckTrue(test_name + " non-negative community", c >= 0);
        CheckTrue(test_name + " community in statistic range",
                  static_cast<std::size_t>(c) < nc);
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

} // namespace

int main()
{
    TestCPMUnweighted();
    TestModularityUnweighted();
    TestWeighted();
    TestSelfLoop();
    TestReproducibility();
    TestDifferentSeeds();

    std::cout << "All MoveNodesFast tests passed.\n";
    return EXIT_SUCCESS;
}
