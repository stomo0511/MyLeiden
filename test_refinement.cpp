#include "Leiden.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kTolerance = 1.0e-10;

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
    const std::vector<bool> subset = {true, true, false, false};

    CheckNear("Test D self-loop excluded from E(v,S-v)",
              1.5,
              EdgeWeightFromNodeToSubset(G, 0, subset));
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
    const std::vector<bool> subset = {true, true, false};

    CheckTrue("Test G CPM node_size threshold false",
              !IsNodeWellConnectedToSubset(G, stats, cpm, 0, subset, 11.0));

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
    const std::vector<bool> subset = {true, true, false};

    CheckTrue("Test H modularity strength threshold",
              IsNodeWellConnectedToSubset(G,
                                          huge_stats,
                                          modularity,
                                          0,
                                          subset,
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

} // namespace

int main()
{
    TestRefinementOfMoveNodesFastPartition();
    TestNoIllegalCrossMerge();
    TestWeighted();
    TestSelfLoop();
    TestThetaReproducibility();
    TestDifferentSeeds();
    TestAggregateNodeSizeCPM();
    TestModularityMass();
    TestThetaValidation();

    std::cout << "All refinement tests passed.\n";
    return EXIT_SUCCESS;
}
