#include "Leiden.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
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

void CheckEqual(const std::string& test_name,
                std::size_t expected,
                std::size_t actual)
{
    if (expected != actual) {
        Fail(test_name,
             std::to_string(expected),
             std::to_string(actual),
             std::to_string(static_cast<long long>(actual) -
                            static_cast<long long>(expected)));
    }
}

void CheckEqual(const std::string& test_name, int expected, int actual)
{
    if (expected != actual) {
        Fail(test_name,
             std::to_string(expected),
             std::to_string(actual),
             std::to_string(actual - expected));
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

double Sum(const std::vector<double>& values)
{
    return std::accumulate(values.begin(), values.end(), 0.0);
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

void CheckValidOriginalPartition(const std::string& test_name,
                                 const Graph& G,
                                 const LeidenPartition& partition)
{
    CheckEqual(test_name + " assignment size",
               static_cast<std::size_t>(num_vertices(G)),
               partition.community_of.size());
    for (Vertex v = 0; v < num_vertices(G); ++v) {
        const Community community = partition.community_of[v];
        CheckTrue(test_name + " assigned", community >= 0);
        CheckTrue(test_name + " in range",
                  static_cast<std::size_t>(community) <
                      partition.community_size.size());
    }
}

void CheckQualityNondecreasing(const std::string& test_name,
                               const Graph& G,
                               const LeidenGraphStats& stats,
                               const LeidenPartition& result_partition,
                               const QualityFunction& qf)
{
    const LeidenPartition singleton = MakeSingletonPartition(G, stats);
    const double before = qf.quality(G, stats, singleton);
    const double after = qf.quality(G, stats, result_partition);
    CheckTrue(test_name, after + kTolerance >= before);
}

Graph MakeNoEdgeGraph(int n)
{
    return MakeGraph(n);
}

Graph MakeWeightedGraph()
{
    Graph G = MakeGraph(6);
    add_undirected_edge(G, 0, 1, 2.0);
    add_undirected_edge(G, 0, 2, 1.0);
    add_undirected_edge(G, 1, 2, 2.5);
    add_undirected_edge(G, 3, 4, 3.0);
    add_undirected_edge(G, 4, 5, 2.0);
    add_undirected_edge(G, 3, 5, 2.5);
    add_undirected_edge(G, 2, 3, 0.5);
    return G;
}

Graph MakeSelfLoopGraph()
{
    Graph G = MakeGraph(5);
    add_undirected_edge(G, 0, 0, 2.0);
    add_undirected_edge(G, 0, 1, 3.0);
    add_undirected_edge(G, 1, 2, 1.5);
    add_undirected_edge(G, 2, 2, 4.0);
    add_undirected_edge(G, 2, 3, 2.0);
    add_undirected_edge(G, 3, 4, 3.5);
    add_undirected_edge(G, 4, 4, 1.0);
    return G;
}

Graph MakeHierarchicalGraph()
{
    Graph G = MakeGraph(4);
    add_undirected_edge(G, 0, 1, 10.0);
    add_undirected_edge(G, 2, 3, 10.0);
    add_undirected_edge(G, 0, 2, 1.5);
    add_undirected_edge(G, 0, 3, 1.5);
    add_undirected_edge(G, 1, 2, 1.5);
    add_undirected_edge(G, 1, 3, 1.5);
    return G;
}

void TestSingleVertex()
{
    const Graph G = MakeGraph(1);
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const CPMQualityFunction cpm(0.5);
    const LeidenResult result = Leiden(G, stats, cpm, LeidenOptions{0.01, 1, 0});

    CheckEqual("Test A levels", 1U, result.num_levels);
    CheckEqual("Test A moves", 0U, result.total_moves);
    CheckEqual("Test A assignment", 0, result.partition.community_of[0]);
    CheckStatsAgainstRecompute("Test A stats", G, stats, result.partition);
}

void TestOneLevelConvergence()
{
    const Graph G = MakeNoEdgeGraph(4);
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const CPMQualityFunction cpm(0.5);
    const LeidenResult result = Leiden(G, stats, cpm, LeidenOptions{0.01, 2, 0});

    CheckEqual("Test B levels", 1U, result.num_levels);
    CheckEqual("Test B moves", 0U, result.total_moves);
    CheckValidOriginalPartition("Test B valid", G, result.partition);
}

void TestAtLeastOneAggregation()
{
    const Graph G = MakeHierarchicalGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const CPMQualityFunction cpm(1.0);
    const LeidenResult result = Leiden(G, stats, cpm, LeidenOptions{0.01, 3, 0});

    CheckTrue("Test C levels", result.num_levels >= 2);
    CheckTrue("Test C moves", result.total_moves > 0);
    CheckValidOriginalPartition("Test C valid", G, result.partition);
    CheckQualityNondecreasing("Test C quality", G, stats, result.partition, cpm);
}

void TestCPMFullLeiden()
{
    const Graph G = MakeWeightedGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const CPMQualityFunction cpm(0.3);
    const LeidenOptions options{0.01, 123, 0};
    const LeidenResult a = Leiden(G, stats, cpm, options);
    const LeidenResult b = Leiden(G, stats, cpm, options);

    CheckValidOriginalPartition("Test D valid", G, a.partition);
    CheckQualityNondecreasing("Test D quality", G, stats, a.partition, cpm);
    CheckTrue("Test D reproducibility assignment",
              a.partition.community_of == b.partition.community_of);
    CheckEqual("Test D reproducibility levels", a.num_levels, b.num_levels);
    CheckEqual("Test D reproducibility moves", a.total_moves, b.total_moves);
}

void TestModularityFullLeiden()
{
    const Graph G = MakeWeightedGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const ModularityQualityFunction modularity(1.0);
    const LeidenOptions options{0.01, 456, 0};
    const LeidenResult a = Leiden(G, stats, modularity, options);
    const LeidenResult b = Leiden(G, stats, modularity, options);

    CheckValidOriginalPartition("Test E valid", G, a.partition);
    CheckQualityNondecreasing("Test E quality", G, stats, a.partition, modularity);
    CheckTrue("Test E reproducibility assignment",
              a.partition.community_of == b.partition.community_of);
    CheckEqual("Test E reproducibility levels", a.num_levels, b.num_levels);
    CheckEqual("Test E reproducibility moves", a.total_moves, b.total_moves);
}

void TestSelfLoop()
{
    const Graph G = MakeSelfLoopGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const CPMQualityFunction cpm(0.4);
    const ModularityQualityFunction modularity(1.2);
    const LeidenResult cpm_result = Leiden(G, stats, cpm, LeidenOptions{0.01, 7, 0});
    const LeidenResult modularity_result =
        Leiden(G, stats, modularity, LeidenOptions{0.01, 7, 0});

    CheckStatsAgainstRecompute("Test F CPM stats", G, stats, cpm_result.partition);
    CheckStatsAgainstRecompute("Test F modularity stats",
                               G,
                               stats,
                               modularity_result.partition);
    CheckQualityNondecreasing("Test F CPM quality", G, stats, cpm_result.partition, cpm);
    CheckQualityNondecreasing("Test F modularity quality",
                              G,
                              stats,
                              modularity_result.partition,
                              modularity);
}

void TestAggregateNodeSize()
{
    const Graph G = MakeHierarchicalGraph();
    const std::vector<double> node_size = {2.0, 3.0, 5.0, 7.0};
    const LeidenGraphStats stats = BuildLeidenGraphStats(G, node_size);
    const CPMQualityFunction cpm(1.0);
    const LeidenResult result = Leiden(G, stats, cpm, LeidenOptions{0.01, 8, 0});

    CheckValidOriginalPartition("Test G valid", G, result.partition);
    CheckNear("Test G node size sum", Sum(stats.node_size),
              Sum(result.partition.community_size));
    CheckStatsAgainstRecompute("Test G stats", G, stats, result.partition);
}

void TestSameSeedReproducibility()
{
    const Graph G = MakeHierarchicalGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const CPMQualityFunction cpm(1.0);
    const LeidenOptions options{0.01, 2026, 0};
    const LeidenResult a = Leiden(G, stats, cpm, options);
    const LeidenResult b = Leiden(G, stats, cpm, options);

    CheckTrue("Test H assignment", a.partition.community_of == b.partition.community_of);
    CheckEqual("Test H levels", a.num_levels, b.num_levels);
    CheckEqual("Test H moves", a.total_moves, b.total_moves);
    CheckNear("Test H quality",
              cpm.quality(G, stats, a.partition),
              cpm.quality(G, stats, b.partition));
}

void TestDifferentSeeds()
{
    const Graph G = MakeWeightedGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const ModularityQualityFunction modularity(1.0);
    const LeidenResult a = Leiden(G, stats, modularity, LeidenOptions{0.01, 1, 0});
    const LeidenResult b = Leiden(G, stats, modularity, LeidenOptions{0.01, 2, 0});

    CheckValidOriginalPartition("Test I seed 1", G, a.partition);
    CheckValidOriginalPartition("Test I seed 2", G, b.partition);
    CheckQualityNondecreasing("Test I quality 1", G, stats, a.partition, modularity);
    CheckQualityNondecreasing("Test I quality 2", G, stats, b.partition, modularity);
}

void TestFlattenMapping()
{
    const Graph G = MakeHierarchicalGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const CPMQualityFunction cpm(1.0);
    const LeidenResult result = Leiden(G, stats, cpm, LeidenOptions{0.01, 11, 0});

    CheckTrue("Test J multi-level", result.num_levels >= 3);
    CheckTrue("Test J pair 0-1 same",
              result.partition.community_of[0] == result.partition.community_of[1]);
    CheckTrue("Test J pair 2-3 same",
              result.partition.community_of[2] == result.partition.community_of[3]);
}

void TestFinalStatistics()
{
    const Graph G = MakeWeightedGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const CPMQualityFunction cpm(0.3);
    const LeidenResult result = Leiden(G, stats, cpm, LeidenOptions{0.01, 12, 0});
    CheckStatsAgainstRecompute("Test K final stats", G, stats, result.partition);
}

void TestLevelQualityPreservation()
{
    Graph current_graph = MakeHierarchicalGraph();
    LeidenGraphStats current_stats = BuildLeidenGraphStats(current_graph);
    const CPMQualityFunction cpm(1.0);
    LeidenPartition current_partition =
        MakeSingletonPartition(current_graph, current_stats);
    std::mt19937 rng(13);

    const MoveNodesFastResult moved =
        MoveNodesFast(current_graph, current_stats, current_partition, cpm, rng);
    const LeidenPartition moved_partition = moved.partition;
    const double before_aggregation =
        cpm.quality(current_graph, current_stats, moved_partition);
    const LeidenPartition refined =
        RefinePartition(current_graph,
                        current_stats,
                        moved_partition,
                        cpm,
                        0.01,
                        rng);
    const AggregateGraphResult aggregate =
        AggregateGraph(current_graph, current_stats, refined);
    const LeidenPartition coarse_partition =
        BuildCoarsePartition(aggregate, moved_partition, refined);
    const double after_aggregation =
        cpm.quality(aggregate.graph, aggregate.stats, coarse_partition);

    CheckNear("Test L level quality preservation",
              before_aggregation,
              after_aggregation);
}

void TestMaxLevels()
{
    const Graph G = MakeHierarchicalGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const CPMQualityFunction cpm(1.0);
    const LeidenResult result = Leiden(G, stats, cpm, LeidenOptions{0.01, 14, 1});

    CheckEqual("Test M max levels", 1U, result.num_levels);
    CheckValidOriginalPartition("Test M valid", G, result.partition);
    CheckQualityNondecreasing("Test M quality", G, stats, result.partition, cpm);
}

void TestEmptyGraph()
{
    const Graph G = MakeGraph(0);
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const CPMQualityFunction cpm(0.5);
    const LeidenResult result = Leiden(G, stats, cpm, LeidenOptions{0.01, 15, 0});

    CheckEqual("Test N levels", 0U, result.num_levels);
    CheckEqual("Test N moves", 0U, result.total_moves);
    CheckTrue("Test N empty assignment", result.partition.community_of.empty());
}

void TestDebugDoesNotChangeResult()
{
    const Graph G = MakeWeightedGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const ModularityQualityFunction modularity(1.0);

    LeidenOptions normal{0.01, 2027, 0};
    normal.debug = false;
    normal.debug_interval = 2;

    std::ostringstream normal_stderr;
    std::streambuf* old_stderr = std::cerr.rdbuf(normal_stderr.rdbuf());
    const LeidenResult normal_result = Leiden(G, stats, modularity, normal);
    std::cerr.rdbuf(old_stderr);

    CheckTrue("Test O normal debug output empty",
              normal_stderr.str().empty());

    LeidenOptions debug = normal;
    debug.debug = true;

    std::ostringstream debug_stderr;
    old_stderr = std::cerr.rdbuf(debug_stderr.rdbuf());
    const LeidenResult debug_result = Leiden(G, stats, modularity, debug);
    std::cerr.rdbuf(old_stderr);

    CheckTrue("Test O assignment",
              normal_result.partition.community_of ==
                  debug_result.partition.community_of);
    CheckEqual("Test O levels",
               normal_result.num_levels,
               debug_result.num_levels);
    CheckEqual("Test O moves",
               normal_result.total_moves,
               debug_result.total_moves);
    CheckNear("Test O quality",
              modularity.quality(G, stats, normal_result.partition),
              modularity.quality(G, stats, debug_result.partition));
    CheckTrue("Test O debug output present",
              debug_stderr.str().find("[Leiden] level 1 start") !=
                  std::string::npos);
}

} // namespace

int main()
{
    TestSingleVertex();
    TestOneLevelConvergence();
    TestAtLeastOneAggregation();
    TestCPMFullLeiden();
    TestModularityFullLeiden();
    TestSelfLoop();
    TestAggregateNodeSize();
    TestSameSeedReproducibility();
    TestDifferentSeeds();
    TestFlattenMapping();
    TestFinalStatistics();
    TestLevelQualityPreservation();
    TestMaxLevels();
    TestEmptyGraph();
    TestDebugDoesNotChangeResult();

    std::cout << "All Leiden tests passed.\n";
    return EXIT_SUCCESS;
}
