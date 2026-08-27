#include "Leiden.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
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
                                const AggregateGraphResult& aggregate,
                                const LeidenPartition& coarse_partition)
{
    const LeidenPartition recomputed =
        MakePartition(aggregate.graph,
                      aggregate.stats,
                      coarse_partition.community_of);
    const std::size_t nc = std::max(coarse_partition.community_size.size(),
                                    recomputed.community_size.size());
    for (std::size_t c = 0; c < nc; ++c) {
        const double expected_size =
            (c < recomputed.community_size.size()) ? recomputed.community_size[c] : 0.0;
        const double actual_size =
            (c < coarse_partition.community_size.size()) ? coarse_partition.community_size[c] : 0.0;
        const double expected_strength =
            (c < recomputed.community_strength.size()) ? recomputed.community_strength[c] : 0.0;
        const double actual_strength =
            (c < coarse_partition.community_strength.size()) ? coarse_partition.community_strength[c] : 0.0;
        const double expected_internal =
            (c < recomputed.internal_edge_weight.size()) ? recomputed.internal_edge_weight[c] : 0.0;
        const double actual_internal =
            (c < coarse_partition.internal_edge_weight.size()) ? coarse_partition.internal_edge_weight[c] : 0.0;

        CheckNear(test_name + " community_size", expected_size, actual_size);
        CheckNear(test_name + " community_strength", expected_strength, actual_strength);
        CheckNear(test_name + " internal_edge_weight", expected_internal, actual_internal);
    }
}

void CheckQualityPreservation(const std::string& test_name,
                              const Graph& G,
                              const LeidenGraphStats& stats,
                              const LeidenPartition& partition,
                              const LeidenPartition& refined,
                              const QualityFunction& quality_function)
{
    const AggregateGraphResult aggregate = AggregateGraph(G, stats, refined);
    const LeidenPartition coarse_partition =
        BuildCoarsePartition(aggregate, partition, refined);

    // Prefined defines the aggregate vertices, while P defines their coarse
    // communities. Since Prefined is a refinement of P, grouping aggregate
    // vertices by P represents exactly P on the original graph.
    CheckNear(test_name,
              quality_function.quality(G, stats, partition),
              quality_function.quality(aggregate.graph,
                                       aggregate.stats,
                                       coarse_partition));
}

Graph MakeSimpleGraph()
{
    Graph G = MakeGraph(6);
    add_undirected_edge(G, 0, 1, 1.0);
    add_undirected_edge(G, 1, 2, 2.0);
    add_undirected_edge(G, 2, 3, 3.0);
    add_undirected_edge(G, 3, 4, 4.0);
    add_undirected_edge(G, 4, 5, 5.0);
    add_undirected_edge(G, 0, 5, 6.0);
    return G;
}

Graph MakeWeightedSelfLoopGraph()
{
    Graph G = MakeGraph(6);
    add_undirected_edge(G, 0, 0, 2.0);
    add_undirected_edge(G, 0, 1, 1.5);
    add_undirected_edge(G, 1, 2, 2.5);
    add_undirected_edge(G, 2, 3, 3.5);
    add_undirected_edge(G, 3, 4, 4.5);
    add_undirected_edge(G, 4, 5, 5.5);
    add_undirected_edge(G, 5, 5, 6.5);
    add_undirected_edge(G, 1, 4, 7.5);
    return G;
}

void TestBasicCoarsePartition()
{
    const Graph G = MakeSimpleGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition partition =
        MakePartition(G, stats, {0, 0, 0, 1, 1, 1});
    const LeidenPartition refined =
        MakePartition(G, stats, {2, 2, 5, 7, 7, 9});
    const AggregateGraphResult aggregate = AggregateGraph(G, stats, refined);
    const LeidenPartition coarse =
        BuildCoarsePartition(aggregate, partition, refined);

    CheckEqual("Test A coarse vertex count", 4, num_vertices(aggregate.graph));
    CheckEqual("Test A coarse assignment 0", 0, coarse.community_of[0]);
    CheckEqual("Test A coarse assignment 1", 0, coarse.community_of[1]);
    CheckEqual("Test A coarse assignment 2", 1, coarse.community_of[2]);
    CheckEqual("Test A coarse assignment 3", 1, coarse.community_of[3]);
    CheckTrue("Test A grouping 0-1", coarse.community_of[0] == coarse.community_of[1]);
    CheckTrue("Test A grouping 2-3", coarse.community_of[2] == coarse.community_of[3]);
    CheckTrue("Test A groups differ", coarse.community_of[0] != coarse.community_of[2]);
    CheckStatsAgainstRecompute("Test A stats", aggregate, coarse);
}

void TestNonContiguousPartitionIds()
{
    const Graph G = MakeSimpleGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition partition =
        MakePartition(G, stats, {4, 4, 4, 9, 9, 9});
    const LeidenPartition refined =
        MakePartition(G, stats, {2, 2, 5, 7, 7, 9});
    const AggregateGraphResult aggregate = AggregateGraph(G, stats, refined);
    const LeidenPartition coarse =
        BuildCoarsePartition(aggregate, partition, refined);

    CheckEqual("Test B compact P id 4", 0, coarse.community_of[0]);
    CheckEqual("Test B compact P id 4 split", 0, coarse.community_of[1]);
    CheckEqual("Test B compact P id 9", 1, coarse.community_of[2]);
    CheckEqual("Test B compact P id 9 split", 1, coarse.community_of[3]);
}

void TestNonContiguousRefinedIds()
{
    const Graph G = MakeSimpleGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition partition =
        MakePartition(G, stats, {0, 0, 0, 1, 1, 1});
    const LeidenPartition refined =
        MakePartition(G, stats, {2, 2, 5, 7, 7, 9});
    const AggregateGraphResult aggregate = AggregateGraph(G, stats, refined);

    CheckEqual("Test C refined 2 coarse 0", 0, aggregate.coarse_of[0]);
    CheckEqual("Test C refined 5 coarse 1", 1, aggregate.coarse_of[2]);
    CheckEqual("Test C refined 7 coarse 2", 2, aggregate.coarse_of[3]);
    CheckEqual("Test C refined 9 coarse 3", 3, aggregate.coarse_of[5]);
    const LeidenPartition coarse =
        BuildCoarsePartition(aggregate, partition, refined);
    CheckEqual("Test C coarse assignment 0", 0, coarse.community_of[0]);
    CheckEqual("Test C coarse assignment 1", 0, coarse.community_of[1]);
    CheckEqual("Test C coarse assignment 2", 1, coarse.community_of[2]);
    CheckEqual("Test C coarse assignment 3", 1, coarse.community_of[3]);
}

void TestInvalidRefinement()
{
    Graph G = MakeGraph(4);
    add_undirected_edge(G, 0, 1, 1.0);
    add_undirected_edge(G, 1, 2, 1.0);
    add_undirected_edge(G, 2, 3, 1.0);
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition partition = MakePartition(G, stats, {0, 0, 1, 1});
    const LeidenPartition refined = MakePartition(G, stats, {5, 5, 5, 9});
    const AggregateGraphResult aggregate = AggregateGraph(G, stats, refined);

    try {
        (void)BuildCoarsePartition(aggregate, partition, refined);
        Fail("Test D invalid refinement",
             "std::invalid_argument",
             "no exception",
             "invalid refinement accepted");
    } catch (const std::invalid_argument&) {
    }
}

void TestSingletonRefinement()
{
    const Graph G = MakeSimpleGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition partition =
        MakePartition(G, stats, {0, 0, 0, 1, 1, 1});
    const LeidenPartition refined =
        MakePartition(G, stats, {0, 1, 2, 3, 4, 5});
    const AggregateGraphResult aggregate = AggregateGraph(G, stats, refined);
    const LeidenPartition coarse =
        BuildCoarsePartition(aggregate, partition, refined);

    CheckEqual("Test E aggregate vertex count", num_vertices(G),
               num_vertices(aggregate.graph));
    for (Vertex v = 0; v < num_vertices(G); ++v) {
        CheckEqual("Test E same grouping",
                   partition.community_of[v],
                   coarse.community_of[v]);
    }
}

void TestPrefinedEqualsPartition()
{
    const Graph G = MakeSimpleGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition partition =
        MakePartition(G, stats, {0, 0, 0, 1, 1, 1});
    const AggregateGraphResult aggregate = AggregateGraph(G, stats, partition);
    const LeidenPartition coarse =
        BuildCoarsePartition(aggregate, partition, partition);

    CheckEqual("Test F aggregate vertex count", 2, num_vertices(aggregate.graph));
    CheckEqual("Test F coarse singleton 0", 0, coarse.community_of[0]);
    CheckEqual("Test F coarse singleton 1", 1, coarse.community_of[1]);
}

void TestAggregateNodeSize()
{
    const Graph G = MakeSimpleGraph();
    const std::vector<double> node_size = {2.0, 3.0, 5.0, 7.0, 11.0, 13.0};
    const LeidenGraphStats stats = BuildLeidenGraphStats(G, node_size);
    const LeidenPartition partition =
        MakePartition(G, stats, {0, 0, 0, 1, 1, 1});
    const LeidenPartition refined =
        MakePartition(G, stats, {2, 2, 5, 7, 7, 9});
    const AggregateGraphResult aggregate = AggregateGraph(G, stats, refined);
    const LeidenPartition coarse =
        BuildCoarsePartition(aggregate, partition, refined);

    CheckNear("Test G community size 0", 10.0, coarse.community_size[0]);
    CheckNear("Test G community size 1", 31.0, coarse.community_size[1]);
    CheckNear("Test G total node size", Sum(stats.node_size),
              Sum(aggregate.stats.node_size));
}

void TestWeightedSelfLoop()
{
    const Graph G = MakeWeightedSelfLoopGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition partition =
        MakePartition(G, stats, {0, 0, 0, 1, 1, 1});
    const LeidenPartition refined =
        MakePartition(G, stats, {2, 2, 5, 7, 7, 9});
    const AggregateGraphResult aggregate = AggregateGraph(G, stats, refined);
    const LeidenPartition coarse =
        BuildCoarsePartition(aggregate, partition, refined);

    CheckStatsAgainstRecompute("Test H weighted self-loop stats", aggregate, coarse);
    CheckTrue("Test H nonzero strength", !coarse.community_strength.empty() &&
              coarse.community_strength[0] > 0.0 &&
              coarse.community_strength[1] > 0.0);
    CheckTrue("Test H nonzero internal edge", !coarse.internal_edge_weight.empty() &&
              coarse.internal_edge_weight[0] > 0.0 &&
              coarse.internal_edge_weight[1] > 0.0);
}

void TestQualityPreservation()
{
    const Graph G = MakeWeightedSelfLoopGraph();
    const std::vector<double> node_size = {2.0, 3.0, 5.0, 7.0, 11.0, 13.0};
    const LeidenGraphStats stats = BuildLeidenGraphStats(G, node_size);
    const LeidenPartition partition =
        MakePartition(G, stats, {0, 0, 0, 1, 1, 1});
    const LeidenPartition refined =
        MakePartition(G, stats, {2, 2, 5, 7, 7, 9});
    const CPMQualityFunction cpm(0.05);
    const ModularityQualityFunction modularity(1.1);

    CheckQualityPreservation("Test I CPM quality preservation",
                             G,
                             stats,
                             partition,
                             refined,
                             cpm);
    CheckQualityPreservation("Test I modularity quality preservation",
                             G,
                             stats,
                             partition,
                             refined,
                             modularity);
}

void TestDeterministicMapping()
{
    const Graph G = MakeSimpleGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition partition =
        MakePartition(G, stats, {4, 4, 4, 9, 9, 9});
    const LeidenPartition refined =
        MakePartition(G, stats, {2, 2, 5, 7, 7, 9});
    const AggregateGraphResult aggregate = AggregateGraph(G, stats, refined);
    const LeidenPartition a = BuildCoarsePartition(aggregate, partition, refined);
    const LeidenPartition b = BuildCoarsePartition(aggregate, partition, refined);

    CheckTrue("Test J deterministic assignment", a.community_of == b.community_of);
}

void TestInvalidInput()
{
    const Graph G = MakeSimpleGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition partition =
        MakePartition(G, stats, {0, 0, 0, 1, 1, 1});
    const LeidenPartition refined =
        MakePartition(G, stats, {2, 2, 5, 7, 7, 9});
    AggregateGraphResult aggregate = AggregateGraph(G, stats, refined);

    try {
        LeidenPartition bad_partition = partition;
        bad_partition.community_of.pop_back();
        (void)BuildCoarsePartition(aggregate, bad_partition, refined);
        Fail("Test K partition size", "std::invalid_argument", "no exception", "accepted");
    } catch (const std::invalid_argument&) {
    }

    try {
        LeidenPartition bad_refined = refined;
        bad_refined.community_of.pop_back();
        (void)BuildCoarsePartition(aggregate, partition, bad_refined);
        Fail("Test K refined size", "std::invalid_argument", "no exception", "accepted");
    } catch (const std::invalid_argument&) {
    }

    try {
        AggregateGraphResult bad_aggregate = aggregate;
        bad_aggregate.coarse_of[0] = -1;
        (void)BuildCoarsePartition(bad_aggregate, partition, refined);
        Fail("Test K coarse negative", "std::invalid_argument", "no exception", "accepted");
    } catch (const std::invalid_argument&) {
    }

    try {
        AggregateGraphResult bad_aggregate = aggregate;
        bad_aggregate.coarse_of[0] = num_vertices(aggregate.graph);
        (void)BuildCoarsePartition(bad_aggregate, partition, refined);
        Fail("Test K coarse out of range", "std::invalid_argument", "no exception", "accepted");
    } catch (const std::invalid_argument&) {
    }

    try {
        LeidenPartition bad_partition = partition;
        bad_partition.community_of[0] = -1;
        (void)BuildCoarsePartition(aggregate, bad_partition, refined);
        Fail("Test K negative P community", "std::invalid_argument", "no exception", "accepted");
    } catch (const std::invalid_argument&) {
    }
}

} // namespace

int main()
{
    TestBasicCoarsePartition();
    TestNonContiguousPartitionIds();
    TestNonContiguousRefinedIds();
    TestInvalidRefinement();
    TestSingletonRefinement();
    TestPrefinedEqualsPartition();
    TestAggregateNodeSize();
    TestWeightedSelfLoop();
    TestQualityPreservation();
    TestDeterministicMapping();
    TestInvalidInput();

    std::cout << "All coarse-partition tests passed.\n";
    return EXIT_SUCCESS;
}
