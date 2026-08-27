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

void CheckEqual(const std::string& test_name, int expected, int actual)
{
    if (expected != actual) {
        Fail(test_name,
             std::to_string(expected),
             std::to_string(actual),
             std::to_string(actual - expected));
    }
}

double EdgeWeightInAggregate(const Graph& G, int u, int v)
{
    double weight = 0.0;
    for (const Edge& e : G.adj[u]) {
        if (e.to == v) {
            weight += e.weight;
        }
    }
    return weight;
}

double Sum(const std::vector<double>& values)
{
    return std::accumulate(values.begin(), values.end(), 0.0);
}

void CheckStatsPreservation(const std::string& test_name,
                            const Graph& G,
                            const LeidenGraphStats& stats,
                            const AggregateGraphResult& aggregate)
{
    CheckNear(test_name + " total edge weight",
              stats.total_edge_weight,
              aggregate.stats.total_edge_weight);
    CheckNear(test_name + " node size sum",
              Sum(stats.node_size),
              Sum(aggregate.stats.node_size));
    CheckNear(test_name + " strength sum",
              2.0 * aggregate.stats.total_edge_weight,
              Sum(aggregate.stats.node_strength));

    std::vector<double> expected_strength(num_vertices(aggregate.graph), 0.0);
    for (Vertex v = 0; v < num_vertices(G); ++v) {
        expected_strength[aggregate.coarse_of[v]] += stats.node_strength[v];
    }
    for (Vertex c = 0; c < num_vertices(aggregate.graph); ++c) {
        CheckNear(test_name + " coarse strength",
                  expected_strength[c],
                  aggregate.stats.node_strength[c]);
    }
}

void CheckQualityPreservation(const std::string& test_name,
                              const Graph& G,
                              const LeidenGraphStats& stats,
                              const LeidenPartition& refined,
                              const QualityFunction& qf)
{
    const AggregateGraphResult aggregate = AggregateGraph(G, stats, refined);
    const LeidenPartition singleton =
        MakeSingletonPartition(aggregate.graph, aggregate.stats);

    // Each aggregate singleton represents one refined community, so this
    // coarse singleton partition is the same partition as Prefined on G.
    CheckNear(test_name,
              qf.quality(G, stats, refined),
              qf.quality(aggregate.graph, aggregate.stats, singleton));
}

Graph MakeSimpleGraph()
{
    Graph G = MakeGraph(4);
    add_undirected_edge(G, 0, 1);
    add_undirected_edge(G, 1, 2);
    add_undirected_edge(G, 2, 3);
    add_undirected_edge(G, 0, 3);
    return G;
}

Graph MakeWeightedGraph()
{
    Graph G = MakeGraph(5);
    add_undirected_edge(G, 0, 1, 2.0);
    add_undirected_edge(G, 0, 2, 1.0);
    add_undirected_edge(G, 1, 3, 2.0);
    add_undirected_edge(G, 2, 4, 3.0);
    add_undirected_edge(G, 3, 4, 4.0);
    return G;
}

Graph MakeSelfLoopGraph()
{
    Graph G = MakeGraph(4);
    add_undirected_edge(G, 0, 0, 5.0);
    add_undirected_edge(G, 0, 1, 2.0);
    add_undirected_edge(G, 1, 1, 7.0);
    add_undirected_edge(G, 1, 2, 3.0);
    add_undirected_edge(G, 2, 3, 4.0);
    add_undirected_edge(G, 3, 3, 11.0);
    return G;
}

void TestSimpleUnweighted()
{
    const Graph G = MakeSimpleGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition refined = MakePartition(G, stats, {0, 0, 1, 1});
    const AggregateGraphResult aggregate = AggregateGraph(G, stats, refined);

    CheckEqual("Test A aggregate vertex count", 2, num_vertices(aggregate.graph));
    CheckTrue("Test A coarse_of 0-1",
              aggregate.coarse_of[0] == aggregate.coarse_of[1]);
    CheckTrue("Test A coarse_of 2-3",
              aggregate.coarse_of[2] == aggregate.coarse_of[3]);
    CheckTrue("Test A coarse groups differ",
              aggregate.coarse_of[0] != aggregate.coarse_of[2]);
    CheckNear("Test A self-loop 0", 1.0,
              EdgeWeightInAggregate(aggregate.graph, 0, 0));
    CheckNear("Test A self-loop 1", 1.0,
              EdgeWeightInAggregate(aggregate.graph, 1, 1));
    CheckNear("Test A inter edge", 2.0,
              EdgeWeightInAggregate(aggregate.graph, 0, 1));
    CheckNear("Test A node size 0", 2.0, aggregate.stats.node_size[0]);
    CheckNear("Test A node size 1", 2.0, aggregate.stats.node_size[1]);
    CheckStatsPreservation("Test A preservation", G, stats, aggregate);
}

void TestWeightedParallelAggregation()
{
    const Graph G = MakeWeightedGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition refined = MakePartition(G, stats, {0, 0, 1, 1, 1});
    const AggregateGraphResult aggregate = AggregateGraph(G, stats, refined);

    CheckNear("Test B internal self-loop coarse 0", 2.0,
              EdgeWeightInAggregate(aggregate.graph, 0, 0));
    CheckNear("Test B internal self-loop coarse 1", 7.0,
              EdgeWeightInAggregate(aggregate.graph, 1, 1));
    CheckNear("Test B parallel edges summed", 3.0,
              EdgeWeightInAggregate(aggregate.graph, 0, 1));
    CheckStatsPreservation("Test B preservation", G, stats, aggregate);
}

void TestInternalEdgesAndOriginalSelfLoops()
{
    const Graph G = MakeSelfLoopGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition refined = MakePartition(G, stats, {0, 0, 1, 1});
    const AggregateGraphResult aggregate = AggregateGraph(G, stats, refined);

    CheckNear("Test C coarse 0 self-loop", 14.0,
              EdgeWeightInAggregate(aggregate.graph, 0, 0));
    CheckNear("Test C coarse 1 self-loop", 15.0,
              EdgeWeightInAggregate(aggregate.graph, 1, 1));
    CheckNear("Test C inter edge", 3.0,
              EdgeWeightInAggregate(aggregate.graph, 0, 1));
    CheckStatsPreservation("Test C preservation", G, stats, aggregate);
}

void TestNonContiguousCommunityIds()
{
    Graph G = MakeGraph(5);
    add_undirected_edge(G, 0, 1, 1.0);
    add_undirected_edge(G, 1, 2, 2.0);
    add_undirected_edge(G, 3, 4, 3.0);
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition refined = MakePartition(G, stats, {2, 2, 5, 5, 9});
    const AggregateGraphResult aggregate = AggregateGraph(G, stats, refined);

    CheckEqual("Test D compact vertex count", 3, num_vertices(aggregate.graph));
    CheckEqual("Test D community 2 compacted", 0, aggregate.coarse_of[0]);
    CheckEqual("Test D community 2 compacted pair", 0, aggregate.coarse_of[1]);
    CheckEqual("Test D community 5 compacted", 1, aggregate.coarse_of[2]);
    CheckEqual("Test D community 5 compacted pair", 1, aggregate.coarse_of[3]);
    CheckEqual("Test D community 9 compacted", 2, aggregate.coarse_of[4]);
}

void TestAggregateNodeSize()
{
    const Graph G = MakeWeightedGraph();
    const std::vector<double> node_size = {2.0, 3.0, 1.0, 4.0, 5.0};
    const LeidenGraphStats stats = BuildLeidenGraphStats(G, node_size);
    const LeidenPartition refined = MakePartition(G, stats, {0, 0, 1, 1, 1});
    const AggregateGraphResult aggregate = AggregateGraph(G, stats, refined);

    CheckNear("Test E aggregate node size 0", 5.0,
              aggregate.stats.node_size[0]);
    CheckNear("Test E aggregate node size 1", 10.0,
              aggregate.stats.node_size[1]);
    CheckStatsPreservation("Test E preservation", G, stats, aggregate);
}

void TestTotalEdgeWeightPreservation()
{
    const Graph unweighted = MakeSimpleGraph();
    const LeidenGraphStats unweighted_stats = BuildLeidenGraphStats(unweighted);
    CheckStatsPreservation("Test F unweighted",
                           unweighted,
                           unweighted_stats,
                           AggregateGraph(unweighted,
                                          unweighted_stats,
                                          MakePartition(unweighted,
                                                        unweighted_stats,
                                                        {0, 0, 1, 1})));

    const Graph weighted = MakeWeightedGraph();
    const LeidenGraphStats weighted_stats = BuildLeidenGraphStats(weighted);
    CheckStatsPreservation("Test F weighted",
                           weighted,
                           weighted_stats,
                           AggregateGraph(weighted,
                                          weighted_stats,
                                          MakePartition(weighted,
                                                        weighted_stats,
                                                        {0, 0, 1, 1, 1})));

    const Graph self_loop = MakeSelfLoopGraph();
    const LeidenGraphStats self_loop_stats = BuildLeidenGraphStats(self_loop);
    CheckStatsPreservation("Test F self-loop",
                           self_loop,
                           self_loop_stats,
                           AggregateGraph(self_loop,
                                          self_loop_stats,
                                          MakePartition(self_loop,
                                                        self_loop_stats,
                                                        {0, 0, 1, 1})));
}

void TestStrengthPreservation()
{
    const Graph G = MakeSelfLoopGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition refined = MakePartition(G, stats, {0, 0, 1, 1});
    CheckStatsPreservation("Test G strength", G, stats,
                           AggregateGraph(G, stats, refined));
}

void TestCPMQualityPreservation()
{
    const Graph G = MakeWeightedGraph();
    const std::vector<double> node_size = {2.0, 3.0, 1.0, 4.0, 5.0};
    const LeidenGraphStats stats = BuildLeidenGraphStats(G, node_size);
    const LeidenPartition refined = MakePartition(G, stats, {0, 0, 1, 1, 1});
    const CPMQualityFunction cpm(0.35);
    CheckQualityPreservation("Test H CPM quality preservation",
                             G,
                             stats,
                             refined,
                             cpm);
}

void TestModularityQualityPreservation()
{
    const Graph G = MakeSelfLoopGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition refined = MakePartition(G, stats, {0, 0, 1, 1});
    const ModularityQualityFunction modularity(1.2);
    CheckQualityPreservation("Test I modularity quality preservation",
                             G,
                             stats,
                             refined,
                             modularity);
}

void TestRepeatedAggregationReadiness()
{
    const Graph G = MakeWeightedGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition refined = MakePartition(G, stats, {0, 0, 1, 1, 1});
    const AggregateGraphResult first = AggregateGraph(G, stats, refined);
    const LeidenPartition second_partition =
        MakePartition(first.graph, first.stats, {2, 2});
    const AggregateGraphResult second =
        AggregateGraph(first.graph, first.stats, second_partition);

    CheckEqual("Test J second aggregate vertex count", 1,
               num_vertices(second.graph));
    CheckNear("Test J total edge preserved",
              stats.total_edge_weight,
              second.stats.total_edge_weight);
    CheckNear("Test J node size preserved",
              Sum(stats.node_size),
              Sum(second.stats.node_size));
}

void TestInvalidInput()
{
    const Graph G = MakeSimpleGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);

    try {
        const LeidenPartition bad_size =
            MakePartition(G, stats, {0, 0, 1, 1});
        LeidenPartition copied = bad_size;
        copied.community_of.pop_back();
        (void)AggregateGraph(G, stats, copied);
        Fail("Test K invalid size", "std::invalid_argument", "no exception", "accepted");
    } catch (const std::invalid_argument&) {
    }

    try {
        LeidenPartition bad_negative =
            MakePartition(G, stats, {0, 0, 1, 1});
        bad_negative.community_of[2] = -1;
        (void)AggregateGraph(G, stats, bad_negative);
        Fail("Test K negative community", "std::invalid_argument", "no exception", "accepted");
    } catch (const std::invalid_argument&) {
    }
}

void TestDeterministicAggregation()
{
    const Graph G = MakeSelfLoopGraph();
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition refined = MakePartition(G, stats, {4, 4, 7, 7});
    const AggregateGraphResult a = AggregateGraph(G, stats, refined);
    const AggregateGraphResult b = AggregateGraph(G, stats, refined);

    CheckTrue("Test L deterministic coarse_of", a.coarse_of == b.coarse_of);
    CheckEqual("Test L deterministic vertex count",
               num_vertices(a.graph),
               num_vertices(b.graph));
    for (Vertex u = 0; u < num_vertices(a.graph); ++u) {
        CheckEqual("Test L deterministic degree",
                   static_cast<int>(a.graph.adj[u].size()),
                   static_cast<int>(b.graph.adj[u].size()));
        for (std::size_t i = 0; i < a.graph.adj[u].size(); ++i) {
            CheckEqual("Test L deterministic edge target",
                       a.graph.adj[u][i].to,
                       b.graph.adj[u][i].to);
            CheckNear("Test L deterministic edge weight",
                      a.graph.adj[u][i].weight,
                      b.graph.adj[u][i].weight);
        }
    }
}

} // namespace

int main()
{
    TestSimpleUnweighted();
    TestWeightedParallelAggregation();
    TestInternalEdgesAndOriginalSelfLoops();
    TestNonContiguousCommunityIds();
    TestAggregateNodeSize();
    TestTotalEdgeWeightPreservation();
    TestStrengthPreservation();
    TestCPMQualityPreservation();
    TestModularityQualityPreservation();
    TestRepeatedAggregationReadiness();
    TestInvalidInput();
    TestDeterministicAggregation();

    std::cout << "All AggregateGraph tests passed.\n";
    return EXIT_SUCCESS;
}
