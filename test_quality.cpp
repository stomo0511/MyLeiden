#include "QualityFunction.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr double kTolerance = 1.0e-10;

void CheckNear(const std::string& test_name, double expected, double actual)
{
    const double diff = std::abs(expected - actual);
    if (diff > kTolerance) {
        std::cerr << test_name << " failed\n"
                  << "expected value: " << expected << "\n"
                  << "actual value:   " << actual << "\n"
                  << "difference:     " << diff << "\n";
        std::exit(EXIT_FAILURE);
    }
}

void CheckEqual(const std::string& test_name, double expected, double actual)
{
    if (expected != actual) {
        std::cerr << test_name << " failed\n"
                  << "expected value: " << expected << "\n"
                  << "actual value:   " << actual << "\n"
                  << "difference:     " << (actual - expected) << "\n";
        std::exit(EXIT_FAILURE);
    }
}

void CheckDeltaByRecompute(const std::string& test_name,
                           const Graph& G,
                           const QualityFunction& qf,
                           const std::vector<Community>& assignment,
                           Vertex v,
                           Community target,
                           const std::vector<double>& node_size = {})
{
    const LeidenGraphStats stats =
        node_size.empty() ? BuildLeidenGraphStats(G)
                          : BuildLeidenGraphStats(G, node_size);
    const LeidenPartition before = MakePartition(G, stats, assignment);

    std::vector<Community> moved_assignment = assignment;
    moved_assignment[v] = target;
    const LeidenPartition after = MakePartition(G, stats, moved_assignment);

    const double expected =
        qf.quality(G, stats, after) - qf.quality(G, stats, before);
    const double actual = qf.deltaMove(G, stats, before, v, target);
    CheckNear(test_name, expected, actual);
}

void CheckActualMove(const std::string& test_name,
                     const Graph& G,
                     const QualityFunction& qf,
                     const std::vector<Community>& assignment,
                     Vertex v,
                     Community target)
{
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    LeidenPartition partition = MakePartition(G, stats, assignment);
    const double before_quality = qf.quality(G, stats, partition);
    const double expected_delta = qf.deltaMove(G, stats, partition, v, target);

    MoveNodeToCommunity(G, stats, partition, v, target);
    const double actual_delta = qf.quality(G, stats, partition) - before_quality;
    CheckNear(test_name, expected_delta, actual_delta);
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

void TestUnweighted()
{
    const Graph G = MakeUnweightedGraph();
    const std::vector<Community> p = {0, 0, 0, 1, 1, 1};
    const CPMQualityFunction cpm(0.4);
    const ModularityQualityFunction modularity(1.0);

    CheckDeltaByRecompute("Test 1 CPM unweighted delta", G, cpm, p, 2, 1);
    CheckDeltaByRecompute("Test 1 modularity unweighted delta", G, modularity, p, 2, 1);
}

void TestWeighted()
{
    const Graph G = MakeWeightedGraph();
    const std::vector<Community> p = {0, 0, 1, 1, 1};
    const CPMQualityFunction cpm(0.35);
    const ModularityQualityFunction modularity(0.8);

    CheckDeltaByRecompute("Test 2 CPM weighted delta", G, cpm, p, 1, 1);
    CheckDeltaByRecompute("Test 2 modularity weighted delta", G, modularity, p, 1, 1);
}

void TestSelfLoop()
{
    const Graph G = MakeSelfLoopGraph();
    const std::vector<Community> p = {0, 0, 1, 1};
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition partition = MakePartition(G, stats, p);
    const CPMQualityFunction cpm(0.2);
    const ModularityQualityFunction modularity(1.2);

    CheckEqual("Test 3 total edge weight with self-loop", 9.5, stats.total_edge_weight);
    CheckEqual("Test 3 vertex strength with self-loop v0", 5.5, stats.node_strength[0]);
    CheckEqual("Test 3 vertex strength with self-loop v2", 9.0, stats.node_strength[2]);
    CheckEqual("Test 3 internal edge weight community 0", 3.5,
               partition.internal_edge_weight[0]);
    CheckEqual("Test 3 internal edge weight community 1", 5.5,
               partition.internal_edge_weight[1]);

    CheckDeltaByRecompute("Test 3 CPM self-loop delta", G, cpm, p, 1, 1);
    CheckDeltaByRecompute("Test 3 modularity self-loop delta", G, modularity, p, 1, 1);
}

void TestActualMove()
{
    const Graph G = MakeWeightedGraph();
    const std::vector<Community> p = {0, 0, 1, 1, 1};
    const CPMQualityFunction cpm(0.3);
    const ModularityQualityFunction modularity(1.0);

    CheckActualMove("Test 4 CPM actual node move", G, cpm, p, 0, 1);
    CheckActualMove("Test 4 modularity actual node move", G, modularity, p, 0, 1);
}

void TestAggregateNodeSize()
{
    const Graph G = MakeWeightedGraph();
    const std::vector<double> node_size = {2.0, 1.0, 3.0, 1.0, 2.0};
    const std::vector<Community> p = {0, 0, 1, 1, 1};
    const CPMQualityFunction cpm(0.15);

    CheckDeltaByRecompute("Aggregate node-size CPM delta", G, cpm, p, 0, 1, node_size);
}

} // namespace

int main()
{
    TestUnweighted();
    TestWeighted();
    TestSelfLoop();
    TestActualMove();
    TestAggregateNodeSize();

    std::cout << "All quality-function tests passed.\n";
    return EXIT_SUCCESS;
}
