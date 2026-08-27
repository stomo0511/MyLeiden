#include "common/Block_Eval.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr double kTolerance = 1.0e-10;

void Fail(const std::string& test_name,
          const std::string& expected,
          const std::string& actual)
{
    std::cerr << test_name << " failed\n"
              << "expected value: " << expected << "\n"
              << "actual value:   " << actual << "\n";
    std::exit(EXIT_FAILURE);
}

void CheckTrue(const std::string& test_name, bool condition)
{
    if (!condition) {
        Fail(test_name, "true", "false");
    }
}

void CheckEqual(const std::string& test_name, int expected, int actual)
{
    if (expected != actual) {
        Fail(test_name, std::to_string(expected), std::to_string(actual));
    }
}

void CheckNear(const std::string& test_name, double expected, double actual)
{
    const double diff = std::abs(expected - actual);
    if (diff > kTolerance) {
        Fail(test_name,
             std::to_string(expected),
             std::to_string(actual));
    }
}

Graph MakeSquareGraph()
{
    Graph G = MakeGraph(4);
    add_undirected_edge(G, 0, 1);
    add_undirected_edge(G, 0, 2);
    add_undirected_edge(G, 1, 3);
    add_undirected_edge(G, 2, 3);
    return G;
}

void TestSimpleGraphMetrics()
{
    const Graph G = MakeSquareGraph();
    const std::vector<int> block_of = {0, 0, 1, 1};

    const std::vector<int> nodes = CountNodesPerBlock(block_of, 2);
    CheckEqual("nodes block 0", 2, nodes[0]);
    CheckEqual("nodes block 1", 2, nodes[1]);

    const std::vector<int> internal = CountInternalEdges(G, block_of, 2);
    CheckEqual("internal block 0", 1, internal[0]);
    CheckEqual("internal block 1", 1, internal[1]);

    int internal_edges = 0;
    int external_edges = 0;
    CountEdgeLocality(G, block_of, 2, internal_edges, external_edges);
    CheckEqual("locality internal", 2, internal_edges);
    CheckEqual("locality external", 2, external_edges);
    CheckNear("internal ratio", 0.5,
              InternalEdgeRatio(internal_edges, external_edges));

    const double modularity = Modularity_Unweighted(G, block_of);
    CheckTrue("unweighted modularity finite", std::isfinite(modularity));

    const BlockEvaluationResult metrics =
        EvaluatePartitioningMetrics(G, block_of, 2);
    CheckEqual("metrics vertices", 4, metrics.num_vertices);
    CheckEqual("metrics edges", 4, metrics.num_edges);
    CheckEqual("metrics block graph edges", 1, metrics.block_graph_edges);
}

void TestWeightedModularity()
{
    Graph G = MakeGraph(4);
    add_undirected_edge(G, 0, 1, -2.0);
    add_undirected_edge(G, 1, 2, 3.0);
    add_undirected_edge(G, 2, 3, -4.0);

    const std::vector<int> block_of = {0, 0, 1, 1};
    const double modularity = Modularity_Weighted(G, block_of);
    CheckTrue("weighted modularity finite", std::isfinite(modularity));
}

void TestBinaryBlockGraph()
{
    const Graph G = MakeSquareGraph();
    const std::vector<int> block_of = {0, 0, 1, 1};
    const Graph T = BuildBlockGraph(G, block_of, BlockEdgeWeight::Binary);
    const std::vector<int> block_degree = BlockDegreesBinary(T);

    CheckEqual("block graph vertices", 2, num_vertices(T));
    CheckEqual("block graph edges", 1, static_cast<int>(num_edges(T)));
    CheckEqual("block degree 0", 1, block_degree[0]);
    CheckEqual("block degree 1", 1, block_degree[1]);
    CheckNear("block graph density", 1.0, BlockGraphDensity(T));
}

} // namespace

int main()
{
    TestSimpleGraphMetrics();
    TestWeightedModularity();
    TestBinaryBlockGraph();

    std::cout << "All Block_Eval tests passed.\n";
    return EXIT_SUCCESS;
}
