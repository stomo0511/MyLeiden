#include "Leiden.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

[[noreturn]] void Fail(const std::string& message)
{
    std::cerr << "Stage4B isolated test failed: " << message << "\n";
    std::exit(EXIT_FAILURE);
}

void ValidateResult(const Graph& graph,
                    const LeidenGraphStats& stats,
                    const QualityFunction& quality,
                    const MoveNodesFastResult& result)
{
    if (result.partition.community_of.size() !=
        static_cast<std::size_t>(num_vertices(graph))) {
        Fail("lost vertex");
    }
    for (Community community : result.partition.community_of) {
        if (community < 0 ||
            static_cast<std::size_t>(community) >=
                result.partition.community_size.size()) {
            Fail("invalid community id");
        }
    }
    const LeidenPartition rebuilt =
        MakePartition(graph, stats, result.partition.community_of);
    if (rebuilt.community_of != result.partition.community_of ||
        rebuilt.community_is_empty != result.partition.community_is_empty) {
        Fail("partition rebuild mismatch");
    }
    const auto equal_values = [](const std::vector<double>& lhs,
                                 const std::vector<double>& rhs) {
        if (lhs.size() != rhs.size()) return false;
        for (std::size_t i = 0; i < lhs.size(); ++i) {
            if (std::abs(lhs[i] - rhs[i]) > 1.0e-9) return false;
        }
        return true;
    };
    if (!equal_values(rebuilt.community_size,
                      result.partition.community_size) ||
        !equal_values(rebuilt.community_strength,
                      result.partition.community_strength) ||
        !equal_values(rebuilt.internal_edge_weight,
                      result.partition.internal_edge_weight)) {
        Fail("community statistics mismatch");
    }
    if (!std::isfinite(quality.quality(graph, stats, result.partition))) {
        Fail("non-finite quality");
    }
    for (Vertex v = 0; v < num_vertices(graph); ++v) {
        const Community source = result.partition.community_of[v];
        const auto weights =
            BuildNeighborCommunityWeights(graph, result.partition, v);
        const auto lookup = [&weights](Community community) {
            const auto it = weights.find(community);
            return it == weights.end() ? 0.0 : it->second;
        };
        for (Community target = 0;
             target <= static_cast<Community>(
                 result.partition.community_size.size()); ++target) {
            if (target == source) continue;
            const double delta = quality.deltaMoveFromWeights(
                stats, result.partition, v, target,
                lookup(source), lookup(target));
            if (delta > 1.0e-10) Fail("positive move after verification sweep");
        }
    }
}

void Run(const Graph& graph,
         const LeidenPartition& initial,
         const QualityFunction& quality,
         std::uint64_t seed)
{
    const LeidenGraphStats stats = BuildLeidenGraphStats(graph);
    const MoveNodesFastResult result = MoveNodesFastParallelStage4B(
        graph, stats, initial, quality, seed, 1);
    ValidateResult(graph, stats, quality, result);
}

Graph Chain(int n)
{
    Graph graph = MakeGraph(n);
    for (Vertex v = 0; v + 1 < n; ++v) {
        add_undirected_edge(graph, v, v + 1, 1.0 + 0.1 * v);
    }
    return graph;
}

void RunSuite(std::uint64_t seed)
{
    const CPMQualityFunction cpm(0.1);
    const ModularityQualityFunction modularity(1.0);

    Graph single = MakeGraph(1);
    LeidenGraphStats single_stats = BuildLeidenGraphStats(single);
    Run(single, MakeSingletonPartition(single, single_stats), cpm, seed);

    Graph no_edges = MakeGraph(8);
    LeidenGraphStats no_edge_stats = BuildLeidenGraphStats(no_edges);
    Run(no_edges,
        MakePartition(no_edges, no_edge_stats, {0, 0, 0, 0, 0, 0, 0, 0}),
        cpm, seed + 1);

    Graph chain = Chain(12);
    LeidenGraphStats chain_stats = BuildLeidenGraphStats(chain);
    Run(chain, MakeSingletonPartition(chain, chain_stats), modularity, seed + 2);

    Graph disconnected = MakeGraph(10);
    add_undirected_edge(disconnected, 0, 1, 2.0);
    add_undirected_edge(disconnected, 1, 2, 1.0);
    add_undirected_edge(disconnected, 3, 4, 3.0);
    add_undirected_edge(disconnected, 5, 6, 1.5);
    add_undirected_edge(disconnected, 7, 8, 0.75);
    LeidenGraphStats disconnected_stats =
        BuildLeidenGraphStats(disconnected);
    Run(disconnected,
        MakeSingletonPartition(disconnected, disconnected_stats),
        modularity, seed + 3);

    // Shared source, shared target, and opposite A<->B proposals are all
    // possible in this dense two-community contention graph.
    Graph contention = MakeGraph(8);
    for (Vertex u = 0; u < 8; ++u) {
        for (Vertex v = u + 1; v < 8; ++v) {
            add_undirected_edge(contention, u, v,
                                ((u + v) % 3 == 0) ? 2.0 : 0.5);
        }
    }
    add_undirected_edge(contention, 0, 0, 0.25);
    add_undirected_edge(contention, 0, 1, 0.75); // parallel edge
    LeidenGraphStats contention_stats = BuildLeidenGraphStats(contention);
    Run(contention,
        MakePartition(contention, contention_stats,
                      {0, 0, 0, 0, 1, 1, 1, 1}),
        cpm, seed + 4);
}

} // namespace

int main()
{
    int iterations = 1;
    if (const char* text = std::getenv("STAGE4B_STRESS_ITERATIONS")) {
        iterations = std::max(1, std::atoi(text));
    }
    for (int i = 0; i < iterations; ++i) {
        RunSuite(static_cast<std::uint64_t>(1000 + i));
    }
    std::cout << "Stage4B isolated tests passed (iterations="
              << iterations << ").\n";
    return EXIT_SUCCESS;
}
