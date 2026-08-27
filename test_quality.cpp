#include "QualityFunction.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
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

double LookupWeight(const std::unordered_map<Community, double>& weights,
                    Community community)
{
    const auto it = weights.find(community);
    return (it == weights.end()) ? 0.0 : it->second;
}

int MaxCommunityId(const std::vector<Community>& assignment)
{
    return *std::max_element(assignment.begin(), assignment.end());
}

void CheckPartitionStatsAgainstRecompute(const std::string& test_name,
                                         const Graph& G,
                                         const LeidenGraphStats& stats,
                                         const LeidenPartition& partition)
{
    const LeidenPartition recomputed =
        MakePartition(G, stats, partition.community_of);

    if (partition.community_of != recomputed.community_of) {
        std::cerr << test_name << " failed\n"
                  << "expected value: recomputed community_of\n"
                  << "actual value:   updated community_of\n"
                  << "difference:     assignment differs\n";
        std::exit(EXIT_FAILURE);
    }

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

void CheckExhaustiveMoves(const std::string& test_name,
                          const Graph& G,
                          const QualityFunction& qf,
                          const std::vector<Community>& assignment,
                          const std::vector<double>& node_size = {})
{
    const LeidenGraphStats stats =
        node_size.empty() ? BuildLeidenGraphStats(G)
                          : BuildLeidenGraphStats(G, node_size);
    const LeidenPartition before = MakePartition(G, stats, assignment);
    const double before_quality = qf.quality(G, stats, before);
    const Community new_community = MaxCommunityId(assignment) + 1;

    for (Vertex v = 0; v < num_vertices(G); ++v) {
        const Community source = before.community_of[v];
        const auto weights = BuildNeighborCommunityWeights(G, before, v);
        const double weight_to_source = LookupWeight(weights, source);

        for (Community target = 0; target <= new_community; ++target) {
            std::vector<Community> moved_assignment = assignment;
            moved_assignment[v] = target;
            const LeidenPartition after =
                MakePartition(G, stats, moved_assignment);
            const double expected =
                qf.quality(G, stats, after) - before_quality;
            const double from_delta =
                qf.deltaMove(G, stats, before, v, target);
            const double from_weights =
                qf.deltaMoveFromWeights(stats,
                                        before,
                                        v,
                                        target,
                                        weight_to_source,
                                        LookupWeight(weights, target));

            const std::string label =
                test_name + " v=" + std::to_string(v)
                + " target=" + std::to_string(target);
            CheckNear(label + " deltaMove", expected, from_delta);
            CheckNear(label + " deltaMoveFromWeights", expected, from_weights);

            if (target == source) {
                CheckNear(label + " target==source", 0.0, from_delta);
                CheckNear(label + " target==source fast", 0.0, from_weights);
            }
        }
    }
}

void CheckActualMoveSequence(const std::string& test_name,
                             const Graph& G,
                             const QualityFunction& qf,
                             const std::vector<Community>& assignment,
                             const std::vector<std::pair<Vertex, Community>>& moves)
{
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    LeidenPartition partition = MakePartition(G, stats, assignment);

    for (std::size_t i = 0; i < moves.size(); ++i) {
        const Vertex v = moves[i].first;
        const Community target = moves[i].second;
        const double before_quality = qf.quality(G, stats, partition);
        const double expected_delta = qf.deltaMove(G, stats, partition, v, target);

        MoveNodeToCommunity(G, stats, partition, v, target);

        const double actual_delta =
            qf.quality(G, stats, partition) - before_quality;
        const std::string label =
            test_name + " move=" + std::to_string(i);
        CheckNear(label, expected_delta, actual_delta);
        CheckPartitionStatsAgainstRecompute(label, G, stats, partition);
    }
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

    CheckExhaustiveMoves("Test A CPM unweighted exhaustive", G, cpm, p);
    CheckExhaustiveMoves("Test A modularity unweighted exhaustive", G, modularity, p);
}

void TestWeighted()
{
    const Graph G = MakeWeightedGraph();
    const std::vector<Community> p = {0, 0, 1, 1, 1};
    const CPMQualityFunction cpm(0.35);
    const ModularityQualityFunction modularity(0.8);

    CheckExhaustiveMoves("Test B CPM weighted exhaustive", G, cpm, p);
    CheckExhaustiveMoves("Test B modularity weighted exhaustive", G, modularity, p);
}

void TestSelfLoop()
{
    const Graph G = MakeSelfLoopGraph();
    const std::vector<Community> p = {0, 0, 1, 1};
    const LeidenGraphStats stats = BuildLeidenGraphStats(G);
    const LeidenPartition partition = MakePartition(G, stats, p);
    const CPMQualityFunction cpm(0.2);
    const ModularityQualityFunction modularity(1.2);

    CheckEqual("Test C total edge weight with self-loop", 9.5,
               stats.total_edge_weight);
    CheckEqual("Test C vertex strength with self-loop v0", 5.5,
               stats.node_strength[0]);
    CheckEqual("Test C vertex strength with self-loop v2", 9.0,
               stats.node_strength[2]);
    CheckEqual("Test C internal edge weight community 0", 3.5,
               partition.internal_edge_weight[0]);
    CheckEqual("Test C internal edge weight community 1", 5.5,
               partition.internal_edge_weight[1]);

    CheckExhaustiveMoves("Test C CPM self-loop exhaustive", G, cpm, p);
    CheckExhaustiveMoves("Test C modularity self-loop exhaustive", G, modularity, p);
}

void TestSingletonSource()
{
    const Graph G = MakeWeightedGraph();
    const std::vector<Community> p = {0, 1, 1, 2, 2};
    const CPMQualityFunction cpm(0.25);
    const ModularityQualityFunction modularity(1.0);

    CheckExhaustiveMoves("Test D CPM singleton-source exhaustive", G, cpm, p);
    CheckExhaustiveMoves("Test D modularity singleton-source exhaustive", G, modularity, p);
    CheckActualMoveSequence("Test D CPM singleton-source move",
                            G,
                            cpm,
                            p,
                            {{0, 1}});
    CheckActualMoveSequence("Test D modularity singleton-source move",
                            G,
                            modularity,
                            p,
                            {{0, 1}});
}

void TestAggregateNodeSize()
{
    const Graph G = MakeWeightedGraph();
    const std::vector<double> node_size = {2.0, 1.0, 3.0, 1.0, 2.0};
    const std::vector<Community> p = {0, 0, 1, 1, 1};
    const CPMQualityFunction cpm(0.15);
    const ModularityQualityFunction modularity(1.1);

    CheckExhaustiveMoves("Test E CPM aggregate-size exhaustive",
                         G,
                         cpm,
                         p,
                         node_size);
    CheckExhaustiveMoves("Test E modularity aggregate-size exhaustive",
                         G,
                         modularity,
                         p,
                         node_size);
}

void TestNewTargetCommunity()
{
    const Graph G = MakeUnweightedGraph();
    const std::vector<Community> p = {0, 0, 1, 1, 2, 2};
    const CPMQualityFunction cpm(0.3);
    const ModularityQualityFunction modularity(0.9);

    CheckExhaustiveMoves("Test F CPM new-target exhaustive", G, cpm, p);
    CheckExhaustiveMoves("Test F modularity new-target exhaustive", G, modularity, p);
}

void TestActualMoveSequence()
{
    const Graph G = MakeWeightedGraph();
    const std::vector<Community> p = {0, 0, 1, 1, 1};
    const std::vector<std::pair<Vertex, Community>> moves = {
        {0, 1},
        {2, 0},
        {4, 2},
        {1, 2}
    };
    const CPMQualityFunction cpm(0.3);
    const ModularityQualityFunction modularity(1.0);

    CheckActualMoveSequence("Test G CPM sequential actual moves", G, cpm, p, moves);
    CheckActualMoveSequence("Test G modularity sequential actual moves",
                            G,
                            modularity,
                            p,
                            moves);
}

} // namespace

int main()
{
    TestUnweighted();
    TestWeighted();
    TestSelfLoop();
    TestSingletonSource();
    TestAggregateNodeSize();
    TestNewTargetCommunity();
    TestActualMoveSequence();

    std::cout << "All quality-function tests passed.\n";
    return EXIT_SUCCESS;
}
