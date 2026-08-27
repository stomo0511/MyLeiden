#include "QualityFunction.hpp"

#include <cmath>
#include <stdexcept>

namespace {

void ValidateStatsSize(const Graph& G, const LeidenGraphStats& stats)
{
    const int n = num_vertices(G);
    if (static_cast<int>(stats.node_size.size()) != n ||
        static_cast<int>(stats.node_strength.size()) != n) {
        throw std::invalid_argument("graph stats size does not match graph");
    }
}

void ValidatePartitionSize(const Graph& G, const LeidenPartition& partition)
{
    if (static_cast<int>(partition.community_of.size()) != num_vertices(G)) {
        throw std::invalid_argument("partition size does not match graph");
    }
    const std::size_t nc = partition.community_size.size();
    if (partition.community_strength.size() != nc ||
        partition.internal_edge_weight.size() != nc) {
        throw std::invalid_argument("community statistic arrays have inconsistent sizes");
    }
}

void EnsureCommunity(LeidenPartition& partition, Community community)
{
    if (community < 0) {
        throw std::out_of_range("community id out of range");
    }
    const std::size_t need = static_cast<std::size_t>(community) + 1;
    if (partition.community_size.size() < need) {
        partition.community_size.resize(need, 0.0);
        partition.community_strength.resize(need, 0.0);
        partition.internal_edge_weight.resize(need, 0.0);
    }
}

double SumInternalEdgeWeight(const Graph& G,
                             const LeidenPartition& partition,
                             Community community)
{
    double sum = 0.0;
    for_each_undirected_edge(G, [&](int u, int v, double w) {
        if (partition.community_of[u] == community &&
            partition.community_of[v] == community) {
            sum += w;
        }
    });
    return sum;
}

} // namespace

LeidenGraphStats BuildLeidenGraphStats(const Graph& G)
{
    std::vector<double> node_size(num_vertices(G), 1.0);
    return BuildLeidenGraphStats(G, node_size);
}

LeidenGraphStats BuildLeidenGraphStats(const Graph& G,
                                       const std::vector<double>& node_size)
{
    const int n = num_vertices(G);
    if (static_cast<int>(node_size.size()) != n) {
        throw std::invalid_argument("node_size length does not match graph");
    }

    LeidenGraphStats stats;
    stats.node_size = node_size;
    stats.node_strength.assign(n, 0.0);

    for_each_undirected_edge(G, [&](int u, int v, double w) {
        stats.total_edge_weight += w;
        if (u == v) {
            stats.node_strength[u] += 2.0 * w;
        } else {
            stats.node_strength[u] += w;
            stats.node_strength[v] += w;
        }
    });

    return stats;
}

LeidenPartition MakeSingletonPartition(const Graph& G,
                                       const LeidenGraphStats& stats)
{
    ValidateStatsSize(G, stats);

    const int n = num_vertices(G);
    std::vector<Community> community_of(n);
    for (int v = 0; v < n; ++v) {
        community_of[v] = v;
    }
    return MakePartition(G, stats, community_of);
}

LeidenPartition MakePartition(const Graph& G,
                              const LeidenGraphStats& stats,
                              const std::vector<Community>& community_of)
{
    ValidateStatsSize(G, stats);
    if (static_cast<int>(community_of.size()) != num_vertices(G)) {
        throw std::invalid_argument("community_of length does not match graph");
    }

    const int nc = NumCommunitiesFromAssignment(community_of);
    LeidenPartition partition;
    partition.community_of = community_of;
    partition.community_size.assign(nc, 0.0);
    partition.community_strength.assign(nc, 0.0);
    partition.internal_edge_weight.assign(nc, 0.0);

    for (int v = 0; v < num_vertices(G); ++v) {
        const Community c = community_of[v];
        if (c < 0) {
            throw std::invalid_argument("community_of contains negative id");
        }
        partition.community_size[c] += stats.node_size[v];
        partition.community_strength[c] += stats.node_strength[v];
    }

    for (Community c = 0; c < nc; ++c) {
        partition.internal_edge_weight[c] = SumInternalEdgeWeight(G, partition, c);
    }

    return partition;
}

double WeightFromNodeToCommunity(const Graph& G,
                                 const LeidenPartition& partition,
                                 Vertex v,
                                 Community community,
                                 bool include_self_loop)
{
    ValidateVertex(v, G);
    ValidatePartitionSize(G, partition);
    if (community < 0) {
        throw std::out_of_range("community id out of range");
    }

    double weight = 0.0;
    for (const Edge& e : G.adj[v]) {
        if (e.to == v && !include_self_loop) {
            continue;
        }
        if (partition.community_of[e.to] == community) {
            weight += e.weight;
        }
    }
    return weight;
}

double SelfLoopWeight(const Graph& G, Vertex v)
{
    ValidateVertex(v, G);
    double weight = 0.0;
    for (const Edge& e : G.adj[v]) {
        if (e.to == v) {
            weight += e.weight;
        }
    }
    return weight;
}

void RemoveNodeFromCommunity(const Graph& G,
                             const LeidenGraphStats& stats,
                             LeidenPartition& partition,
                             Vertex v)
{
    ValidateStatsSize(G, stats);
    ValidatePartitionSize(G, partition);
    ValidateVertex(v, G);

    const Community source = partition.community_of[v];
    if (source < 0 || static_cast<std::size_t>(source) >= partition.community_size.size()) {
        throw std::invalid_argument("node is not assigned to a valid community");
    }

    const double internal_incident =
        WeightFromNodeToCommunity(G, partition, v, source, true);

    partition.community_size[source] -= stats.node_size[v];
    partition.community_strength[source] -= stats.node_strength[v];
    partition.internal_edge_weight[source] -= internal_incident;
    partition.community_of[v] = -1;
}

void InsertNodeIntoCommunity(const Graph& G,
                             const LeidenGraphStats& stats,
                             LeidenPartition& partition,
                             Vertex v,
                             Community community)
{
    ValidateStatsSize(G, stats);
    ValidatePartitionSize(G, partition);
    ValidateVertex(v, G);
    if (partition.community_of[v] != -1) {
        throw std::invalid_argument("node is already assigned to a community");
    }

    EnsureCommunity(partition, community);
    const double internal_incident =
        WeightFromNodeToCommunity(G, partition, v, community, false)
        + SelfLoopWeight(G, v);

    partition.community_of[v] = community;
    partition.community_size[community] += stats.node_size[v];
    partition.community_strength[community] += stats.node_strength[v];
    partition.internal_edge_weight[community] += internal_incident;
}

void MoveNodeToCommunity(const Graph& G,
                         const LeidenGraphStats& stats,
                         LeidenPartition& partition,
                         Vertex v,
                         Community community)
{
    ValidatePartitionSize(G, partition);
    ValidateVertex(v, G);
    if (partition.community_of[v] == community) {
        return;
    }
    RemoveNodeFromCommunity(G, stats, partition, v);
    InsertNodeIntoCommunity(G, stats, partition, v, community);
}

CPMQualityFunction::CPMQualityFunction(double gamma)
    : gamma_(gamma)
{
}

double CPMQualityFunction::quality(const Graph& G,
                                   const LeidenGraphStats& stats,
                                   const LeidenPartition& partition) const
{
    (void)G;
    (void)stats;
    double h = 0.0;

    // CPM of Traag et al. is sum_C [m_C - gamma * choose(n_C, 2)].
    // m_C is the total weight of undirected internal edges, including
    // self-loops once. n_C is the sum of aggregate node sizes, so the
    // null-model term still counts pairs of original nodes after aggregation.
    for (std::size_t c = 0; c < partition.community_size.size(); ++c) {
        const double size = partition.community_size[c];
        h += partition.internal_edge_weight[c]
             - gamma_ * size * (size - 1.0) / 2.0;
    }
    return h;
}

double CPMQualityFunction::deltaMove(const Graph& G,
                                     const LeidenGraphStats& stats,
                                     const LeidenPartition& partition,
                                     Vertex v,
                                     Community target_community) const
{
    ValidateStatsSize(G, stats);
    ValidatePartitionSize(G, partition);
    ValidateVertex(v, G);
    if (target_community < 0) {
        throw std::out_of_range("target community id out of range");
    }

    const Community source = partition.community_of[v];
    if (source == target_community) {
        return 0.0;
    }

    const double source_size_without_v =
        partition.community_size[source] - stats.node_size[v];
    const double target_size =
        (static_cast<std::size_t>(target_community) < partition.community_size.size())
            ? partition.community_size[target_community]
            : 0.0;
    const double w_to_source =
        WeightFromNodeToCommunity(G, partition, v, source, false);
    const double w_to_target =
        WeightFromNodeToCommunity(G, partition, v, target_community, false);

    return (w_to_target - w_to_source)
           - gamma_ * stats.node_size[v] * (target_size - source_size_without_v);
}

ModularityQualityFunction::ModularityQualityFunction(double gamma)
    : gamma_(gamma)
{
}

double ModularityQualityFunction::quality(const Graph& G,
                                          const LeidenGraphStats& stats,
                                          const LeidenPartition& partition) const
{
    (void)G;
    if (stats.total_edge_weight == 0.0) {
        return 0.0;
    }

    double h = 0.0;

    // This uses H = 2m * Q_gamma, where
    // Q_gamma = sum_C [m_C / m - gamma * (K_C / (2m))^2].
    // Here m is total undirected edge weight, m_C counts internal undirected
    // edge weight including self-loops once, and K_C is total vertex strength
    // with each self-loop contributing 2w. Since 2m is a positive constant,
    // deltaMove() has the same sign as the usual modularity change.
    for (std::size_t c = 0; c < partition.community_strength.size(); ++c) {
        const double strength = partition.community_strength[c];
        h += 2.0 * partition.internal_edge_weight[c]
             - gamma_ * strength * strength / (2.0 * stats.total_edge_weight);
    }
    return h;
}

double ModularityQualityFunction::deltaMove(const Graph& G,
                                            const LeidenGraphStats& stats,
                                            const LeidenPartition& partition,
                                            Vertex v,
                                            Community target_community) const
{
    ValidateStatsSize(G, stats);
    ValidatePartitionSize(G, partition);
    ValidateVertex(v, G);
    if (target_community < 0) {
        throw std::out_of_range("target community id out of range");
    }
    if (stats.total_edge_weight == 0.0) {
        return 0.0;
    }

    const Community source = partition.community_of[v];
    if (source == target_community) {
        return 0.0;
    }

    const double kv = stats.node_strength[v];
    const double source_strength_without_v =
        partition.community_strength[source] - kv;
    const double target_strength =
        (static_cast<std::size_t>(target_community) < partition.community_strength.size())
            ? partition.community_strength[target_community]
            : 0.0;
    const double w_to_source =
        WeightFromNodeToCommunity(G, partition, v, source, false);
    const double w_to_target =
        WeightFromNodeToCommunity(G, partition, v, target_community, false);

    return 2.0 * (w_to_target - w_to_source)
           - gamma_ * kv * (target_strength - source_strength_without_v)
                 / stats.total_edge_weight;
}
