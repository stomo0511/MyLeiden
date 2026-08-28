#include "QualityFunction.hpp"

#ifdef ENABLE_MOVENODESFAST_PROFILE
#include <chrono>
#endif
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

#ifdef ENABLE_MOVENODESFAST_PROFILE
// TEMPORARY MOVENODESFAST PERFORMANCE PROFILING
using Clock = std::chrono::steady_clock;

double ElapsedSeconds(Clock::time_point begin, Clock::time_point end)
{
    return std::chrono::duration<double>(end - begin).count();
}
// END TEMPORARY MOVENODESFAST PERFORMANCE PROFILING
#endif

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
        partition.internal_edge_weight.size() != nc ||
        partition.community_is_empty.size() != nc) {
        throw std::invalid_argument("community statistic arrays have inconsistent sizes");
    }
    if (partition.smallest_empty_community < 0 ||
        static_cast<std::size_t>(partition.smallest_empty_community) > nc) {
        throw std::invalid_argument("smallest empty community is out of range");
    }
}

void ValidateStatsPartitionVertex(const LeidenGraphStats& stats,
                                  const LeidenPartition& partition,
                                  Vertex v)
{
    if (v < 0 || static_cast<std::size_t>(v) >= stats.node_size.size() ||
        static_cast<std::size_t>(v) >= stats.node_strength.size() ||
        static_cast<std::size_t>(v) >= partition.community_of.size()) {
        throw std::out_of_range("vertex id out of range");
    }
    const std::size_t nc = partition.community_size.size();
    if (partition.community_strength.size() != nc ||
        partition.internal_edge_weight.size() != nc ||
        partition.community_is_empty.size() != nc) {
        throw std::invalid_argument("community statistic arrays have inconsistent sizes");
    }
    if (partition.smallest_empty_community < 0 ||
        static_cast<std::size_t>(partition.smallest_empty_community) > nc) {
        throw std::invalid_argument("smallest empty community is out of range");
    }
    const Community source = partition.community_of[v];
    if (source < 0 || static_cast<std::size_t>(source) >= nc) {
        throw std::invalid_argument("node is not assigned to a valid community");
    }
}

void RefreshSmallestEmptyCommunity(LeidenPartition& partition)
{
    const Community nc = static_cast<Community>(partition.community_is_empty.size());
    while (partition.smallest_empty_community < nc &&
           !partition.community_is_empty[partition.smallest_empty_community]) {
        ++partition.smallest_empty_community;
    }
}

void MarkCommunityEmpty(LeidenPartition& partition, Community community)
{
    partition.community_is_empty[community] = 1;
    if (community < partition.smallest_empty_community) {
        partition.smallest_empty_community = community;
    }
}

void MarkCommunityNonEmpty(LeidenPartition& partition, Community community)
{
    partition.community_is_empty[community] = 0;
    if (community == partition.smallest_empty_community) {
        RefreshSmallestEmptyCommunity(partition);
    }
}

void EnsureCommunity(LeidenPartition& partition, Community community)
{
    if (community < 0) {
        throw std::out_of_range("community id out of range");
    }
    const std::size_t old_size = partition.community_size.size();
    const std::size_t need = static_cast<std::size_t>(community) + 1;
    if (old_size < need) {
        partition.community_size.resize(need, 0.0);
        partition.community_strength.resize(need, 0.0);
        partition.internal_edge_weight.resize(need, 0.0);
        partition.community_is_empty.resize(need, 1);
        if (static_cast<Community>(old_size) <
            partition.smallest_empty_community) {
            partition.smallest_empty_community =
                static_cast<Community>(old_size);
        }
    }
}

// double SumInternalEdgeWeight(const Graph& G,
//                              const LeidenPartition& partition,
//                              Community community)
// {
//     double sum = 0.0;
//     for_each_undirected_edge(G, [&](int u, int v, double w) {
//         if (partition.community_of[u] == community &&
//             partition.community_of[v] == community) {
//             sum += w;
//         }
//     });
//     return sum;
// }

} // namespace

std::unordered_map<Community, double>
BuildNeighborCommunityWeights(const Graph& G,
                              const LeidenPartition& partition,
                              Vertex v)
{
    ValidateVertex(v, G);
    ValidatePartitionSize(G, partition);

    std::unordered_map<Community, double> weights;
    weights.reserve(G.adj[v].size());
    for (const Edge& e : G.adj[v]) {
        if (e.to == v) {
            continue;
        }
        const Community community = partition.community_of[e.to];
        if (community >= 0) {
            weights[community] += e.weight;
        }
    }
    return weights;
}

void BuildNeighborCommunityWeights(const Graph& G,
                                   const LeidenPartition& partition,
                                   Vertex v,
                                   NeighborCommunityScratch& scratch)
{
    ValidateVertex(v, G);
    ValidatePartitionSize(G, partition);

    const std::size_t community_count = partition.community_size.size();
    if (scratch.weights.size() < community_count) {
        scratch.weights.resize(community_count, 0.0);
        scratch.marks.resize(community_count, 0);
    }

    if (scratch.generation == std::numeric_limits<std::size_t>::max()) {
        std::fill(scratch.marks.begin(), scratch.marks.end(), 0);
        scratch.generation = 1;
    } else {
        ++scratch.generation;
    }
    scratch.touched.clear();

    for (const Edge& e : G.adj[v]) {
        if (e.to == v) {
            continue;
        }
        const Community community = partition.community_of[e.to];
        if (community < 0) {
            continue;
        }

        const std::size_t c = static_cast<std::size_t>(community);
        if (scratch.marks[c] != scratch.generation) {
            scratch.marks[c] = scratch.generation;
            scratch.weights[c] = 0.0;
            scratch.touched.push_back(community);
        }
        scratch.weights[c] += e.weight;
    }
}

double LookupNeighborCommunityWeight(
    const NeighborCommunityScratch& scratch,
    Community community)
{
    if (community < 0 ||
        static_cast<std::size_t>(community) >= scratch.weights.size()) {
        return 0.0;
    }
    const std::size_t c = static_cast<std::size_t>(community);
    return (scratch.marks[c] == scratch.generation)
               ? scratch.weights[c]
               : 0.0;
}

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
    partition.community_is_empty.assign(nc, 0);
    partition.smallest_empty_community = nc;

    for (int v = 0; v < num_vertices(G); ++v) {
        const Community c = community_of[v];
        if (c < 0) {
            throw std::invalid_argument("community_of contains negative id");
        }
        partition.community_size[c] += stats.node_size[v];
        partition.community_strength[c] += stats.node_strength[v];
    }

    for (Community c = 0; c < nc; ++c) {
        if (partition.community_size[c] == 0.0) {
            MarkCommunityEmpty(partition, c);
        }
    }

    // for (Community c = 0; c < nc; ++c) {
    //     partition.internal_edge_weight[c] = SumInternalEdgeWeight(G, partition, c);
    // }
    for_each_undirected_edge(G, [&](int u, int v, double w) {
        const Community cu = partition.community_of[u];
        const Community cv = partition.community_of[v];

        if (cu == cv) {
            partition.internal_edge_weight[cu] += w;
        }
    }); 

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
    if (partition.community_size[source] == 0.0) {
        MarkCommunityEmpty(partition, source);
    }
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
    MarkCommunityNonEmpty(partition, community);
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

void MoveNodeToCommunityFromWeights(const Graph& G,
                                    const LeidenGraphStats& stats,
                                    LeidenPartition& partition,
                                    Vertex v,
                                    Community community,
                                    double weight_to_source,
                                    double weight_to_target,
                                    double self_loop_weight,
                                    MoveNodesFastProfile* profile)
{
#ifdef ENABLE_MOVENODESFAST_PROFILE
    // TEMPORARY MOVENODESFAST PERFORMANCE PROFILING
    const Clock::time_point validation_begin = Clock::now();
#else
    (void)profile;
#endif
    ValidateStatsSize(G, stats);
    ValidatePartitionSize(G, partition);
    ValidateVertex(v, G);
#ifdef ENABLE_MOVENODESFAST_PROFILE
    // TEMPORARY MOVENODESFAST PERFORMANCE PROFILING
    if (profile != nullptr) {
        profile->move_validation +=
            ElapsedSeconds(validation_begin, Clock::now());
    }
#endif

    const Community source = partition.community_of[v];
    if (source < 0 ||
        static_cast<std::size_t>(source) >= partition.community_size.size()) {
        throw std::invalid_argument("node is not assigned to a valid community");
    }
    if (source == community) {
        return;
    }

#ifdef ENABLE_MOVENODESFAST_PROFILE
    // TEMPORARY MOVENODESFAST PERFORMANCE PROFILING
    const Clock::time_point ensure_begin = Clock::now();
#endif
    EnsureCommunity(partition, community);
#ifdef ENABLE_MOVENODESFAST_PROFILE
    // TEMPORARY MOVENODESFAST PERFORMANCE PROFILING
    if (profile != nullptr) {
        profile->move_ensure_community +=
            ElapsedSeconds(ensure_begin, Clock::now());
    }
#endif

#ifdef ENABLE_MOVENODESFAST_PROFILE
    // TEMPORARY MOVENODESFAST PERFORMANCE PROFILING
    const Clock::time_point statistics_begin = Clock::now();
#endif
    const double source_internal_incident =
        weight_to_source + self_loop_weight;
    const double target_internal_incident =
        weight_to_target + self_loop_weight;

    partition.community_size[source] -= stats.node_size[v];
    partition.community_strength[source] -= stats.node_strength[v];
    partition.internal_edge_weight[source] -= source_internal_incident;
    partition.community_of[v] = -1;
#ifdef ENABLE_MOVENODESFAST_PROFILE
    // TEMPORARY MOVENODESFAST PERFORMANCE PROFILING
    const Clock::time_point statistics_first_end = Clock::now();
    if (profile != nullptr) {
        profile->move_statistics_update +=
            ElapsedSeconds(statistics_begin, statistics_first_end);
    }
    const Clock::time_point empty_first_begin = statistics_first_end;
#endif

    if (partition.community_size[source] == 0.0) {
        MarkCommunityEmpty(partition, source);
    }
#ifdef ENABLE_MOVENODESFAST_PROFILE
    // TEMPORARY MOVENODESFAST PERFORMANCE PROFILING
    const Clock::time_point empty_first_end = Clock::now();
    if (profile != nullptr) {
        profile->move_empty_community +=
            ElapsedSeconds(empty_first_begin, empty_first_end);
    }
    const Clock::time_point statistics_second_begin = empty_first_end;
#endif

    partition.community_of[v] = community;
    partition.community_size[community] += stats.node_size[v];
    partition.community_strength[community] += stats.node_strength[v];
    partition.internal_edge_weight[community] += target_internal_incident;
#ifdef ENABLE_MOVENODESFAST_PROFILE
    // TEMPORARY MOVENODESFAST PERFORMANCE PROFILING
    const Clock::time_point statistics_second_end = Clock::now();
    if (profile != nullptr) {
        profile->move_statistics_update +=
            ElapsedSeconds(statistics_second_begin, statistics_second_end);
    }
    const Clock::time_point empty_second_begin = statistics_second_end;
#endif

    MarkCommunityNonEmpty(partition, community);
#ifdef ENABLE_MOVENODESFAST_PROFILE
    // TEMPORARY MOVENODESFAST PERFORMANCE PROFILING
    if (profile != nullptr) {
        profile->move_empty_community +=
            ElapsedSeconds(empty_second_begin, Clock::now());
    }
#endif
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

    const double w_to_source =
        WeightFromNodeToCommunity(G, partition, v, source, false);
    const double w_to_target =
        WeightFromNodeToCommunity(G, partition, v, target_community, false);

    return deltaMoveFromWeights(stats,
                                partition,
                                v,
                                target_community,
                                w_to_source,
                                w_to_target);
}

double CPMQualityFunction::deltaMoveFromWeights(const LeidenGraphStats& stats,
                                                const LeidenPartition& partition,
                                                Vertex v,
                                                Community target_community,
                                                double weight_to_source,
                                                double weight_to_target) const
{
    ValidateStatsPartitionVertex(stats, partition, v);
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

    return (weight_to_target - weight_to_source)
           - gamma_ * stats.node_size[v] * (target_size - source_size_without_v);
}

double CPMQualityFunction::refinementNodeMass(const LeidenGraphStats& stats,
                                              Vertex v) const
{
    if (v < 0 || static_cast<std::size_t>(v) >= stats.node_size.size()) {
        throw std::out_of_range("vertex id out of range");
    }
    return stats.node_size[v];
}

double CPMQualityFunction::refinementResolution(const LeidenGraphStats& stats) const
{
    (void)stats;
    return gamma_;
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

    const double w_to_source =
        WeightFromNodeToCommunity(G, partition, v, source, false);
    const double w_to_target =
        WeightFromNodeToCommunity(G, partition, v, target_community, false);

    return deltaMoveFromWeights(stats,
                                partition,
                                v,
                                target_community,
                                w_to_source,
                                w_to_target);
}

double ModularityQualityFunction::deltaMoveFromWeights(
    const LeidenGraphStats& stats,
    const LeidenPartition& partition,
    Vertex v,
    Community target_community,
    double weight_to_source,
    double weight_to_target) const
{
    ValidateStatsPartitionVertex(stats, partition, v);
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

    return 2.0 * (weight_to_target - weight_to_source)
           - gamma_ * kv * (target_strength - source_strength_without_v)
                 / stats.total_edge_weight;
}

double ModularityQualityFunction::refinementNodeMass(
    const LeidenGraphStats& stats,
    Vertex v) const
{
    if (v < 0 || static_cast<std::size_t>(v) >= stats.node_strength.size()) {
        throw std::out_of_range("vertex id out of range");
    }
    return stats.node_strength[v];
}

double ModularityQualityFunction::refinementResolution(
    const LeidenGraphStats& stats) const
{
    if (stats.total_edge_weight == 0.0) {
        return 0.0;
    }
    return gamma_ / (2.0 * stats.total_edge_weight);
}
