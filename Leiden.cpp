#include "Leiden.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <deque>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

struct AggregateScratch {
    std::vector<double> weights;
    std::vector<std::size_t> marks;
    std::vector<Community> touched;
    std::size_t generation = 1;

    explicit AggregateScratch(int n_coarse)
        : weights(n_coarse, 0.0),
          marks(n_coarse, 0)
    {
    }

    void reset()
    {
        touched.clear();
        if (generation == std::numeric_limits<std::size_t>::max()) {
            std::fill(marks.begin(), marks.end(), 0);
            generation = 1;
        } else {
            ++generation;
        }
    }

    void add(Community community, double weight)
    {
        const std::size_t c = static_cast<std::size_t>(community);
        if (marks[c] != generation) {
            marks[c] = generation;
            weights[c] = 0.0;
            touched.push_back(community);
        }
        weights[c] += weight;
    }
};

struct AggregateCommunityEdges {
    bool has_self_loop = false;
    double self_loop_weight = 0.0;
    std::vector<Edge> upper_edges;
};

struct LocalRefinementResult {
    std::vector<Vertex> vertices;
    std::vector<Community> local_assignment;
    Community local_community_count = 0;
};

#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
// TEMPORARY REFINEPARTITION DETAILED PROFILING
struct RefineSubsetProfile {
    Community parent = -1;
    std::size_t vertices = 0;
    std::size_t adjacency_volume = 0;
    std::size_t candidates = 0;
    std::size_t processed_candidates = 0;
    std::size_t singleton_source_skips = 0;
    std::size_t neighbor_builds = 0;
    std::size_t neighbor_adjacency_scans = 0;
    std::size_t touched_communities = 0;
    std::size_t active_community_iterations = 0;
    std::size_t delta_evaluations = 0;
    std::size_t nondecreasing_candidates = 0;
    std::size_t stochastic_selections = 0;
    std::size_t refinement_moves = 0;
    std::size_t stats_build_vertices = 0;
    std::size_t stats_build_adjacency_scans = 0;
    std::size_t stats_update_calls = 0;
    std::size_t stats_update_adjacency_scans = 0;
    Community local_community_count = 0;
    double reset_time = 0.0;
    double subset_mass_time = 0.0;
    double candidate_scan_time = 0.0;
    double shuffle_time = 0.0;
    double scratch_init_time = 0.0;
    double stats_build_time = 0.0;
    double candidate_processing_time = 0.0;
    double result_compaction_time = 0.0;
    double total_time = 0.0;
};

struct RefineThreadProfile {
    std::size_t subsets = 0;
    std::size_t vertices = 0;
    std::size_t adjacency_volume = 0;
    double initialization_time = 0.0;
    double subset_work_time = 0.0;
};
#endif

struct MoveProposal {
    bool active = false;
    bool positive = false;
    bool target_was_empty = false;
    Vertex vertex = -1;
    Community source = -1;
    Community target = -1;
    double delta = 0.0;
};

std::uint64_t MixStage4A(std::uint64_t value)
{
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
}

std::uint64_t Stage4APriority(std::uint64_t seed,
                              std::uint64_t level,
                              Vertex vertex)
{
    std::uint64_t value = MixStage4A(seed);
    value ^= MixStage4A(level + 0x632be59bd9b4e019ULL);
    value ^= MixStage4A(static_cast<std::uint64_t>(vertex) +
                        0x85157af5ULL);
    return MixStage4A(value);
}

using Clock = std::chrono::steady_clock;

double ElapsedSeconds(Clock::time_point begin, Clock::time_point end)
{
    return std::chrono::duration<double>(end - begin).count();
}

bool DebugEnabled(const LeidenOptions* options)
{
    return options != nullptr && options->debug;
}

std::size_t DebugInterval(const LeidenOptions* options)
{
    if (options == nullptr || options->debug_interval == 0) {
        return 100000;
    }
    return options->debug_interval;
}

std::size_t CountActiveCommunities(const LeidenPartition& partition)
{
    std::vector<Community> active = partition.community_of;
    std::sort(active.begin(), active.end());
    active.erase(std::unique(active.begin(), active.end()), active.end());
    return active.size();
}

double LookupWeight(const std::unordered_map<Community, double>& weights,
                    Community community)
{
    const auto it = weights.find(community);
    return (it == weights.end()) ? 0.0 : it->second;
}

double SubsetMass(const LeidenGraphStats& stats,
                  const QualityFunction& quality_function,
                  const std::vector<Vertex>& subset)
{
    double mass = 0.0;
    for (Vertex v : subset) {
        mass += quality_function.refinementNodeMass(stats, v);
    }
    return mass;
}

std::mt19937 MakeSubsetRng(std::uint64_t global_seed,
                           std::uint64_t level,
                           std::uint64_t parent)
{
    const std::uint64_t seed = MakeSubsetSeed(global_seed, level, parent);
    std::seed_seq seq{
        static_cast<std::uint32_t>(seed),
        static_cast<std::uint32_t>(seed >> 32U),
        static_cast<std::uint32_t>(parent),
        static_cast<std::uint32_t>(parent >> 32U),
        static_cast<std::uint32_t>(level),
        static_cast<std::uint32_t>(level >> 32U)};
    return std::mt19937(seq);
}

void ResetLocalSingletonsForSubset(const Graph& G,
                                   const LeidenGraphStats& stats,
                                   const std::vector<double>& singleton_internal_edge_weight,
                                   LeidenPartition& local_partition,
                                   const std::vector<Vertex>& subset)
{
    for (Vertex v : subset) {
        ValidateVertex(v, G);
        local_partition.community_of[v] = v;
        local_partition.community_size[v] = stats.node_size[v];
        local_partition.community_strength[v] = stats.node_strength[v];
        local_partition.internal_edge_weight[v] =
            singleton_internal_edge_weight[v];
        local_partition.community_is_empty[v] = 0;
    }
    local_partition.smallest_empty_community =
        static_cast<Community>(local_partition.community_size.size());
}

LeidenPartition MakeThreadLocalSingletonPartition(const Graph& G,
                                                  const LeidenGraphStats& stats,
                                                  const std::vector<double>& singleton_internal_edge_weight)
{
    const int n = num_vertices(G);
    LeidenPartition partition;
    partition.community_of.resize(n);
    std::iota(partition.community_of.begin(), partition.community_of.end(), 0);
    partition.community_size = stats.node_size;
    partition.community_strength = stats.node_strength;
    partition.internal_edge_weight = singleton_internal_edge_weight;
    partition.community_is_empty.assign(n, 0);
    partition.smallest_empty_community = n;
    (void)G;
    return partition;
}

double EdgeWeightFromNodeToParentSubset(const Graph& G,
                                        const LeidenPartition& partition,
                                        Vertex v,
                                        Community parent)
{
    ValidateVertex(v, G);
    double weight = 0.0;
    for (const Edge& e : G.adj[v]) {
        if (e.to != v && partition.community_of[e.to] == parent) {
            weight += e.weight;
        }
    }
    return weight;
}

bool IsNodeWellConnectedToParentSubset(
    const Graph& G,
    const LeidenGraphStats& stats,
    const LeidenPartition& partition,
    const QualityFunction& quality_function,
    Vertex v,
    Community parent,
    double subset_mass)
{
    const double node_mass = quality_function.refinementNodeMass(stats, v);
    const double edge_weight =
        EdgeWeightFromNodeToParentSubset(G, partition, v, parent);
    const double threshold =
        quality_function.refinementResolution(stats) *
        node_mass *
        (subset_mass - node_mass);
    return edge_weight >= threshold;
}

RefinementCommunityEntry& EnsureRefinementCommunityEntry(
    RefinementCommunityStats& community_stats,
    Community community)
{
    if (community < 0) {
        throw std::out_of_range("community id out of range");
    }
    return community_stats.entries[community];
}

void RefreshActiveCommunities(RefinementCommunityStats& community_stats)
{
    community_stats.active_communities.clear();
    community_stats.active_communities.reserve(community_stats.entries.size());
    for (const auto& item : community_stats.entries) {
        if (item.second.member_count > 0) {
            community_stats.active_communities.push_back(item.first);
        }
    }
    std::sort(community_stats.active_communities.begin(),
              community_stats.active_communities.end());
}

std::vector<std::vector<Vertex>> BuildPartitionMembers(
    const LeidenPartition& partition)
{
    std::vector<std::vector<Vertex>> members(partition.community_size.size());
    for (Vertex v = 0; v < static_cast<Vertex>(partition.community_of.size()); ++v) {
        const Community community = partition.community_of[v];
        if (community >= 0) {
            if (static_cast<std::size_t>(community) >= members.size()) {
                members.resize(static_cast<std::size_t>(community) + 1);
            }
            members[community].push_back(v);
        }
    }
    return members;
}

std::vector<Community> BuildCandidateCommunities(
    const LeidenPartition& partition,
    const NeighborCommunityScratch& scratch)
{
    std::vector<Community> candidates;
    candidates.reserve(scratch.touched.size() + 1);

    for (Community community : scratch.touched) {
        if (community >= 0) {
            candidates.push_back(community);
        }
    }

    candidates.push_back(EmptyCommunityForMove(partition));
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()),
                     candidates.end());
    return candidates;
}

void ValidateAggregateInput(const Graph& G,
                            const LeidenGraphStats& stats,
                            const LeidenPartition& refined)
{
    const int n = num_vertices(G);
    if (static_cast<int>(stats.node_size.size()) != n ||
        static_cast<int>(stats.node_strength.size()) != n ||
        static_cast<int>(refined.community_of.size()) != n) {
        throw std::invalid_argument("aggregate input size does not match graph");
    }
    for (Community community : refined.community_of) {
        if (community < 0) {
            throw std::invalid_argument("refined partition contains negative community id");
        }
    }
}

void ValidateLeidenInput(const Graph& G,
                         const LeidenGraphStats& stats,
                         const LeidenOptions& options)
{
    const int n = num_vertices(G);
    if (static_cast<int>(stats.node_size.size()) != n ||
        static_cast<int>(stats.node_strength.size()) != n) {
        throw std::invalid_argument("Leiden input stats size does not match graph");
    }
    if (options.theta <= 0.0) {
        throw std::invalid_argument("theta must be positive");
    }
    if (options.debug_interval == 0) {
        throw std::invalid_argument("debug_interval must be positive");
    }
}

std::vector<Community> BuildCommunityToCoarse(const LeidenPartition& refined)
{
    std::vector<Community> active;
    active.reserve(refined.community_size.size());
    for (Community c : refined.community_of) {
        active.push_back(c);
    }
    std::sort(active.begin(), active.end());
    active.erase(std::unique(active.begin(), active.end()), active.end());

    const Community max_community = active.empty() ? -1 : active.back();
    std::vector<Community> community_to_coarse(
        static_cast<std::size_t>(max_community + 1),
        -1);
    for (Community coarse = 0; coarse < static_cast<Community>(active.size()); ++coarse) {
        community_to_coarse[active[coarse]] = coarse;
    }
    return community_to_coarse;
}

std::vector<Community> BuildCompactCommunityMap(
    const std::vector<Community>& community_ids)
{
    std::vector<Community> active = community_ids;
    std::sort(active.begin(), active.end());
    active.erase(std::unique(active.begin(), active.end()), active.end());

    const Community max_community = active.empty() ? -1 : active.back();
    std::vector<Community> community_map(
        static_cast<std::size_t>(max_community + 1),
        -1);
    for (Community compact = 0;
         compact < static_cast<Community>(active.size());
         ++compact) {
        community_map[active[compact]] = compact;
    }
    return community_map;
}

std::vector<Community> CompactAssignment(
    const std::vector<Community>& assignment)
{
    if (assignment.empty()) {
        return {};
    }
    const std::vector<Community> compact_map =
        BuildCompactCommunityMap(assignment);

    std::vector<Community> compacted(assignment.size(), -1);
    for (std::size_t i = 0; i < assignment.size(); ++i) {
        const Community community = assignment[i];
        if (community < 0) {
            throw std::invalid_argument("assignment contains negative community id");
        }
        compacted[i] = compact_map[community];
    }
    return compacted;
}

} // namespace

Stage4A1ReactivationStats UpdateStage4A1AffectedNeighbors(
    const Graph& G,
    const LeidenPartition& committed_partition,
    Vertex moved_vertex,
    Community committed_target,
    std::vector<unsigned char>& affected_next)
{
    ValidateVertex(moved_vertex, G);
    if (affected_next.size() != static_cast<std::size_t>(num_vertices(G))) {
        throw std::invalid_argument("affected_next size does not match graph");
    }

    Stage4A1ReactivationStats stats;
    for (const Edge& e : G.adj[moved_vertex]) {
        const Vertex u = e.to;
        if (u == moved_vertex) {
            continue;
        }
        ++stats.neighbor_scans;
        if (committed_partition.community_of[u] == committed_target) {
            ++stats.target_community_exclusions;
            continue;
        }
        unsigned char& mark = affected_next[static_cast<std::size_t>(u)];
        if (mark == 0) {
            mark = 1;
            ++stats.newly_activated;
        } else {
            ++stats.duplicate_attempts;
        }
    }
    return stats;
}

double EdgeWeightFromNodeToSubset(const Graph& G,
                                  Vertex v,
                                  const std::vector<std::size_t>& subset_mark,
                                  std::size_t subset_generation)
{
    ValidateVertex(v, G);
    if (subset_mark.size() != static_cast<std::size_t>(num_vertices(G))) {
        throw std::invalid_argument("subset mark size does not match graph");
    }

    double weight = 0.0;
    for (const Edge& e : G.adj[v]) {
        if (e.to != v && subset_mark[e.to] == subset_generation) {
            weight += e.weight;
        }
    }
    return weight;
}

bool IsNodeWellConnectedToSubset(const Graph& G,
                                 const LeidenGraphStats& stats,
                                 const QualityFunction& quality_function,
                                 Vertex v,
                                 const std::vector<std::size_t>& subset_mark,
                                 std::size_t subset_generation,
                                 double subset_mass)
{
    const double node_mass = quality_function.refinementNodeMass(stats, v);
    const double edge_weight =
        EdgeWeightFromNodeToSubset(G, v, subset_mark, subset_generation);
    const double threshold =
        quality_function.refinementResolution(stats) *
        node_mass *
        (subset_mass - node_mass);
    return edge_weight >= threshold;
}

RefinementCommunityStats BuildRefinementCommunityStats(
    const Graph& G,
    const LeidenGraphStats& stats,
    const LeidenPartition& refined,
    const QualityFunction& quality_function,
    const std::vector<Vertex>& subset,
    const std::vector<std::size_t>& subset_mark,
    std::size_t subset_generation)
{
    if (subset_mark.size() != static_cast<std::size_t>(num_vertices(G))) {
        throw std::invalid_argument("subset mark size does not match graph");
    }

    RefinementCommunityStats community_stats;
    community_stats.entries.reserve(subset.size());

    for (Vertex v : subset) {
        ValidateVertex(v, G);
        const Community community = refined.community_of[v];
        RefinementCommunityEntry& entry =
            EnsureRefinementCommunityEntry(community_stats, community);
        ++entry.member_count;
        entry.mass += quality_function.refinementNodeMass(stats, v);
    }

    for (Vertex u : subset) {
        ValidateVertex(u, G);
        const Community cu = refined.community_of[u];
        RefinementCommunityEntry& cu_entry =
            EnsureRefinementCommunityEntry(community_stats, cu);

        for (const Edge& e : G.adj[u]) {
            if (e.to == u ||
                subset_mark[e.to] != subset_generation ||
                u > e.to) {
                continue;
            }

            const Community cv = refined.community_of[e.to];
            RefinementCommunityEntry& cv_entry =
                EnsureRefinementCommunityEntry(community_stats, cv);
            if (cu != cv) {
                cu_entry.external_weight += e.weight;
                cv_entry.external_weight += e.weight;
            }
        }
    }

    RefreshActiveCommunities(community_stats);
    return community_stats;
}

bool IsCommunityWellConnectedFromStats(
    const LeidenGraphStats& stats,
    const QualityFunction& quality_function,
    Community community,
    double subset_mass,
    const RefinementCommunityStats& community_stats)
{
    if (community < 0) {
        return false;
    }
    const auto it = community_stats.entries.find(community);
    if (it == community_stats.entries.end() ||
        it->second.member_count <= 0) {
        return false;
    }

    const double community_mass = it->second.mass;
    const double threshold =
        quality_function.refinementResolution(stats) *
        community_mass *
        (subset_mass - community_mass);
    return it->second.external_weight >= threshold;
}

void UpdateRefinementCommunityStatsForMove(
    const Graph& G,
    const LeidenGraphStats& stats,
    const LeidenPartition& refined_before_move,
    const QualityFunction& quality_function,
    const std::vector<std::size_t>& subset_mark,
    std::size_t subset_generation,
    Vertex v,
    Community target,
    RefinementCommunityStats& community_stats)
{
    ValidateVertex(v, G);
    if (subset_mark.size() != static_cast<std::size_t>(num_vertices(G))) {
        throw std::invalid_argument("subset mark size does not match graph");
    }

    const Community source = refined_before_move.community_of[v];
    if (source == target) {
        return;
    }

    EnsureRefinementCommunityEntry(community_stats, source);
    EnsureRefinementCommunityEntry(community_stats, target);
    RefinementCommunityEntry& source_entry = community_stats.entries[source];
    RefinementCommunityEntry& target_entry = community_stats.entries[target];

    const double node_mass =
        quality_function.refinementNodeMass(stats, v);
    --source_entry.member_count;
    ++target_entry.member_count;
    source_entry.mass -= node_mass;
    target_entry.mass += node_mass;

    // For each subset-internal edge (v,u), update E(C,S-C) using communities
    // before the move A -> B:
    //   C == A: internal A edge becomes A-B external, so A += w and B += w.
    //   C == B: A-B external becomes internal B edge, so A -= w and B -= w.
    //   otherwise: A-C external becomes B-C external, so A -= w and B += w.
    // The third community C is external both before and after, so its value is
    // unchanged. Self-loops and subset-external edges are ignored.
    for (const Edge& e : G.adj[v]) {
        const Vertex u = e.to;
        if (u == v || subset_mark[u] != subset_generation) {
            continue;
        }

        const Community neighbor_community = refined_before_move.community_of[u];
        if (neighbor_community == source) {
            source_entry.external_weight += e.weight;
            target_entry.external_weight += e.weight;
        } else if (neighbor_community == target) {
            source_entry.external_weight -= e.weight;
            target_entry.external_weight -= e.weight;
        } else {
            source_entry.external_weight -= e.weight;
            target_entry.external_weight += e.weight;
        }
    }

    RefreshActiveCommunities(community_stats);
}

RefinementCommunityStats BuildRefinementCommunityStatsForParent(
    const Graph& G,
    const LeidenGraphStats& stats,
    const LeidenPartition& parent_partition,
    const LeidenPartition& local_refined,
    const QualityFunction& quality_function,
    const std::vector<Vertex>& subset,
    Community parent)
{
    RefinementCommunityStats community_stats;
    community_stats.entries.reserve(subset.size());

    for (Vertex v : subset) {
        ValidateVertex(v, G);
        const Community community = local_refined.community_of[v];
        RefinementCommunityEntry& entry =
            EnsureRefinementCommunityEntry(community_stats, community);
        ++entry.member_count;
        entry.mass += quality_function.refinementNodeMass(stats, v);
    }

    for (Vertex u : subset) {
        ValidateVertex(u, G);
        const Community cu = local_refined.community_of[u];
        RefinementCommunityEntry& cu_entry =
            EnsureRefinementCommunityEntry(community_stats, cu);

        for (const Edge& e : G.adj[u]) {
            if (e.to == u ||
                parent_partition.community_of[e.to] != parent ||
                u > e.to) {
                continue;
            }

            const Community cv = local_refined.community_of[e.to];
            RefinementCommunityEntry& cv_entry =
                EnsureRefinementCommunityEntry(community_stats, cv);
            if (cu != cv) {
                cu_entry.external_weight += e.weight;
                cv_entry.external_weight += e.weight;
            }
        }
    }

    RefreshActiveCommunities(community_stats);
    return community_stats;
}

void UpdateRefinementCommunityStatsForParentMove(
    const Graph& G,
    const LeidenGraphStats& stats,
    const LeidenPartition& parent_partition,
    const LeidenPartition& local_refined_before_move,
    const QualityFunction& quality_function,
    Vertex v,
    Community parent,
    Community target,
    RefinementCommunityStats& community_stats)
{
    ValidateVertex(v, G);
    const Community source = local_refined_before_move.community_of[v];
    if (source == target) {
        return;
    }

    EnsureRefinementCommunityEntry(community_stats, source);
    EnsureRefinementCommunityEntry(community_stats, target);
    RefinementCommunityEntry& source_entry = community_stats.entries[source];
    RefinementCommunityEntry& target_entry = community_stats.entries[target];

    const double node_mass =
        quality_function.refinementNodeMass(stats, v);
    --source_entry.member_count;
    ++target_entry.member_count;
    source_entry.mass -= node_mass;
    target_entry.mass += node_mass;

    for (const Edge& e : G.adj[v]) {
        const Vertex u = e.to;
        if (u == v || parent_partition.community_of[u] != parent) {
            continue;
        }

        const Community neighbor_community =
            local_refined_before_move.community_of[u];
        if (neighbor_community == source) {
            source_entry.external_weight += e.weight;
            target_entry.external_weight += e.weight;
        } else if (neighbor_community == target) {
            source_entry.external_weight -= e.weight;
            target_entry.external_weight -= e.weight;
        } else {
            source_entry.external_weight -= e.weight;
            target_entry.external_weight += e.weight;
        }
    }

    RefreshActiveCommunities(community_stats);
}

// -----------------------------------------------------------------------------
// TEMPORARY MOVENODESFAST PERFORMANCE PROFILING
// These local macros keep profiling removable and compile to nothing when the
// profile is disabled. Each measured region uses only one begin/end pair.
// -----------------------------------------------------------------------------
#ifdef ENABLE_MOVENODESFAST_PROFILE
#define MNF_PROFILE_BEGIN(name) \
    const Clock::time_point name##_begin = Clock::now()
#define MNF_PROFILE_END(name, field) \
    result.profile.field += ElapsedSeconds(name##_begin, Clock::now())
#define MNF_PROFILE_CANDIDATES(count) \
    result.profile.total_candidates += (count)
#define MNF_PROFILE_ADD(field, begin, end) \
    result.profile.field += ElapsedSeconds((begin), (end))
#define MNF_PROFILE_ARG (&result.profile)
#else
#define MNF_PROFILE_BEGIN(name)
#define MNF_PROFILE_END(name, field)
#define MNF_PROFILE_CANDIDATES(count)
#define MNF_PROFILE_ADD(field, begin, end)
#define MNF_PROFILE_ARG nullptr
#endif

MoveNodesFastResult MoveNodesFast(const Graph& G,
                                  const LeidenGraphStats& stats,
                                  LeidenPartition partition,
                                  const QualityFunction& quality_function,
                                  std::mt19937& rng,
                                  const LeidenOptions* options)
{
    const int n = num_vertices(G);
    std::vector<Vertex> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), rng);

    std::deque<Vertex> queue;
    std::vector<bool> in_queue(n, false);
    for (Vertex v : order) {
        queue.push_back(v);
        in_queue[v] = true;
    }

    MoveNodesFastResult result;
    result.partition = std::move(partition);
    NeighborCommunityScratch neighbor_scratch;
    neighbor_scratch.weights.assign(result.partition.community_size.size(), 0.0);
    neighbor_scratch.marks.assign(result.partition.community_size.size(), 0);
    neighbor_scratch.touched.reserve(32);
    const bool debug = DebugEnabled(options);
    const std::size_t debug_interval = DebugInterval(options);

    while (!queue.empty()) {
        const Vertex v = queue.front();
        queue.pop_front();
        in_queue[v] = false;
        ++result.num_visits;

        const Community source = result.partition.community_of[v];
        MNF_PROFILE_BEGIN(neighbor_weights);
        BuildNeighborCommunityWeights(
            G, result.partition, v, neighbor_scratch);
        MNF_PROFILE_END(neighbor_weights, neighbor_weights);
        const double weight_to_source =
            LookupNeighborCommunityWeight(neighbor_scratch, source);
        MNF_PROFILE_BEGIN(candidate_build);
        const std::vector<Community> candidates =
            BuildCandidateCommunities(result.partition, neighbor_scratch);
        MNF_PROFILE_END(candidate_build, candidate_build);
        MNF_PROFILE_CANDIDATES(candidates.size());

        Community best_community = source;
        double best_delta = 0.0;

        // Deterministic tie-breaking: if deltas are exactly equal, choose
        // the smaller community id. The initial node order remains randomized
        // by the caller-supplied rng.
        MNF_PROFILE_BEGIN(delta_evaluation);
        for (Community candidate : candidates) {
            const double delta =
                (candidate == source)
                    ? 0.0
                    : quality_function.deltaMoveFromWeights(
                          stats,
                          result.partition,
                          v,
                          candidate,
                          weight_to_source,
                          LookupNeighborCommunityWeight(neighbor_scratch,
                                                        candidate));

            if (delta > best_delta ||
                (delta == best_delta && candidate < best_community)) {
                best_delta = delta;
                best_community = candidate;
            }
        }
        MNF_PROFILE_END(delta_evaluation, delta_evaluation);

        if (best_delta > 0.0 && best_community != source) {
            MNF_PROFILE_BEGIN(move_node);
            const double weight_to_target =
                LookupNeighborCommunityWeight(neighbor_scratch,
                                              best_community);
#ifdef ENABLE_MOVENODESFAST_PROFILE
            // TEMPORARY MOVENODESFAST PERFORMANCE PROFILING
            const Clock::time_point self_loop_begin = Clock::now();
#endif
            const double self_loop_weight = SelfLoopWeight(G, v);
#ifdef ENABLE_MOVENODESFAST_PROFILE
            // TEMPORARY MOVENODESFAST PERFORMANCE PROFILING
            const Clock::time_point self_loop_end = Clock::now();
            MNF_PROFILE_ADD(move_self_loop_scan,
                            self_loop_begin,
                            self_loop_end);
#endif
            MoveNodeToCommunityFromWeights(G,
                                           stats,
                                           result.partition,
                                           v,
                                           best_community,
                                           weight_to_source,
                                           weight_to_target,
                                           self_loop_weight,
                                           MNF_PROFILE_ARG);
            MNF_PROFILE_END(move_node, move_node);
            ++result.num_moves;

            MNF_PROFILE_BEGIN(neighbor_requeue);
            for (const Edge& e : G.adj[v]) {
                const Vertex u = e.to;
                if (u == v) {
                    continue;
                }
                if (result.partition.community_of[u] != best_community &&
                    !in_queue[u]) {
                    queue.push_back(u);
                    in_queue[u] = true;
                }
            }
            MNF_PROFILE_END(neighbor_requeue, neighbor_requeue);
        }

        if (debug && result.num_visits % debug_interval == 0) {
            std::cerr << "[MoveNodesFast] visits=" << result.num_visits
                      << " moves=" << result.num_moves
                      << " queue=" << queue.size() << "\n";
        }
    }

#ifdef ENABLE_MOVENODESFAST_PROFILE
    result.profile.num_visits = result.num_visits;
#endif
    return result;
}

MoveNodesFastResult MoveNodesFastParallelStage4A(
    const Graph& G,
    const LeidenGraphStats& stats,
    LeidenPartition partition,
    const QualityFunction& quality_function,
    std::uint64_t global_seed,
    std::uint64_t leiden_level,
    const LeidenOptions* options)
{
    // TEMPORARY MOVENODESFAST STAGE4A1 PROFILING
    const Clock::time_point total_begin = Clock::now();
    const int n = num_vertices(G);
    MoveNodesFastResult result;
    result.partition = std::move(partition);

    std::vector<unsigned char> affected_current(
        static_cast<std::size_t>(n), 1);
    std::vector<unsigned char> affected_next(
        static_cast<std::size_t>(n), 0);
    std::vector<MoveProposal> proposals(static_cast<std::size_t>(n));
    std::vector<Vertex> deterministic_order(static_cast<std::size_t>(n));
    std::iota(deterministic_order.begin(), deterministic_order.end(), 0);
    std::sort(deterministic_order.begin(), deterministic_order.end(),
              [global_seed, leiden_level](Vertex lhs, Vertex rhs) {
                  const std::uint64_t lhs_priority =
                      Stage4APriority(global_seed, leiden_level, lhs);
                  const std::uint64_t rhs_priority =
                      Stage4APriority(global_seed, leiden_level, rhs);
                  return lhs_priority < rhs_priority ||
                         (lhs_priority == rhs_priority && lhs < rhs);
              });
    const bool debug = DebugEnabled(options);
    std::size_t round = 0;
    while (true) {
        ++round;
        std::size_t round_visits = 0;
        std::size_t positive_proposals = 0;
        double neighbor_seconds = 0.0;
        double candidate_seconds = 0.0;
        const Clock::time_point proposal_begin = Clock::now();

        // The partition is immutable throughout this region. Each iteration
        // writes only proposals[v], and scratch storage belongs to its thread.
#pragma omp parallel reduction(+:round_visits, positive_proposals, neighbor_seconds, candidate_seconds)
        {
            NeighborCommunityScratch scratch;
            scratch.weights.assign(
                result.partition.community_size.size(), 0.0);
            scratch.marks.assign(
                result.partition.community_size.size(), 0);
            scratch.touched.reserve(32);

#pragma omp for schedule(dynamic, 256)
            for (Vertex v = 0; v < n; ++v) {
                if (affected_current[static_cast<std::size_t>(v)] == 0) {
                    continue;
                }
                ++round_visits;
                MoveProposal proposal;
                proposal.active = true;
                proposal.vertex = v;
                proposal.source = result.partition.community_of[v];
                proposal.target = proposal.source;

                const Clock::time_point neighbor_begin = Clock::now();
                BuildNeighborCommunityWeights(
                    G, result.partition, v, scratch);
                neighbor_seconds +=
                    ElapsedSeconds(neighbor_begin, Clock::now());

                const Clock::time_point candidate_begin = Clock::now();
                const double weight_to_source =
                    LookupNeighborCommunityWeight(scratch, proposal.source);
                Community best = proposal.source;
                double best_delta = 0.0;
                for (Community candidate : scratch.touched) {
                    if (candidate == proposal.source) {
                        continue;
                    }
                    const double delta =
                        quality_function.deltaMoveFromWeights(
                            stats,
                            result.partition,
                            v,
                            candidate,
                            weight_to_source,
                            LookupNeighborCommunityWeight(scratch, candidate));
                    if (delta > best_delta ||
                        (delta == best_delta && candidate < best)) {
                        best_delta = delta;
                        best = candidate;
                    }
                }
                const Community empty =
                    EmptyCommunityForMove(result.partition);
                if (empty != proposal.source) {
                    const double delta =
                        quality_function.deltaMoveFromWeights(
                            stats,
                            result.partition,
                            v,
                            empty,
                            weight_to_source,
                            LookupNeighborCommunityWeight(scratch, empty));
                    if (delta > best_delta ||
                        (delta == best_delta && empty < best)) {
                        best_delta = delta;
                        best = empty;
                    }
                }
                proposal.target = best;
                proposal.delta = best_delta;
                proposal.positive =
                    best_delta > 0.0 && best != proposal.source;
                proposal.target_was_empty =
                    proposal.positive && best == empty;
                if (proposal.positive) {
                    ++positive_proposals;
                }
                proposals[static_cast<std::size_t>(v)] = proposal;
                candidate_seconds +=
                    ElapsedSeconds(candidate_begin, Clock::now());
            }
        }
        const Clock::time_point proposal_end = Clock::now();

        std::size_t committed = 0;
        std::size_t rejected = 0;
        std::size_t next_active = 0;
        Stage4A1ReactivationStats reactivation_stats;
        double revalidation_seconds = 0.0;
        double affected_seconds = 0.0;
        NeighborCommunityScratch commit_scratch;
        commit_scratch.weights.assign(
            result.partition.community_size.size(), 0.0);
        commit_scratch.marks.assign(
            result.partition.community_size.size(), 0);
        commit_scratch.touched.reserve(32);
        const Clock::time_point commit_begin = Clock::now();

        for (Vertex v : deterministic_order) {
            if (affected_current[static_cast<std::size_t>(v)] == 0) {
                continue;
            }
            const MoveProposal& proposal =
                proposals[static_cast<std::size_t>(v)];
            if (!proposal.active || !proposal.positive) {
                continue;
            }

            const Clock::time_point revalidation_begin = Clock::now();
            bool valid = result.partition.community_of[v] == proposal.source;
            if (valid && proposal.target_was_empty) {
                valid = proposal.target >= 0 &&
                        (static_cast<std::size_t>(proposal.target) >=
                             result.partition.community_size.size() ||
                         result.partition.community_is_empty[
                             static_cast<std::size_t>(proposal.target)] != 0);
            }
            double weight_to_source = 0.0;
            double weight_to_target = 0.0;
            double current_delta = 0.0;
            if (valid) {
                BuildNeighborCommunityWeights(
                    G, result.partition, v, commit_scratch);
                weight_to_source = LookupNeighborCommunityWeight(
                    commit_scratch, proposal.source);
                weight_to_target = LookupNeighborCommunityWeight(
                    commit_scratch, proposal.target);
                current_delta = quality_function.deltaMoveFromWeights(
                    stats,
                    result.partition,
                    v,
                    proposal.target,
                    weight_to_source,
                    weight_to_target);
                valid = current_delta > 0.0;
            }
            revalidation_seconds +=
                ElapsedSeconds(revalidation_begin, Clock::now());

            if (!valid) {
                ++rejected;
                // Its snapshot suggestion may have been invalidated by an
                // earlier serial commit. Reconsider it in the next round.
                unsigned char& mark =
                    affected_next[static_cast<std::size_t>(v)];
                if (mark == 0) {
                    mark = 1;
                    ++next_active;
                }
                continue;
            }

            MoveNodeToCommunityFromWeights(
                G,
                stats,
                result.partition,
                v,
                proposal.target,
                weight_to_source,
                weight_to_target,
                SelfLoopWeight(G, v),
                MNF_PROFILE_ARG);
            ++committed;
            ++result.num_moves;

            const Clock::time_point affected_begin = Clock::now();
            // Match serial MoveNodesFast: after the commit, reactivate only
            // non-self neighbors that are not in the committed target. The
            // moved vertex is reconsidered only if a later move affects it.
            const Stage4A1ReactivationStats update =
                UpdateStage4A1AffectedNeighbors(
                    G, result.partition, v, proposal.target, affected_next);
            reactivation_stats.neighbor_scans += update.neighbor_scans;
            reactivation_stats.target_community_exclusions +=
                update.target_community_exclusions;
            reactivation_stats.newly_activated += update.newly_activated;
            reactivation_stats.duplicate_attempts += update.duplicate_attempts;
            next_active += update.newly_activated;
            affected_seconds +=
                ElapsedSeconds(affected_begin, Clock::now());
        }
        const Clock::time_point commit_end = Clock::now();
        result.num_visits += round_visits;

#ifdef ENABLE_MOVENODESFAST_PROFILE
        result.profile.stage4a_active_scan +=
            std::max(0.0,
                     ElapsedSeconds(proposal_begin, proposal_end) -
                         neighbor_seconds - candidate_seconds);
        result.profile.stage4a_neighbor_weights += neighbor_seconds;
        result.profile.stage4a_candidate_evaluation += candidate_seconds;
        result.profile.stage4a_proposal_generation +=
            ElapsedSeconds(proposal_begin, proposal_end);
        result.profile.stage4a_deterministic_commit +=
            ElapsedSeconds(commit_begin, commit_end);
        result.profile.stage4a_commit_revalidation += revalidation_seconds;
        result.profile.stage4a_affected_next_update += affected_seconds;
        ++result.profile.stage4a_rounds;
        result.profile.stage4a_active_vertices += round_visits;
        result.profile.stage4a_positive_proposals += positive_proposals;
        result.profile.stage4a_committed_moves += committed;
        result.profile.stage4a_rejected_proposals += rejected;
        result.profile.stage4a_max_active_vertices = std::max(
            result.profile.stage4a_max_active_vertices, round_visits);
        result.profile.stage4a_reactivation_neighbor_scans +=
            reactivation_stats.neighbor_scans;
        result.profile.stage4a_reactivation_target_exclusions +=
            reactivation_stats.target_community_exclusions;
        result.profile.stage4a_reactivation_new_activations +=
            reactivation_stats.newly_activated;
        result.profile.stage4a_reactivation_duplicate_attempts +=
            reactivation_stats.duplicate_attempts;
#else
        (void)proposal_begin;
        (void)proposal_end;
        (void)commit_begin;
        (void)commit_end;
#endif
        const bool needs_full_verification =
            committed == 0 && round_visits != static_cast<std::size_t>(n);
        if (needs_full_verification) {
            next_active = static_cast<std::size_t>(n);
        }
        if (debug) {
            std::cerr << "[MoveNodesFast Stage4A.1] round=" << round
                      << " active=" << round_visits
                      << " positive=" << positive_proposals
                      << " committed=" << committed
                      << " rejected=" << rejected
                      << " next_active=" << next_active << "\n";
        }
        if (committed == 0 && round_visits == static_cast<std::size_t>(n)) {
            break;
        }

        // A partial affected set cannot prove global local optimality: a
        // target-community neighbor may have become improving without being
        // reactivated by the serial-compatible rule. Before termination, run
        // one deterministic full verification round. If it finds moves, the
        // ordinary double-buffer progression resumes.
        if (needs_full_verification) {
            std::fill(affected_next.begin(), affected_next.end(), 1);
        }

        const Clock::time_point clear_begin = Clock::now();
        affected_current.swap(affected_next);
        std::fill(affected_next.begin(), affected_next.end(), 0);
#ifdef ENABLE_MOVENODESFAST_PROFILE
        result.profile.stage4a_buffer_clear +=
            ElapsedSeconds(clear_begin, Clock::now());
#else
        (void)clear_begin;
#endif
    }

#ifdef ENABLE_MOVENODESFAST_PROFILE
    result.profile.num_visits = result.num_visits;
    result.profile.stage4a_total =
        ElapsedSeconds(total_begin, Clock::now());
#else
    (void)total_begin;
#endif
    return result;
}

#undef MNF_PROFILE_BEGIN
#undef MNF_PROFILE_END
#undef MNF_PROFILE_CANDIDATES
#undef MNF_PROFILE_ADD
#undef MNF_PROFILE_ARG
// END TEMPORARY MOVENODESFAST PERFORMANCE PROFILING

void MergeNodesSubset(const Graph& G,
                      const LeidenGraphStats& stats,
                      LeidenPartition& refined,
                      const std::vector<Vertex>& subset,
                      const QualityFunction& quality_function,
                      double theta,
                      std::mt19937& rng,
                      const std::vector<std::size_t>& subset_mark,
                      std::size_t subset_generation,
                      const LeidenOptions* options)
{
    if (theta <= 0.0) {
        throw std::invalid_argument("theta must be positive");
    }
    if (subset.empty()) {
        return;
    }

    const bool debug = DebugEnabled(options);
    const std::size_t debug_interval = DebugInterval(options);
    std::size_t next_vertex_report = debug_interval;
    std::size_t processed_vertices = 0;

    // if (debug) {
    //     std::cerr << "[MergeNodesSubset] start subset_size=" << subset.size()
    //               << "\n";
    // }

    const double subset_mass = SubsetMass(stats, quality_function, subset);
    std::vector<Vertex> candidates;
    candidates.reserve(subset.size());
    for (Vertex v : subset) {
        if (IsNodeWellConnectedToSubset(G,
                                        stats,
                                        quality_function,
                                        v,
                                        subset_mark,
                                        subset_generation,
                                        subset_mass)) {
            candidates.push_back(v);
        }
    }
    std::shuffle(candidates.begin(), candidates.end(), rng);

    RefinementCommunityStats community_stats =
        BuildRefinementCommunityStats(G,
                                      stats,
                                      refined,
                                      quality_function,
                                      subset,
                                      subset_mark,
                                      subset_generation);
    for (Vertex v : candidates) {
        ++processed_vertices;
        const Community source = refined.community_of[v];
        const auto source_entry = community_stats.entries.find(source);
        if (source < 0 ||
            source_entry == community_stats.entries.end() ||
            source_entry->second.member_count != 1) {
            if (debug && processed_vertices >= next_vertex_report) {
                std::cerr << "[MergeNodesSubset] processed_vertices="
                          << processed_vertices
                          << " candidate_vertices=" << candidates.size()
                          << " current_subset_size=" << subset.size()
                          << " active_refined_communities="
                          << community_stats.active_communities.size() << "\n";
                next_vertex_report += debug_interval;
            }
            continue;
        }

        const auto neighbor_weights =
            BuildNeighborCommunityWeights(G, refined, v);
        const double weight_to_source = LookupWeight(neighbor_weights, source);

        struct WeightedCandidate {
            Community community;
            double delta;
        };
        std::vector<WeightedCandidate> nondecreasing_candidates;
        nondecreasing_candidates.reserve(community_stats.active_communities.size());

        // Community-level full subset rescans are avoided here. The subset
        // stats are built once in O(|S| + |E_S|); after each merge they are
        // updated from G.adj[v]. Each well-connectedness query below is O(1).
        for (Community community : community_stats.active_communities) {
            if (!IsCommunityWellConnectedFromStats(stats,
                                                  quality_function,
                                                  community,
                                                  subset_mass,
                                                  community_stats)) {
                continue;
            }

            const double delta =
                (community == source)
                    ? 0.0
                    : quality_function.deltaMoveFromWeights(
                          stats,
                          refined,
                          v,
                          community,
                          weight_to_source,
                          LookupWeight(neighbor_weights, community));

            // Algorithm A.2 admits refinement merges with Delta H >= 0,
            // unlike MoveNodesFast, which uses strictly positive moves.
            if (delta >= 0.0) {
                nondecreasing_candidates.push_back({community, delta});
            }
        }

        if (nondecreasing_candidates.empty()) {
            continue;
        }

        double max_delta = nondecreasing_candidates.front().delta;
        for (const WeightedCandidate& candidate : nondecreasing_candidates) {
            max_delta = std::max(max_delta, candidate.delta);
        }

        std::vector<double> weights;
        weights.reserve(nondecreasing_candidates.size());
        for (const WeightedCandidate& candidate : nondecreasing_candidates) {
            weights.push_back(std::exp((candidate.delta - max_delta) / theta));
        }

        std::discrete_distribution<std::size_t> distribution(weights.begin(),
                                                             weights.end());
        const Community target =
            nondecreasing_candidates[distribution(rng)].community;

        if (target != source) {
            UpdateRefinementCommunityStatsForMove(G,
                                                  stats,
                                                  refined,
                                                  quality_function,
                                                  subset_mark,
                                                  subset_generation,
                                                  v,
                                                  target,
                                                  community_stats);
            MoveNodeToCommunity(G, stats, refined, v, target);
        }

        if (debug && processed_vertices >= next_vertex_report) {
            std::cerr << "[MergeNodesSubset] processed_vertices="
                      << processed_vertices
                      << " candidate_vertices=" << candidates.size()
                      << " current_subset_size=" << subset.size()
                      << " active_refined_communities="
                      << community_stats.active_communities.size() << "\n";
            next_vertex_report += debug_interval;
        }
    }

    // if (debug) {
    //     std::cerr << "[MergeNodesSubset] done processed_vertices="
    //               << processed_vertices
    //               << " candidate_vertices=" << candidates.size()
    //               << " current_subset_size=" << subset.size()
    //               << " active_refined_communities="
    //               << community_stats.active_communities.size() << "\n";
    // }
}

LocalRefinementResult MergeNodesSubsetLocal(
    const Graph& G,
    const LeidenGraphStats& stats,
    const LeidenPartition& parent_partition,
    const std::vector<double>& singleton_internal_edge_weight,
    LeidenPartition& local_refined,
    const std::vector<Vertex>& subset,
    Community parent,
    const QualityFunction& quality_function,
    double theta,
    std::mt19937& local_rng
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
    , RefineSubsetProfile* detail
#endif
    )
{
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
    const Clock::time_point detail_total_begin = Clock::now();
    detail->parent = parent;
    detail->vertices = subset.size();
    for (Vertex v : subset) {
        detail->adjacency_volume += G.adj[v].size();
    }
#endif
    if (theta <= 0.0) {
        throw std::invalid_argument("theta must be positive");
    }

    LocalRefinementResult result;
    result.vertices = subset;
    result.local_assignment.assign(subset.size(), 0);
    if (subset.empty()) {
        return result;
    }

#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
    Clock::time_point detail_begin = Clock::now();
#endif
    ResetLocalSingletonsForSubset(G,
                                  stats,
                                  singleton_internal_edge_weight,
                                  local_refined,
                                  subset);
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
    detail->reset_time = ElapsedSeconds(detail_begin, Clock::now());
    detail_begin = Clock::now();
#endif

    const double subset_mass = SubsetMass(stats, quality_function, subset);
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
    detail->subset_mass_time = ElapsedSeconds(detail_begin, Clock::now());
    detail_begin = Clock::now();
#endif
    std::vector<Vertex> candidates;
    candidates.reserve(subset.size());
    for (Vertex v : subset) {
        if (IsNodeWellConnectedToParentSubset(G,
                                              stats,
                                              parent_partition,
                                              quality_function,
                                              v,
                                              parent,
                                              subset_mass)) {
            candidates.push_back(v);
        }
    }
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
    detail->candidates = candidates.size();
    detail->candidate_scan_time = ElapsedSeconds(detail_begin, Clock::now());
    detail_begin = Clock::now();
#endif
    std::shuffle(candidates.begin(), candidates.end(), local_rng);
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
    detail->shuffle_time = ElapsedSeconds(detail_begin, Clock::now());
    detail_begin = Clock::now();
#endif

    NeighborCommunityScratch neighbor_scratch;
    neighbor_scratch.weights.assign(local_refined.community_size.size(), 0.0);
    neighbor_scratch.marks.assign(local_refined.community_size.size(), 0);
    neighbor_scratch.touched.reserve(32);
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
    detail->scratch_init_time = ElapsedSeconds(detail_begin, Clock::now());
    detail_begin = Clock::now();
#endif

    RefinementCommunityStats community_stats =
        BuildRefinementCommunityStatsForParent(G,
                                               stats,
                                               parent_partition,
                                               local_refined,
                                               quality_function,
                                               subset,
                                               parent);
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
    detail->stats_build_time = ElapsedSeconds(detail_begin, Clock::now());
    detail->stats_build_vertices = subset.size();
    detail->stats_build_adjacency_scans = detail->adjacency_volume;
    detail_begin = Clock::now();
#endif

    for (Vertex v : candidates) {
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
        ++detail->processed_candidates;
#endif
        const Community source = local_refined.community_of[v];
        const auto source_entry = community_stats.entries.find(source);
        if (source < 0 ||
            source_entry == community_stats.entries.end() ||
            source_entry->second.member_count != 1) {
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
            ++detail->singleton_source_skips;
#endif
            continue;
        }

        BuildNeighborCommunityWeights(G, local_refined, v, neighbor_scratch);
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
        ++detail->neighbor_builds;
        detail->neighbor_adjacency_scans += G.adj[v].size();
        detail->touched_communities += neighbor_scratch.touched.size();
        detail->active_community_iterations +=
            community_stats.active_communities.size();
#endif
        const double weight_to_source =
            LookupNeighborCommunityWeight(neighbor_scratch, source);

        struct WeightedCandidate {
            Community community;
            double delta;
        };
        std::vector<WeightedCandidate> nondecreasing_candidates;
        nondecreasing_candidates.reserve(
            community_stats.active_communities.size());

        for (Community community : community_stats.active_communities) {
            if (!IsCommunityWellConnectedFromStats(stats,
                                                  quality_function,
                                                  community,
                                                  subset_mass,
                                                  community_stats)) {
                continue;
            }

            const double delta =
                (community == source)
                    ? 0.0
                    : quality_function.deltaMoveFromWeights(
                          stats,
                          local_refined,
                          v,
                          community,
                          weight_to_source,
                          LookupNeighborCommunityWeight(neighbor_scratch,
                                                        community));
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
            if (community != source) {
                ++detail->delta_evaluations;
            }
#endif

            if (delta >= 0.0) {
                nondecreasing_candidates.push_back({community, delta});
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
                ++detail->nondecreasing_candidates;
#endif
            }
        }

        if (nondecreasing_candidates.empty()) {
            continue;
        }

        double max_delta = nondecreasing_candidates.front().delta;
        for (const WeightedCandidate& candidate : nondecreasing_candidates) {
            max_delta = std::max(max_delta, candidate.delta);
        }

        std::vector<double> weights;
        weights.reserve(nondecreasing_candidates.size());
        for (const WeightedCandidate& candidate : nondecreasing_candidates) {
            weights.push_back(std::exp((candidate.delta - max_delta) / theta));
        }

        std::discrete_distribution<std::size_t> distribution(weights.begin(),
                                                             weights.end());
        const Community target =
            nondecreasing_candidates[distribution(local_rng)].community;
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
        ++detail->stochastic_selections;
#endif

        if (target != source) {
            UpdateRefinementCommunityStatsForParentMove(G,
                                                        stats,
                                                        parent_partition,
                                                        local_refined,
                                                        quality_function,
                                                        v,
                                                        parent,
                                                        target,
                                                        community_stats);
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
            ++detail->stats_update_calls;
            detail->stats_update_adjacency_scans += G.adj[v].size();
            ++detail->refinement_moves;
#endif
            MoveNodeToCommunity(G, stats, local_refined, v, target);
        }
    }
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
    detail->candidate_processing_time =
        ElapsedSeconds(detail_begin, Clock::now());
    detail_begin = Clock::now();
#endif

    std::vector<Community> local_communities;
    local_communities.reserve(subset.size());
    for (Vertex v : subset) {
        local_communities.push_back(local_refined.community_of[v]);
    }
    std::vector<Community> compact_map = BuildCompactCommunityMap(local_communities);
    result.local_community_count =
        static_cast<Community>(compact_map.empty() ? 0 :
            *std::max_element(compact_map.begin(), compact_map.end()) + 1);
    for (std::size_t i = 0; i < subset.size(); ++i) {
        result.local_assignment[i] =
            compact_map[local_refined.community_of[subset[i]]];
    }
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
    detail->local_community_count = result.local_community_count;
    detail->result_compaction_time =
        ElapsedSeconds(detail_begin, Clock::now());
    detail->total_time =
        ElapsedSeconds(detail_total_begin, Clock::now());
#endif
    return result;
}

LeidenPartition RefinePartition(const Graph& G,
                                const LeidenGraphStats& stats,
                                const LeidenPartition& partition,
                                const QualityFunction& quality_function,
                                double theta,
                                std::uint64_t global_seed,
                                std::uint64_t level,
                                const LeidenOptions* options)
{
    if (theta <= 0.0) {
        throw std::invalid_argument("theta must be positive");
    }

#ifdef ENABLE_REFINEPARTITION_PROFILE
    // TEMPORARY REFINEPARTITION PERFORMANCE PROFILING
    const Clock::time_point refine_total_begin = Clock::now();
    double member_construction_time = 0.0;
    double subset_processing_time = 0.0;
    double global_id_prefix_time = 0.0;
    double global_assignment_time = 0.0;
    double make_partition_time = 0.0;
    const Clock::time_point member_construction_begin = Clock::now();
#endif
    const std::vector<std::vector<Vertex>> members =
        BuildPartitionMembers(partition);
#ifdef ENABLE_REFINEPARTITION_PROFILE
    member_construction_time =
        ElapsedSeconds(member_construction_begin, Clock::now());
#endif
    const bool debug = DebugEnabled(options);
    const std::size_t total_parent_communities =
        static_cast<std::size_t>(
            std::count_if(members.begin(),
                          members.end(),
                          [](const std::vector<Vertex>& subset) {
                              return !subset.empty();
                          }));

    std::vector<LocalRefinementResult> local_results(members.size());
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
    // TEMPORARY REFINEPARTITION DETAILED PROFILING
    std::vector<RefineSubsetProfile> subset_profiles(members.size());
#ifdef _OPENMP
    const int refine_profile_threads = omp_get_max_threads();
#else
    const int refine_profile_threads = 1;
#endif
    std::vector<RefineThreadProfile> thread_profiles(
        static_cast<std::size_t>(refine_profile_threads));
#endif
    std::vector<double> singleton_internal_edge_weight(num_vertices(G), 0.0);
    for (Vertex v = 0; v < num_vertices(G); ++v) {
        singleton_internal_edge_weight[v] = SelfLoopWeight(G, v);
    }

    // Parent-community parallel refinement: each worker mutates only its
    // thread-local partition scratch and stores one result per parent.
#ifdef ENABLE_REFINEPARTITION_PROFILE
    // TEMPORARY REFINEPARTITION PERFORMANCE PROFILING
    const Clock::time_point subset_processing_begin = Clock::now();
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
#ifdef _OPENMP
        const int refine_thread_id = omp_get_thread_num();
#else
        const int refine_thread_id = 0;
#endif
        RefineThreadProfile& thread_profile =
            thread_profiles[static_cast<std::size_t>(refine_thread_id)];
        const Clock::time_point thread_init_begin = Clock::now();
#endif
        LeidenPartition local_refined =
            MakeThreadLocalSingletonPartition(G,
                                              stats,
                                              singleton_internal_edge_weight);
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
        thread_profile.initialization_time =
            ElapsedSeconds(thread_init_begin, Clock::now());
#endif

#ifdef _OPENMP
#ifdef REFINEMENT_DYNAMIC_SCHEDULE
#pragma omp for schedule(dynamic, 1)
#else
#pragma omp for schedule(static)
#endif
#endif
        for (Community parent = 0;
             parent < static_cast<Community>(members.size());
             ++parent) {
            const std::vector<Vertex>& subset = members[parent];
            if (subset.empty()) {
                continue;
            }

            std::mt19937 local_rng =
                MakeSubsetRng(global_seed,
                              level,
                              static_cast<std::uint64_t>(parent));
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
            const Clock::time_point subset_work_begin = Clock::now();
#endif
            local_results[parent] =
                MergeNodesSubsetLocal(G,
                                      stats,
                                      partition,
                                      singleton_internal_edge_weight,
                                      local_refined,
                                      subset,
                                      parent,
                                      quality_function,
                                      theta,
                                      local_rng
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
                                      , &subset_profiles[
                                            static_cast<std::size_t>(parent)]
#endif
                                      );
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
            const double subset_work =
                ElapsedSeconds(subset_work_begin, Clock::now());
            ++thread_profile.subsets;
            thread_profile.vertices += subset.size();
            thread_profile.adjacency_volume +=
                subset_profiles[static_cast<std::size_t>(parent)]
                    .adjacency_volume;
            thread_profile.subset_work_time += subset_work;
#endif
        }
    }
#ifdef ENABLE_REFINEPARTITION_PROFILE
    subset_processing_time =
        ElapsedSeconds(subset_processing_begin, Clock::now());
#endif

#ifdef ENABLE_REFINEPARTITION_PROFILE
    // TEMPORARY REFINEPARTITION PERFORMANCE PROFILING
    const Clock::time_point global_id_prefix_begin = Clock::now();
#endif
    std::vector<Community> community_offsets(members.size() + 1, 0);
    for (std::size_t parent = 0; parent < members.size(); ++parent) {
        community_offsets[parent + 1] =
            community_offsets[parent] +
            local_results[parent].local_community_count;
    }
#ifdef ENABLE_REFINEPARTITION_PROFILE
    global_id_prefix_time =
        ElapsedSeconds(global_id_prefix_begin, Clock::now());
#endif

#ifdef ENABLE_REFINEPARTITION_PROFILE
    // TEMPORARY REFINEPARTITION PERFORMANCE PROFILING
    const Clock::time_point global_assignment_begin = Clock::now();
#endif
    std::vector<Community> global_assignment(num_vertices(G), -1);
    std::size_t processed_vertices = 0;
    for (std::size_t parent = 0; parent < members.size(); ++parent) {
        const LocalRefinementResult& local = local_results[parent];
        const Community offset = community_offsets[parent];
        for (std::size_t i = 0; i < local.vertices.size(); ++i) {
            global_assignment[local.vertices[i]] =
                offset + local.local_assignment[i];
        }
        processed_vertices += local.vertices.size();
    }
#ifdef ENABLE_REFINEPARTITION_PROFILE
    global_assignment_time =
        ElapsedSeconds(global_assignment_begin, Clock::now());
#endif

    if (debug) {
        std::cerr << "[RefinePartition] done processed_parent_communities="
                  << total_parent_communities
                  << " total_parent_communities=" << total_parent_communities
                  << " processed_vertices=" << processed_vertices
                  << "\n";
    }

#ifdef ENABLE_REFINEPARTITION_PROFILE
    // TEMPORARY REFINEPARTITION PERFORMANCE PROFILING
    const Clock::time_point make_partition_begin = Clock::now();
#endif
    LeidenPartition refined = MakePartition(G, stats, global_assignment);
#ifdef ENABLE_REFINEPARTITION_PROFILE
    make_partition_time = ElapsedSeconds(make_partition_begin, Clock::now());
    const double refine_total_time =
        ElapsedSeconds(refine_total_begin, Clock::now());
    if (debug) {
        std::cerr << "[RefinePartition profile] member construction: "
                  << member_construction_time << " s\n"
                  << "[RefinePartition profile] subset processing: "
                  << subset_processing_time << " s\n"
                  << "[RefinePartition profile] local result storage: "
                  << subset_processing_time << " s (included in subset processing)\n"
                  << "[RefinePartition profile] global ID prefix: "
                  << global_id_prefix_time << " s\n"
                  << "[RefinePartition profile] global assignment reconstruction: "
                  << global_assignment_time << " s\n"
                  << "[RefinePartition profile] MakePartition rebuild: "
                  << make_partition_time << " s\n"
                  << "[RefinePartition profile] total: "
                  << refine_total_time << " s\n";
    }
    // END TEMPORARY REFINEPARTITION PERFORMANCE PROFILING
#endif
#ifdef ENABLE_REFINEPARTITION_DETAILED_PROFILE
    // TEMPORARY REFINEPARTITION DETAILED PROFILING
    RefineSubsetProfile aggregate_detail;
    std::size_t subset_size_1 = 0;
    std::size_t subset_size_2_10 = 0;
    std::size_t subset_size_11_100 = 0;
    std::size_t subset_size_101_1000 = 0;
    std::size_t subset_size_gt_1000 = 0;
    std::size_t max_subset_size = 0;
    std::size_t max_adjacency_volume = 0;
    std::size_t min_subset_size = std::numeric_limits<std::size_t>::max();
    Community max_local_community_count = 0;
    std::size_t total_local_community_count = 0;
    std::vector<std::size_t> heavy_indices;
    heavy_indices.reserve(total_parent_communities);
    for (std::size_t parent = 0; parent < subset_profiles.size(); ++parent) {
        const RefineSubsetProfile& p = subset_profiles[parent];
        if (p.vertices == 0) {
            continue;
        }
        heavy_indices.push_back(parent);
        max_subset_size = std::max(max_subset_size, p.vertices);
        min_subset_size = std::min(min_subset_size, p.vertices);
        max_adjacency_volume =
            std::max(max_adjacency_volume, p.adjacency_volume);
        if (p.vertices == 1) ++subset_size_1;
        else if (p.vertices <= 10) ++subset_size_2_10;
        else if (p.vertices <= 100) ++subset_size_11_100;
        else if (p.vertices <= 1000) ++subset_size_101_1000;
        else ++subset_size_gt_1000;
#define REFINE_SUM(field) aggregate_detail.field += p.field
        REFINE_SUM(vertices);
        REFINE_SUM(adjacency_volume);
        REFINE_SUM(candidates);
        REFINE_SUM(processed_candidates);
        REFINE_SUM(singleton_source_skips);
        REFINE_SUM(neighbor_builds);
        REFINE_SUM(neighbor_adjacency_scans);
        REFINE_SUM(touched_communities);
        REFINE_SUM(active_community_iterations);
        REFINE_SUM(delta_evaluations);
        REFINE_SUM(nondecreasing_candidates);
        REFINE_SUM(stochastic_selections);
        REFINE_SUM(refinement_moves);
        REFINE_SUM(stats_build_vertices);
        REFINE_SUM(stats_build_adjacency_scans);
        REFINE_SUM(stats_update_calls);
        REFINE_SUM(stats_update_adjacency_scans);
        REFINE_SUM(reset_time);
        REFINE_SUM(subset_mass_time);
        REFINE_SUM(candidate_scan_time);
        REFINE_SUM(shuffle_time);
        REFINE_SUM(scratch_init_time);
        REFINE_SUM(stats_build_time);
        REFINE_SUM(candidate_processing_time);
        REFINE_SUM(result_compaction_time);
        REFINE_SUM(total_time);
#undef REFINE_SUM
        total_local_community_count +=
            static_cast<std::size_t>(p.local_community_count);
        max_local_community_count =
            std::max(max_local_community_count, p.local_community_count);
    }
    std::sort(heavy_indices.begin(), heavy_indices.end(),
              [&subset_profiles](std::size_t lhs, std::size_t rhs) {
                  return subset_profiles[lhs].total_time >
                         subset_profiles[rhs].total_time;
              });

    double init_total = 0.0;
    double init_max = 0.0;
    double thread_work_total = 0.0;
    double thread_work_max = 0.0;
    double thread_work_min = std::numeric_limits<double>::max();
    for (const RefineThreadProfile& p : thread_profiles) {
        init_total += p.initialization_time;
        init_max = std::max(init_max, p.initialization_time);
        thread_work_total += p.subset_work_time;
        thread_work_max = std::max(thread_work_max, p.subset_work_time);
        thread_work_min = std::min(thread_work_min, p.subset_work_time);
    }
    const double thread_count = static_cast<double>(thread_profiles.size());
    const double thread_work_average = thread_work_total / thread_count;
    const std::size_t community_slots = partition.community_size.size();
    const std::size_t empty_slots = static_cast<std::size_t>(std::count(
        partition.community_is_empty.begin(),
        partition.community_is_empty.end(),
        static_cast<unsigned char>(1)));
    const auto minmax_id = std::minmax_element(partition.community_of.begin(),
                                               partition.community_of.end());
    const std::size_t scratch_bytes_per_thread =
        static_cast<std::size_t>(num_vertices(G)) *
        (sizeof(Community) + 3 * sizeof(double) + sizeof(unsigned char));

    std::cout
        << "[RefineDetailed] level=" << level
        << " n=" << num_vertices(G)
        << " m=" << num_edges(G)
        << " community_slots=" << community_slots
        << " active_parents=" << total_parent_communities
        << " empty_slots=" << empty_slots
        << " min_id=" << (partition.community_of.empty() ? -1 : *minmax_id.first)
        << " max_id=" << (partition.community_of.empty() ? -1 : *minmax_id.second)
        << " member_vectors=" << members.size() << "\n"
#ifdef ENABLE_REFINEPARTITION_PROFILE
        << "[RefineDetailed] wall_times member_construction="
        << member_construction_time
        << " subset_processing=" << subset_processing_time
        << " global_id_prefix=" << global_id_prefix_time
        << " global_assignment=" << global_assignment_time
        << " make_partition=" << make_partition_time
        << " total=" << refine_total_time << "\n"
#endif
        << "[RefineDetailed] subsets=" << total_parent_communities
        << " vertices=" << aggregate_detail.vertices
        << " min_size=" << (total_parent_communities == 0 ? 0 : min_subset_size)
        << " max_size=" << max_subset_size
        << " avg_size=" << (total_parent_communities == 0 ? 0.0 :
             static_cast<double>(aggregate_detail.vertices) /
             static_cast<double>(total_parent_communities))
        << " adjacency_volume=" << aggregate_detail.adjacency_volume
        << " max_adjacency_volume=" << max_adjacency_volume
        << " avg_adjacency_volume=" << (total_parent_communities == 0 ? 0.0 :
             static_cast<double>(aggregate_detail.adjacency_volume) /
             static_cast<double>(total_parent_communities)) << "\n"
        << "[RefineDetailed] size_bins=" << subset_size_1 << ","
        << subset_size_2_10 << "," << subset_size_11_100 << ","
        << subset_size_101_1000 << "," << subset_size_gt_1000 << "\n"
        << "[RefineDetailed] candidates=" << aggregate_detail.candidates
        << " processed=" << aggregate_detail.processed_candidates
        << " singleton_skips=" << aggregate_detail.singleton_source_skips
        << " neighbor_builds=" << aggregate_detail.neighbor_builds
        << " neighbor_adjacency_scans="
        << aggregate_detail.neighbor_adjacency_scans
        << " touched=" << aggregate_detail.touched_communities
        << " active_iterations=" << aggregate_detail.active_community_iterations
        << " delta_evaluations=" << aggregate_detail.delta_evaluations
        << " nondecreasing=" << aggregate_detail.nondecreasing_candidates
        << " selections=" << aggregate_detail.stochastic_selections
        << " moves=" << aggregate_detail.refinement_moves << "\n"
        << "[RefineDetailed] stats_build_calls=" << total_parent_communities
        << " stats_build_vertices=" << aggregate_detail.stats_build_vertices
        << " stats_build_adjacency_scans="
        << aggregate_detail.stats_build_adjacency_scans
        << " stats_update_calls=" << aggregate_detail.stats_update_calls
        << " stats_update_adjacency_scans="
        << aggregate_detail.stats_update_adjacency_scans << "\n"
        << "[RefineDetailed] cpu_times reset=" << aggregate_detail.reset_time
        << " mass=" << aggregate_detail.subset_mass_time
        << " candidate_scan=" << aggregate_detail.candidate_scan_time
        << " shuffle=" << aggregate_detail.shuffle_time
        << " scratch_init=" << aggregate_detail.scratch_init_time
        << " stats_build=" << aggregate_detail.stats_build_time
        << " candidate_processing=" << aggregate_detail.candidate_processing_time
        << " compaction=" << aggregate_detail.result_compaction_time
        << " subset_total=" << aggregate_detail.total_time << "\n"
        << "[RefineDetailed] local_communities_total="
        << total_local_community_count
        << " average=" << (total_parent_communities == 0 ? 0.0 :
             static_cast<double>(total_local_community_count) /
             static_cast<double>(total_parent_communities))
        << " max=" << max_local_community_count
        << " local_result_vertices=" << processed_vertices
        << " local_assignment_entries=" << processed_vertices << "\n"
        << "[RefineDetailed] threads=" << thread_profiles.size()
        << " init_total_cpu=" << init_total
        << " init_max_wall_proxy=" << init_max
        << " init_avg=" << init_total / thread_count
        << " scratch_bytes_per_thread=" << scratch_bytes_per_thread
        << " scratch_bytes_total="
        << scratch_bytes_per_thread * thread_profiles.size() << "\n"
        << "[RefineDetailed] thread_work_min=" << thread_work_min
        << " max=" << thread_work_max
        << " avg=" << thread_work_average
        << " max_avg=" << (thread_work_average == 0.0 ? 0.0 :
             thread_work_max / thread_work_average) << "\n";
    for (std::size_t rank = 0;
         rank < std::min<std::size_t>(5, heavy_indices.size()); ++rank) {
        const RefineSubsetProfile& p = subset_profiles[heavy_indices[rank]];
        std::cout << "[RefineDetailed top] rank=" << rank + 1
                  << " parent=" << p.parent
                  << " size=" << p.vertices
                  << " adjacency_volume=" << p.adjacency_volume
                  << " time=" << p.total_time
                  << " local_communities=" << p.local_community_count << "\n";
    }
#endif
    return refined;
}

LeidenPartition RefinePartition(const Graph& G,
                                const LeidenGraphStats& stats,
                                const LeidenPartition& partition,
                                const QualityFunction& quality_function,
                                double theta,
                                std::mt19937& rng,
                                const LeidenOptions* options)
{
    return RefinePartition(G,
                           stats,
                           partition,
                           quality_function,
                           theta,
                           static_cast<std::uint64_t>(rng()),
                           0,
                           options);
}

std::uint64_t MakeSubsetSeed(std::uint64_t global_seed,
                             std::uint64_t level,
                             std::uint64_t parent)
{
    auto mix = [](std::uint64_t x) {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30U)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27U)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31U);
    };

    std::uint64_t seed = mix(global_seed);
    seed ^= mix(level + 0x632be59bd9b4e019ULL);
    seed ^= mix(parent + 0x85157af5ULL);
    return mix(seed);
}

AggregateGraphResult AggregateGraph(const Graph& G,
                                    const LeidenGraphStats& stats,
                                    const LeidenPartition& refined)
{
    ValidateAggregateInput(G, stats, refined);

    AggregateGraphResult result;
    const std::vector<Community> community_to_coarse =
        BuildCommunityToCoarse(refined);
    const int n_coarse =
        community_to_coarse.empty()
            ? 0
            : static_cast<int>(
                  *std::max_element(community_to_coarse.begin(),
                                    community_to_coarse.end()) + 1);

    result.graph = MakeGraph(n_coarse);
    result.coarse_of.assign(num_vertices(G), -1);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (Vertex v = 0; v < num_vertices(G); ++v) {
        result.coarse_of[v] = community_to_coarse[refined.community_of[v]];
    }

    std::vector<std::size_t> coarse_counts(n_coarse, 0);
    for (Vertex v = 0; v < num_vertices(G); ++v) {
        ++coarse_counts[result.coarse_of[v]];
    }

    std::vector<std::size_t> coarse_offsets(n_coarse + 1, 0);
    for (Community c = 0; c < n_coarse; ++c) {
        coarse_offsets[static_cast<std::size_t>(c) + 1] =
            coarse_offsets[c] + coarse_counts[c];
    }

    std::vector<Vertex> coarse_vertices(num_vertices(G));
    std::vector<std::size_t> next_offset = coarse_offsets;
    for (Vertex v = 0; v < num_vertices(G); ++v) {
        const Community coarse = result.coarse_of[v];
        coarse_vertices[next_offset[coarse]++] = v;
    }

    std::vector<double> coarse_node_size(n_coarse, 0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (Community c = 0; c < n_coarse; ++c) {
        double size = 0.0;
        for (std::size_t i = coarse_offsets[c]; i < coarse_offsets[c + 1];
             ++i) {
            size += stats.node_size[coarse_vertices[i]];
        }
        coarse_node_size[c] = size;
    }

    std::vector<AggregateCommunityEdges> community_edges(n_coarse);

    // Community-parallel aggregation: each thread reuses private
    // neighbor-community scratch and each community owns its output slot.
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        AggregateScratch scratch(n_coarse);

#ifdef _OPENMP
#ifdef AGGREGATEGRAPH_STATIC_SCHEDULE
#pragma omp for schedule(static)
#else
#pragma omp for schedule(dynamic, 1)
#endif
#endif
        for (Community c = 0; c < n_coarse; ++c) {
            scratch.reset();

            bool has_internal_edge = false;
            double internal_self_loop = 0.0;
            double internal_nonself_sum = 0.0;

            for (std::size_t i = coarse_offsets[c]; i < coarse_offsets[c + 1];
                 ++i) {
                const Vertex u = coarse_vertices[i];
                for (const Edge& e : G.adj[u]) {
                    const Community d = result.coarse_of[e.to];
                    if (d == c) {
                        has_internal_edge = true;
                        if (e.to == u) {
                            internal_self_loop += e.weight;
                        } else {
                            internal_nonself_sum += e.weight;
                        }
                    } else if (c < d) {
                        scratch.add(d, e.weight);
                    }
                }
            }

            AggregateCommunityEdges& output = community_edges[c];
            output.has_self_loop = has_internal_edge;
            output.self_loop_weight =
                internal_self_loop + 0.5 * internal_nonself_sum;

            std::sort(scratch.touched.begin(), scratch.touched.end());
            output.upper_edges.reserve(scratch.touched.size());
            for (Community d : scratch.touched) {
                output.upper_edges.push_back({d, scratch.weights[d]});
            }
        }
    }

    std::vector<std::size_t> coarse_degrees(n_coarse, 0);
    for (Community c = 0; c < n_coarse; ++c) {
        if (community_edges[c].has_self_loop) {
            ++coarse_degrees[c];
        }
        coarse_degrees[c] += community_edges[c].upper_edges.size();
        for (const Edge& e : community_edges[c].upper_edges) {
            ++coarse_degrees[e.to];
        }
    }

    for (Community c = 0; c < n_coarse; ++c) {
        result.graph.adj[c].reserve(coarse_degrees[c]);
    }
    for (Community c = 0; c < n_coarse; ++c) {
        const AggregateCommunityEdges& edges = community_edges[c];
        if (edges.has_self_loop) {
            result.graph.adj[c].push_back({c, edges.self_loop_weight});
        }
        for (const Edge& e : edges.upper_edges) {
            result.graph.adj[c].push_back(e);
            result.graph.adj[e.to].push_back({c, e.weight});
        }
    }

    result.stats = BuildLeidenGraphStats(result.graph, coarse_node_size);
    return result;
}

LeidenPartition BuildCoarsePartition(const AggregateGraphResult& aggregate,
                                     const LeidenPartition& partition,
                                     const LeidenPartition& refined)
{
    const int n_coarse = num_vertices(aggregate.graph);
    const std::size_t n = aggregate.coarse_of.size();
    if (partition.community_of.size() != n ||
        refined.community_of.size() != n) {
        throw std::invalid_argument("coarse partition input size mismatch");
    }

    std::vector<Community> coarse_original_community(n_coarse, -1);
    std::vector<Vertex> refined_coarse(refined.community_size.size(), -1);

    // Aggregation is based on Prefined, but the partition on the aggregate
    // graph represents the original non-refined P. Each coarse vertex
    // corresponds to one refined community, and all original vertices in that
    // coarse vertex must therefore have the same P community.
    for (Vertex v = 0; v < static_cast<Vertex>(n); ++v) {
        const Vertex coarse = aggregate.coarse_of[v];
        const Community p_community = partition.community_of[v];
        const Community refined_community = refined.community_of[v];

        if (coarse < 0 || coarse >= n_coarse) {
            throw std::invalid_argument("coarse_of contains out-of-range vertex id");
        }
        if (p_community < 0 || refined_community < 0) {
            throw std::invalid_argument("partition contains negative community id");
        }

        if (static_cast<std::size_t>(refined_community) >= refined_coarse.size()) {
            refined_coarse.resize(static_cast<std::size_t>(refined_community) + 1,
                                  -1);
        }
        if (refined_coarse[refined_community] == -1) {
            refined_coarse[refined_community] = coarse;
        } else if (refined_coarse[refined_community] != coarse) {
            throw std::invalid_argument("coarse_of is inconsistent with refined partition");
        }

        if (coarse_original_community[coarse] == -1) {
            coarse_original_community[coarse] = p_community;
        } else if (coarse_original_community[coarse] != p_community) {
            throw std::invalid_argument(
                "refined partition is not a refinement of partition");
        }
    }

    for (Vertex coarse = 0; coarse < n_coarse; ++coarse) {
        if (coarse_original_community[coarse] < 0) {
            throw std::invalid_argument("coarse vertex has no original vertex");
        }
    }

    const std::vector<Community> compact_map =
        BuildCompactCommunityMap(coarse_original_community);
    std::vector<Community> coarse_assignment(n_coarse, -1);
    for (Vertex coarse = 0; coarse < n_coarse; ++coarse) {
        coarse_assignment[coarse] =
            compact_map[coarse_original_community[coarse]];
    }

    return MakePartition(aggregate.graph, aggregate.stats, coarse_assignment);
}

LeidenResult Leiden(const Graph& G,
                    const LeidenGraphStats& stats,
                    const QualityFunction& quality_function,
                    const LeidenOptions& options)
{
    ValidateLeidenInput(G, stats, options);

    // -------------------------------------------------------------------------
    // TEMPORARY PERFORMANCE INSTRUMENTATION
    // Phase durations below reuse the existing debug timestamps and accumulate
    // only the four profiled calls. Keep these additions local to Leiden().
    // -------------------------------------------------------------------------
    const Clock::time_point leiden_begin = Clock::now();

    LeidenResult result;
    const int n_original = num_vertices(G);
    if (n_original == 0) {
        result.partition = MakePartition(G, stats, {});
        if (options.debug) {
            std::cerr << "[Leiden] converged because num_moves == 0\n";
        }
        result.timing.total =
            ElapsedSeconds(leiden_begin, Clock::now());
        return result;
    }

    if (options.debug)
        std::cerr << "[Leiden] graph copy start\n";

    const auto t0 = Clock::now();
    Graph current_graph = G;
    const auto t1 = Clock::now();

    if (options.debug)
        std::cerr << "[Leiden] graph copy done elapsed="
                << ElapsedSeconds(t0, t1)
                << " seconds\n";


    if (options.debug)
        std::cerr << "[Leiden] stats copy start\n";

    const auto t2 = Clock::now();
    LeidenGraphStats current_stats = stats;
    const auto t3 = Clock::now();

    if (options.debug)
        std::cerr << "[Leiden] stats copy done elapsed="
                << ElapsedSeconds(t2, t3)
                << " seconds\n";


    if (options.debug)
        std::cerr << "[Leiden] singleton partition start\n";

    const auto t4 = Clock::now();
    LeidenPartition current_partition =
        MakeSingletonPartition(current_graph, current_stats);
    const auto t5 = Clock::now();

    if (options.debug)
        std::cerr << "[Leiden] singleton partition done elapsed="
                << ElapsedSeconds(t4, t5)
                << " seconds\n";


    if (options.debug)
        std::cerr << "[Leiden] original mapping initialization start\n";

    const auto t6 = Clock::now();
    std::vector<Vertex> original_to_current(n_original);
    std::iota(original_to_current.begin(),
            original_to_current.end(), 0);
    const auto t7 = Clock::now();

    if (options.debug)
        std::cerr << "[Leiden] original mapping initialization done elapsed="
                << ElapsedSeconds(t6, t7)
                << " seconds\n";

#ifndef ENABLE_MOVENODESFAST_STAGE4A
    std::mt19937 rng(options.seed);
#endif
    bool stopped_by_max_levels = false;
    bool converged = false;

    while (options.max_levels == 0 || result.num_levels < options.max_levels) {
        const std::size_t level = result.num_levels + 1;
        const Clock::time_point level_begin = Clock::now();
        const int n_before_aggregation = num_vertices(current_graph);
        if (options.debug) {
            std::cerr << "[Leiden] level " << level << " start\n"
                      << "[Leiden] current number of vertices: "
                      << num_vertices(current_graph) << "\n"
                      << "[Leiden] current number of edges: "
                      << num_edges(current_graph) << "\n";
        }

        if (options.debug) {
            std::cerr << "[Leiden] MoveNodesFast start\n";
        }
        const Clock::time_point move_begin = Clock::now();
        MoveNodesFastResult moved =
#ifdef ENABLE_MOVENODESFAST_STAGE4A
            MoveNodesFastParallelStage4A(
                current_graph,
                current_stats,
                current_partition,
                quality_function,
                static_cast<std::uint64_t>(options.seed),
                static_cast<std::uint64_t>(level),
                &options);
#else
            MoveNodesFast(current_graph,
                          current_stats,
                          current_partition,
                          quality_function,
                          rng,
                          &options);
#endif
        const Clock::time_point move_end = Clock::now();
        result.timing.move_nodes_fast +=
            ElapsedSeconds(move_begin, move_end);
#ifdef ENABLE_MOVENODESFAST_PROFILE
        // TEMPORARY MOVENODESFAST PERFORMANCE PROFILING
        result.move_nodes_fast_profile.neighbor_weights +=
            moved.profile.neighbor_weights;
        result.move_nodes_fast_profile.candidate_build +=
            moved.profile.candidate_build;
        result.move_nodes_fast_profile.delta_evaluation +=
            moved.profile.delta_evaluation;
        result.move_nodes_fast_profile.move_node += moved.profile.move_node;
        result.move_nodes_fast_profile.move_self_loop_scan +=
            moved.profile.move_self_loop_scan;
        result.move_nodes_fast_profile.move_validation +=
            moved.profile.move_validation;
        result.move_nodes_fast_profile.move_ensure_community +=
            moved.profile.move_ensure_community;
        result.move_nodes_fast_profile.move_statistics_update +=
            moved.profile.move_statistics_update;
        result.move_nodes_fast_profile.move_empty_community +=
            moved.profile.move_empty_community;
        result.move_nodes_fast_profile.neighbor_requeue +=
            moved.profile.neighbor_requeue;
        result.move_nodes_fast_profile.num_visits += moved.profile.num_visits;
        result.move_nodes_fast_profile.total_candidates +=
            moved.profile.total_candidates;
        result.move_nodes_fast_profile.stage4a_active_scan +=
            moved.profile.stage4a_active_scan;
        result.move_nodes_fast_profile.stage4a_neighbor_weights +=
            moved.profile.stage4a_neighbor_weights;
        result.move_nodes_fast_profile.stage4a_candidate_evaluation +=
            moved.profile.stage4a_candidate_evaluation;
        result.move_nodes_fast_profile.stage4a_proposal_generation +=
            moved.profile.stage4a_proposal_generation;
        result.move_nodes_fast_profile.stage4a_deterministic_commit +=
            moved.profile.stage4a_deterministic_commit;
        result.move_nodes_fast_profile.stage4a_commit_revalidation +=
            moved.profile.stage4a_commit_revalidation;
        result.move_nodes_fast_profile.stage4a_affected_next_update +=
            moved.profile.stage4a_affected_next_update;
        result.move_nodes_fast_profile.stage4a_buffer_clear +=
            moved.profile.stage4a_buffer_clear;
        result.move_nodes_fast_profile.stage4a_total +=
            moved.profile.stage4a_total;
        result.move_nodes_fast_profile.stage4a_rounds +=
            moved.profile.stage4a_rounds;
        result.move_nodes_fast_profile.stage4a_active_vertices +=
            moved.profile.stage4a_active_vertices;
        result.move_nodes_fast_profile.stage4a_positive_proposals +=
            moved.profile.stage4a_positive_proposals;
        result.move_nodes_fast_profile.stage4a_committed_moves +=
            moved.profile.stage4a_committed_moves;
        result.move_nodes_fast_profile.stage4a_rejected_proposals +=
            moved.profile.stage4a_rejected_proposals;
        result.move_nodes_fast_profile.stage4a_max_active_vertices = std::max(
            result.move_nodes_fast_profile.stage4a_max_active_vertices,
            moved.profile.stage4a_max_active_vertices);
        result.move_nodes_fast_profile.stage4a_reactivation_neighbor_scans +=
            moved.profile.stage4a_reactivation_neighbor_scans;
        result.move_nodes_fast_profile.stage4a_reactivation_target_exclusions +=
            moved.profile.stage4a_reactivation_target_exclusions;
        result.move_nodes_fast_profile.stage4a_reactivation_new_activations +=
            moved.profile.stage4a_reactivation_new_activations;
        result.move_nodes_fast_profile.stage4a_reactivation_duplicate_attempts +=
            moved.profile.stage4a_reactivation_duplicate_attempts;
#endif
        current_partition = std::move(moved.partition);
        ++result.num_levels;
        result.total_moves += moved.num_moves;
        if (options.debug) {
            std::cerr << "[Leiden] MoveNodesFast done\n"
                      << "[Leiden] num_moves: " << moved.num_moves << "\n"
                      << "[Leiden] num_visits: " << moved.num_visits << "\n"
                      << "[Leiden] elapsed time: "
                      << ElapsedSeconds(move_begin, move_end) << " seconds\n";
        }

        if (moved.num_moves == 0 ||
            (options.max_levels != 0 && result.num_levels >= options.max_levels)) {
            converged = (moved.num_moves == 0);
            stopped_by_max_levels =
                !converged &&
                options.max_levels != 0 &&
                result.num_levels >= options.max_levels;
            if (options.debug) {
                std::cerr << "[Leiden] total level time: "
                          << ElapsedSeconds(level_begin, Clock::now())
                          << " seconds\n";
            }
            break;
        }

        if (options.debug) {
            std::cerr << "[Leiden] RefinePartition start\n";
        }
        const Clock::time_point refine_begin = Clock::now();
        LeidenPartition refined =
            RefinePartition(current_graph,
                            current_stats,
                            current_partition,
                            quality_function,
                            options.theta,
                            static_cast<std::uint64_t>(options.seed),
                            level,
                            &options);
        const Clock::time_point refine_end = Clock::now();
        result.timing.refine_partition +=
            ElapsedSeconds(refine_begin, refine_end);
        if (options.debug) {
            std::cerr << "[Leiden] RefinePartition done\n"
                      << "[Leiden] number of refined communities: "
                      << CountActiveCommunities(refined) << "\n"
                      << "[Leiden] elapsed time: "
                      << ElapsedSeconds(refine_begin, refine_end)
                      << " seconds\n";
        }

        if (options.debug) {
            std::cerr << "[Leiden] AggregateGraph start\n";
        }
        const Clock::time_point aggregate_begin = Clock::now();
        AggregateGraphResult aggregate =
            AggregateGraph(current_graph, current_stats, refined);
        const Clock::time_point aggregate_end = Clock::now();
        result.timing.aggregate_graph +=
            ElapsedSeconds(aggregate_begin, aggregate_end);
        if (options.debug) {
            std::cerr << "[Leiden] AggregateGraph done\n"
                      << "[Leiden] aggregate vertices: "
                      << num_vertices(aggregate.graph) << "\n"
                      << "[Leiden] aggregate edges: "
                      << num_edges(aggregate.graph) << "\n"
                      << "[Leiden] elapsed time: "
                      << ElapsedSeconds(aggregate_begin, aggregate_end)
                      << " seconds\n";
        }

        if (options.debug) {
            std::cerr << "[Leiden] BuildCoarsePartition start\n";
        }
        const Clock::time_point coarse_begin = Clock::now();
        LeidenPartition coarse_partition =
            BuildCoarsePartition(aggregate, current_partition, refined);
        const Clock::time_point coarse_end = Clock::now();
        result.timing.build_coarse_partition +=
            ElapsedSeconds(coarse_begin, coarse_end);
        if (options.debug) {
            std::cerr << "[Leiden] BuildCoarsePartition done\n"
                      << "[Leiden] elapsed time: "
                      << ElapsedSeconds(coarse_begin, coarse_end)
                      << " seconds\n";
        }

        if (num_vertices(aggregate.graph) > n_before_aggregation) {
            throw std::logic_error("aggregation increased graph size");
        }

        for (Vertex original = 0; original < n_original; ++original) {
            const Vertex current = original_to_current[original];
            if (current < 0 ||
                static_cast<std::size_t>(current) >= aggregate.coarse_of.size()) {
                throw std::logic_error("original_to_current is out of range");
            }
            original_to_current[original] = aggregate.coarse_of[current];
        }

        current_graph = std::move(aggregate.graph);
        current_stats = std::move(aggregate.stats);
        current_partition = std::move(coarse_partition);

        if (options.debug) {
            std::cerr << "[Leiden] total level time: "
                      << ElapsedSeconds(level_begin, Clock::now())
                      << " seconds\n";
        }
    }

    if (!converged &&
        options.max_levels != 0 &&
        result.num_levels >= options.max_levels) {
        stopped_by_max_levels = true;
    }

    std::vector<Community> final_assignment(n_original, -1);
    for (Vertex original = 0; original < n_original; ++original) {
        const Vertex current = original_to_current[original];
        if (current < 0 ||
            current >= static_cast<Vertex>(current_partition.community_of.size())) {
            throw std::logic_error("final mapping is out of range");
        }
        final_assignment[original] = current_partition.community_of[current];
    }

    result.partition =
        MakePartition(G, stats, CompactAssignment(final_assignment));
    if (options.debug) {
        if (converged) {
            std::cerr << "[Leiden] converged because num_moves == 0\n";
        } else if (stopped_by_max_levels) {
            std::cerr << "[Leiden] stopped because max_levels reached\n";
        } else {
            std::cerr << "[Leiden] converged because num_moves == 0\n";
        }
    }
    result.timing.total =
        ElapsedSeconds(leiden_begin, Clock::now());
    // END TEMPORARY PERFORMANCE INSTRUMENTATION
    return result;
}
