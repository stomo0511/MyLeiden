#include "Leiden.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <deque>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

struct EdgePairHash {
    std::size_t operator()(const std::pair<int, int>& p) const
    {
        const std::size_t a = std::hash<int>{}(p.first);
        const std::size_t b = std::hash<int>{}(p.second);
        return a ^ (b + 0x9e3779b97f4a7c15ULL + (a << 6U) + (a >> 2U));
    }
};

struct WeightedEdge {
    int u;
    int v;
    double weight;
};

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

void MarkSubset(const Graph& G,
                const std::vector<Vertex>& subset,
                std::vector<std::size_t>& subset_mark,
                std::size_t& subset_generation)
{
    if (subset_mark.size() != static_cast<std::size_t>(num_vertices(G))) {
        throw std::invalid_argument("subset mark size does not match graph");
    }

    if (subset_generation == std::numeric_limits<std::size_t>::max()) {
        std::fill(subset_mark.begin(), subset_mark.end(), 0);
        subset_generation = 1;
    } else {
        ++subset_generation;
    }

    for (Vertex v : subset) {
        ValidateVertex(v, G);
        subset_mark[v] = subset_generation;
    }
}

std::vector<Community> BuildCandidateCommunities(
    const LeidenPartition& partition,
    const std::unordered_map<Community, double>& neighbor_weights)
{
    std::vector<Community> candidates;
    candidates.reserve(neighbor_weights.size() + 1);

    for (const auto& item : neighbor_weights) {
        if (item.first >= 0) {
            candidates.push_back(item.first);
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
#else
#define MNF_PROFILE_BEGIN(name)
#define MNF_PROFILE_END(name, field)
#define MNF_PROFILE_CANDIDATES(count)
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
    const bool debug = DebugEnabled(options);
    const std::size_t debug_interval = DebugInterval(options);

    while (!queue.empty()) {
        const Vertex v = queue.front();
        queue.pop_front();
        in_queue[v] = false;
        ++result.num_visits;

        const Community source = result.partition.community_of[v];
        MNF_PROFILE_BEGIN(neighbor_weights);
        const auto neighbor_weights =
            BuildNeighborCommunityWeights(G, result.partition, v);
        MNF_PROFILE_END(neighbor_weights, neighbor_weights);
        const double weight_to_source = LookupWeight(neighbor_weights, source);
        MNF_PROFILE_BEGIN(candidate_build);
        const std::vector<Community> candidates =
            BuildCandidateCommunities(result.partition, neighbor_weights);
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
                          LookupWeight(neighbor_weights, candidate));

            if (delta > best_delta ||
                (delta == best_delta && candidate < best_community)) {
                best_delta = delta;
                best_community = candidate;
            }
        }
        MNF_PROFILE_END(delta_evaluation, delta_evaluation);

        if (best_delta > 0.0 && best_community != source) {
            MNF_PROFILE_BEGIN(move_node);
            MoveNodeToCommunity(G, stats, result.partition, v, best_community);
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

#undef MNF_PROFILE_BEGIN
#undef MNF_PROFILE_END
#undef MNF_PROFILE_CANDIDATES
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

LeidenPartition RefinePartition(const Graph& G,
                                const LeidenGraphStats& stats,
                                const LeidenPartition& partition,
                                const QualityFunction& quality_function,
                                double theta,
                                std::mt19937& rng,
                                const LeidenOptions* options)
{
    if (theta <= 0.0) {
        throw std::invalid_argument("theta must be positive");
    }

    LeidenPartition refined = MakeSingletonPartition(G, stats);
    const std::vector<std::vector<Vertex>> members =
        BuildPartitionMembers(partition);
    const bool debug = DebugEnabled(options);
    const std::size_t debug_interval = DebugInterval(options);
    const std::size_t total_parent_communities =
        static_cast<std::size_t>(
            std::count_if(members.begin(),
                          members.end(),
                          [](const std::vector<Vertex>& subset) {
                              return !subset.empty();
                          }));
    std::size_t processed_parent_communities = 0;
    std::size_t processed_vertices = 0;
    std::size_t next_parent_report = debug_interval;
    std::size_t next_vertex_report = debug_interval;
    std::vector<std::size_t> subset_mark(num_vertices(G), 0);
    std::size_t subset_generation = 0;

    for (const std::vector<Vertex>& subset : members) {
        if (!subset.empty()) {
            MarkSubset(G, subset, subset_mark, subset_generation);
            MergeNodesSubset(G,
                             stats,
                             refined,
                             subset,
                             quality_function,
                             theta,
                             rng,
                             subset_mark,
                             subset_generation,
                             options);
            ++processed_parent_communities;
            processed_vertices += subset.size();
            if (debug &&
                (processed_parent_communities >= next_parent_report ||
                 processed_vertices >= next_vertex_report)) {
                std::cerr << "[RefinePartition] processed_parent_communities="
                          << processed_parent_communities
                          << " total_parent_communities="
                          << total_parent_communities
                          << " processed_vertices=" << processed_vertices
                          << " current_subset_size=" << subset.size()
                          << "\n";
                while (processed_parent_communities >= next_parent_report) {
                    next_parent_report += debug_interval;
                }
                while (processed_vertices >= next_vertex_report) {
                    next_vertex_report += debug_interval;
                }
            }
        }
    }
    if (debug) {
        std::cerr << "[RefinePartition] done processed_parent_communities="
                  << processed_parent_communities
                  << " total_parent_communities=" << total_parent_communities
                  << " processed_vertices=" << processed_vertices
                  << "\n";
    }
    return refined;
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
    std::vector<double> coarse_node_size(n_coarse, 0.0);

    for (Vertex v = 0; v < num_vertices(G); ++v) {
        const Community coarse = community_to_coarse[refined.community_of[v]];
        result.coarse_of[v] = coarse;
        coarse_node_size[coarse] += stats.node_size[v];
    }

    std::unordered_map<std::pair<int, int>, double, EdgePairHash> edge_weight;
    edge_weight.reserve(num_edges(G));

    // Each undirected input edge is processed once. Edges inside a refined
    // community become aggregate self-loops; inter-community and parallel
    // edges are summed by canonical (min,max) coarse pair.
    for_each_undirected_edge(G, [&](int u, int v, double w) {
        int cu = result.coarse_of[u];
        int cv = result.coarse_of[v];
        if (cu > cv) {
            std::swap(cu, cv);
        }
        edge_weight[std::make_pair(cu, cv)] += w;
    });

    std::vector<WeightedEdge> edges;
    edges.reserve(edge_weight.size());
    for (const auto& item : edge_weight) {
        edges.push_back({item.first.first, item.first.second, item.second});
    }
    std::sort(edges.begin(), edges.end(), [](const WeightedEdge& a,
                                             const WeightedEdge& b) {
        if (a.u != b.u) {
            return a.u < b.u;
        }
        return a.v < b.v;
    });

    for (const WeightedEdge& e : edges) {
        add_undirected_edge(result.graph, e.u, e.v, e.weight);
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

    std::mt19937 rng(options.seed);
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
            MoveNodesFast(current_graph,
                          current_stats,
                          current_partition,
                          quality_function,
                          rng,
                          &options);
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
        result.move_nodes_fast_profile.neighbor_requeue +=
            moved.profile.neighbor_requeue;
        result.move_nodes_fast_profile.num_visits += moved.profile.num_visits;
        result.move_nodes_fast_profile.total_candidates +=
            moved.profile.total_candidates;
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
                            rng,
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
