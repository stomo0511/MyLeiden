#include "Leiden.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <functional>
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

double LookupWeight(const std::unordered_map<Community, double>& weights,
                    Community community)
{
    const auto it = weights.find(community);
    return (it == weights.end()) ? 0.0 : it->second;
}

Community EmptyCommunityForMove(const LeidenPartition& partition)
{
    for (Community c = 0;
         c < static_cast<Community>(partition.community_size.size());
         ++c) {
        if (partition.community_size[c] == 0.0) {
            return c;
        }
    }
    return static_cast<Community>(partition.community_size.size());
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

double RefinedCommunityMassInSubset(const LeidenGraphStats& stats,
                                    const LeidenPartition& refined,
                                    const QualityFunction& quality_function,
                                    Community community,
                                    const std::vector<Vertex>& subset)
{
    double mass = 0.0;
    for (Vertex v : subset) {
        if (refined.community_of[v] == community) {
            mass += quality_function.refinementNodeMass(stats, v);
        }
    }
    return mass;
}

void EnsureRefinementCommunityStatsSize(RefinementCommunityStats& community_stats,
                                        Community community)
{
    if (community < 0) {
        throw std::out_of_range("community id out of range");
    }
    const std::size_t need = static_cast<std::size_t>(community) + 1;
    if (community_stats.member_count.size() < need) {
        community_stats.member_count.resize(need, 0);
        community_stats.mass.resize(need, 0.0);
        community_stats.external_weight.resize(need, 0.0);
    }
}

void RefreshActiveCommunities(RefinementCommunityStats& community_stats)
{
    community_stats.active_communities.clear();
    for (Community c = 0;
         c < static_cast<Community>(community_stats.member_count.size());
         ++c) {
        if (community_stats.member_count[c] > 0) {
            community_stats.active_communities.push_back(c);
        }
    }
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

std::vector<bool> BuildSubsetMask(const Graph& G,
                                  const std::vector<Vertex>& subset)
{
    std::vector<bool> in_subset(num_vertices(G), false);
    for (Vertex v : subset) {
        ValidateVertex(v, G);
        in_subset[v] = true;
    }
    return in_subset;
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
                                  const std::vector<bool>& in_subset)
{
    ValidateVertex(v, G);
    if (static_cast<int>(in_subset.size()) != num_vertices(G)) {
        throw std::invalid_argument("subset mask size does not match graph");
    }

    double weight = 0.0;
    for (const Edge& e : G.adj[v]) {
        if (e.to != v && in_subset[e.to]) {
            weight += e.weight;
        }
    }
    return weight;
}

bool IsNodeWellConnectedToSubset(const Graph& G,
                                 const LeidenGraphStats& stats,
                                 const QualityFunction& quality_function,
                                 Vertex v,
                                 const std::vector<bool>& in_subset,
                                 double subset_mass)
{
    const double node_mass = quality_function.refinementNodeMass(stats, v);
    const double edge_weight = EdgeWeightFromNodeToSubset(G, v, in_subset);
    const double threshold =
        quality_function.refinementResolution(stats) *
        node_mass *
        (subset_mass - node_mass);
    return edge_weight >= threshold;
}

bool IsCommunityWellConnectedToSubset(
    const Graph& G,
    const LeidenGraphStats& stats,
    const LeidenPartition& refined,
    const QualityFunction& quality_function,
    Community community,
    const std::vector<Vertex>& subset,
    const std::vector<bool>& in_subset)
{
    if (community < 0 ||
        static_cast<std::size_t>(community) >= refined.community_size.size()) {
        return false;
    }

    double edge_weight = 0.0;
    for (Vertex v : subset) {
        if (refined.community_of[v] != community) {
            continue;
        }
        for (const Edge& e : G.adj[v]) {
            if (e.to != v && in_subset[e.to] &&
                refined.community_of[e.to] != community) {
                edge_weight += e.weight;
            }
        }
    }

    const double subset_mass = SubsetMass(stats, quality_function, subset);
    const double community_mass =
        RefinedCommunityMassInSubset(stats,
                                     refined,
                                     quality_function,
                                     community,
                                     subset);
    const double threshold =
        quality_function.refinementResolution(stats) *
        community_mass *
        (subset_mass - community_mass);
    return edge_weight >= threshold;
}

RefinementCommunityStats BuildRefinementCommunityStats(
    const Graph& G,
    const LeidenGraphStats& stats,
    const LeidenPartition& refined,
    const QualityFunction& quality_function,
    const std::vector<Vertex>& subset,
    const std::vector<bool>& in_subset)
{
    if (static_cast<int>(in_subset.size()) != num_vertices(G)) {
        throw std::invalid_argument("subset mask size does not match graph");
    }

    RefinementCommunityStats community_stats;
    community_stats.member_count.assign(refined.community_size.size(), 0);
    community_stats.mass.assign(refined.community_size.size(), 0.0);
    community_stats.external_weight.assign(refined.community_size.size(), 0.0);

    for (Vertex v : subset) {
        ValidateVertex(v, G);
        const Community community = refined.community_of[v];
        EnsureRefinementCommunityStatsSize(community_stats, community);
        ++community_stats.member_count[community];
        community_stats.mass[community] +=
            quality_function.refinementNodeMass(stats, v);
    }

    for (Vertex u : subset) {
        ValidateVertex(u, G);
        const Community cu = refined.community_of[u];
        EnsureRefinementCommunityStatsSize(community_stats, cu);

        for (const Edge& e : G.adj[u]) {
            if (e.to == u || !in_subset[e.to] || u > e.to) {
                continue;
            }

            const Community cv = refined.community_of[e.to];
            EnsureRefinementCommunityStatsSize(community_stats, cv);
            if (cu != cv) {
                community_stats.external_weight[cu] += e.weight;
                community_stats.external_weight[cv] += e.weight;
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
    if (community < 0 ||
        static_cast<std::size_t>(community) >= community_stats.member_count.size() ||
        community_stats.member_count[community] <= 0) {
        return false;
    }

    const double community_mass = community_stats.mass[community];
    const double threshold =
        quality_function.refinementResolution(stats) *
        community_mass *
        (subset_mass - community_mass);
    return community_stats.external_weight[community] >= threshold;
}

void UpdateRefinementCommunityStatsForMove(
    const Graph& G,
    const LeidenGraphStats& stats,
    const LeidenPartition& refined_before_move,
    const QualityFunction& quality_function,
    const std::vector<bool>& in_subset,
    Vertex v,
    Community target,
    RefinementCommunityStats& community_stats)
{
    ValidateVertex(v, G);
    if (static_cast<int>(in_subset.size()) != num_vertices(G)) {
        throw std::invalid_argument("subset mask size does not match graph");
    }

    const Community source = refined_before_move.community_of[v];
    if (source == target) {
        return;
    }

    EnsureRefinementCommunityStatsSize(community_stats, source);
    EnsureRefinementCommunityStatsSize(community_stats, target);

    const double node_mass =
        quality_function.refinementNodeMass(stats, v);
    --community_stats.member_count[source];
    ++community_stats.member_count[target];
    community_stats.mass[source] -= node_mass;
    community_stats.mass[target] += node_mass;

    // For each subset-internal edge (v,u), update E(C,S-C) using communities
    // before the move A -> B:
    //   C == A: internal A edge becomes A-B external, so A += w and B += w.
    //   C == B: A-B external becomes internal B edge, so A -= w and B -= w.
    //   otherwise: A-C external becomes B-C external, so A -= w and B += w.
    // The third community C is external both before and after, so its value is
    // unchanged. Self-loops and subset-external edges are ignored.
    for (const Edge& e : G.adj[v]) {
        const Vertex u = e.to;
        if (u == v || !in_subset[u]) {
            continue;
        }

        const Community neighbor_community = refined_before_move.community_of[u];
        if (neighbor_community == source) {
            community_stats.external_weight[source] += e.weight;
            community_stats.external_weight[target] += e.weight;
        } else if (neighbor_community == target) {
            community_stats.external_weight[source] -= e.weight;
            community_stats.external_weight[target] -= e.weight;
        } else {
            community_stats.external_weight[source] -= e.weight;
            community_stats.external_weight[target] += e.weight;
        }
    }

    RefreshActiveCommunities(community_stats);
}

MoveNodesFastResult MoveNodesFast(const Graph& G,
                                  const LeidenGraphStats& stats,
                                  LeidenPartition partition,
                                  const QualityFunction& quality_function,
                                  std::mt19937& rng)
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

    while (!queue.empty()) {
        const Vertex v = queue.front();
        queue.pop_front();
        in_queue[v] = false;
        ++result.num_visits;

        const Community source = result.partition.community_of[v];
        const auto neighbor_weights =
            BuildNeighborCommunityWeights(G, result.partition, v);
        const double weight_to_source = LookupWeight(neighbor_weights, source);
        const std::vector<Community> candidates =
            BuildCandidateCommunities(result.partition, neighbor_weights);

        Community best_community = source;
        double best_delta = 0.0;

        // Deterministic tie-breaking: if deltas are exactly equal, choose
        // the smaller community id. The initial node order remains randomized
        // by the caller-supplied rng.
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

        if (best_delta > 0.0 && best_community != source) {
            MoveNodeToCommunity(G, stats, result.partition, v, best_community);
            ++result.num_moves;

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
        }
    }

    return result;
}

void MergeNodesSubset(const Graph& G,
                      const LeidenGraphStats& stats,
                      LeidenPartition& refined,
                      const std::vector<Vertex>& subset,
                      const QualityFunction& quality_function,
                      double theta,
                      std::mt19937& rng)
{
    if (theta <= 0.0) {
        throw std::invalid_argument("theta must be positive");
    }
    if (subset.empty()) {
        return;
    }

    const std::vector<bool> in_subset = BuildSubsetMask(G, subset);
    const double subset_mass = SubsetMass(stats, quality_function, subset);
    std::vector<Vertex> candidates;
    candidates.reserve(subset.size());
    for (Vertex v : subset) {
        if (IsNodeWellConnectedToSubset(G,
                                        stats,
                                        quality_function,
                                        v,
                                        in_subset,
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
                                      in_subset);
    for (Vertex v : candidates) {
        const Community source = refined.community_of[v];
        if (source < 0 ||
            static_cast<std::size_t>(source) >= community_stats.member_count.size() ||
            community_stats.member_count[source] != 1) {
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
                                                  in_subset,
                                                  v,
                                                  target,
                                                  community_stats);
            MoveNodeToCommunity(G, stats, refined, v, target);
        }
    }
}

LeidenPartition RefinePartition(const Graph& G,
                                const LeidenGraphStats& stats,
                                const LeidenPartition& partition,
                                const QualityFunction& quality_function,
                                double theta,
                                std::mt19937& rng)
{
    if (theta <= 0.0) {
        throw std::invalid_argument("theta must be positive");
    }

    LeidenPartition refined = MakeSingletonPartition(G, stats);
    const std::vector<std::vector<Vertex>> members =
        BuildPartitionMembers(partition);
    for (const std::vector<Vertex>& subset : members) {
        if (!subset.empty()) {
            MergeNodesSubset(G,
                             stats,
                             refined,
                             subset,
                             quality_function,
                             theta,
                             rng);
        }
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

    LeidenResult result;
    const int n_original = num_vertices(G);
    if (n_original == 0) {
        result.partition = MakePartition(G, stats, {});
        return result;
    }

    Graph current_graph = G;
    LeidenGraphStats current_stats = stats;
    LeidenPartition current_partition =
        MakeSingletonPartition(current_graph, current_stats);

    // original_to_current[v] is the current coarse-graph vertex that
    // represents original vertex v at the start of each level.
    std::vector<Vertex> original_to_current(n_original);
    std::iota(original_to_current.begin(), original_to_current.end(), 0);

    std::mt19937 rng(options.seed);

    while (options.max_levels == 0 || result.num_levels < options.max_levels) {
        const int n_before_aggregation = num_vertices(current_graph);
        MoveNodesFastResult moved =
            MoveNodesFast(current_graph,
                          current_stats,
                          current_partition,
                          quality_function,
                          rng);
        current_partition = std::move(moved.partition);
        ++result.num_levels;
        result.total_moves += moved.num_moves;

        if (moved.num_moves == 0 ||
            (options.max_levels != 0 && result.num_levels >= options.max_levels)) {
            break;
        }

        LeidenPartition refined =
            RefinePartition(current_graph,
                            current_stats,
                            current_partition,
                            quality_function,
                            options.theta,
                            rng);
        AggregateGraphResult aggregate =
            AggregateGraph(current_graph, current_stats, refined);
        LeidenPartition coarse_partition =
            BuildCoarsePartition(aggregate, current_partition, refined);

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
    return result;
}
