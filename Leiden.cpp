#include "Leiden.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <stdexcept>
#include <numeric>
#include <unordered_map>
#include <vector>

namespace {

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

std::vector<int> CountRefinedCommunityMembers(const LeidenPartition& refined)
{
    std::vector<int> member_count(refined.community_size.size(), 0);
    for (Community community : refined.community_of) {
        if (community >= 0) {
            if (static_cast<std::size_t>(community) >= member_count.size()) {
                member_count.resize(static_cast<std::size_t>(community) + 1, 0);
            }
            ++member_count[community];
        }
    }
    return member_count;
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

    std::vector<int> member_count = CountRefinedCommunityMembers(refined);
    for (Vertex v : candidates) {
        const Community source = refined.community_of[v];
        if (source < 0 ||
            static_cast<std::size_t>(source) >= member_count.size() ||
            member_count[source] != 1) {
            continue;
        }

        const auto neighbor_weights =
            BuildNeighborCommunityWeights(G, refined, v);
        const double weight_to_source = LookupWeight(neighbor_weights, source);

        std::vector<Community> target_communities;
        target_communities.reserve(subset.size());
        for (Vertex u : subset) {
            const Community community = refined.community_of[u];
            if (community >= 0 &&
                IsCommunityWellConnectedToSubset(G,
                                                 stats,
                                                 refined,
                                                 quality_function,
                                                 community,
                                                 subset,
                                                 in_subset)) {
                target_communities.push_back(community);
            }
        }
        std::sort(target_communities.begin(), target_communities.end());
        target_communities.erase(
            std::unique(target_communities.begin(), target_communities.end()),
            target_communities.end());

        struct WeightedCandidate {
            Community community;
            double delta;
        };
        std::vector<WeightedCandidate> nondecreasing_candidates;
        nondecreasing_candidates.reserve(target_communities.size());

        for (Community community : target_communities) {
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
            MoveNodeToCommunity(G, stats, refined, v, target);
            --member_count[source];
            if (static_cast<std::size_t>(target) >= member_count.size()) {
                member_count.resize(static_cast<std::size_t>(target) + 1, 0);
            }
            ++member_count[target];
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
