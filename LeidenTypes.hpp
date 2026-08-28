#pragma once

#include <algorithm>
#include <set>
#include <stdexcept>
#include <vector>

#include "common/Types.hpp"

using Community = int;

struct LeidenGraphStats {
    std::vector<double> node_size;
    std::vector<double> node_strength;
    double total_edge_weight = 0.0;
};

struct LeidenPartition {
    std::vector<Community> community_of;
    std::vector<double> community_size;
    std::vector<double> community_strength;
    std::vector<double> internal_edge_weight;
    std::set<Community> empty_communities;
};

inline Community EmptyCommunityForMove(const LeidenPartition& partition)
{
    if (!partition.empty_communities.empty()) {
        return *partition.empty_communities.begin();
    }
    return static_cast<Community>(partition.community_size.size());
}

inline void ValidateVertex(Vertex v, const Graph& G)
{
    if (v < 0 || v >= num_vertices(G)) {
        throw std::out_of_range("vertex id out of range");
    }
}

inline int NumCommunitiesFromAssignment(const std::vector<Community>& community_of)
{
    if (community_of.empty()) {
        return 0;
    }
    const auto max_it = std::max_element(community_of.begin(), community_of.end());
    return (*max_it < 0) ? 0 : (*max_it + 1);
}

LeidenGraphStats BuildLeidenGraphStats(const Graph& G);
LeidenGraphStats BuildLeidenGraphStats(const Graph& G,
                                       const std::vector<double>& node_size);

LeidenPartition MakeSingletonPartition(const Graph& G,
                                       const LeidenGraphStats& stats);
LeidenPartition MakePartition(const Graph& G,
                              const LeidenGraphStats& stats,
                              const std::vector<Community>& community_of);

double WeightFromNodeToCommunity(const Graph& G,
                                 const LeidenPartition& partition,
                                 Vertex v,
                                 Community community,
                                 bool include_self_loop);

double SelfLoopWeight(const Graph& G, Vertex v);

void RemoveNodeFromCommunity(const Graph& G,
                             const LeidenGraphStats& stats,
                             LeidenPartition& partition,
                             Vertex v);

void InsertNodeIntoCommunity(const Graph& G,
                             const LeidenGraphStats& stats,
                             LeidenPartition& partition,
                             Vertex v,
                             Community community);

void MoveNodeToCommunity(const Graph& G,
                         const LeidenGraphStats& stats,
                         LeidenPartition& partition,
                         Vertex v,
                         Community community);

// Optimized MoveNodesFast path: reuse precomputed neighbor-community weights
// to avoid rescanning G.adj[v] for the source and target communities.
void MoveNodeToCommunityFromWeights(const Graph& G,
                                    const LeidenGraphStats& stats,
                                    LeidenPartition& partition,
                                    Vertex v,
                                    Community community,
                                    double weight_to_source,
                                    double weight_to_target,
                                    double self_loop_weight);
