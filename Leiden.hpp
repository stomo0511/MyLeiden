#pragma once

#include <cstddef>
#include <random>
#include <vector>

#include "QualityFunction.hpp"

struct MoveNodesFastResult {
    LeidenPartition partition;
    std::size_t num_moves = 0;
    std::size_t num_visits = 0;
};

struct RefinementCommunityStats {
    std::vector<int> member_count;
    std::vector<double> mass;
    std::vector<double> external_weight;
    std::vector<Community> active_communities;
};

struct AggregateGraphResult {
    Graph graph;
    LeidenGraphStats stats;
    std::vector<Vertex> coarse_of;
};

MoveNodesFastResult MoveNodesFast(const Graph& G,
                                  const LeidenGraphStats& stats,
                                  LeidenPartition partition,
                                  const QualityFunction& quality_function,
                                  std::mt19937& rng);

double EdgeWeightFromNodeToSubset(const Graph& G,
                                  Vertex v,
                                  const std::vector<bool>& in_subset);

bool IsNodeWellConnectedToSubset(const Graph& G,
                                 const LeidenGraphStats& stats,
                                 const QualityFunction& quality_function,
                                 Vertex v,
                                 const std::vector<bool>& in_subset,
                                 double subset_mass);

bool IsCommunityWellConnectedToSubset(
    const Graph& G,
    const LeidenGraphStats& stats,
    const LeidenPartition& refined,
    const QualityFunction& quality_function,
    Community community,
    const std::vector<Vertex>& subset,
    const std::vector<bool>& in_subset);

RefinementCommunityStats BuildRefinementCommunityStats(
    const Graph& G,
    const LeidenGraphStats& stats,
    const LeidenPartition& refined,
    const QualityFunction& quality_function,
    const std::vector<Vertex>& subset,
    const std::vector<bool>& in_subset);

bool IsCommunityWellConnectedFromStats(
    const LeidenGraphStats& stats,
    const QualityFunction& quality_function,
    Community community,
    double subset_mass,
    const RefinementCommunityStats& community_stats);

void UpdateRefinementCommunityStatsForMove(
    const Graph& G,
    const LeidenGraphStats& stats,
    const LeidenPartition& refined_before_move,
    const QualityFunction& quality_function,
    const std::vector<bool>& in_subset,
    Vertex v,
    Community target,
    RefinementCommunityStats& community_stats);

void MergeNodesSubset(const Graph& G,
                      const LeidenGraphStats& stats,
                      LeidenPartition& refined,
                      const std::vector<Vertex>& subset,
                      const QualityFunction& quality_function,
                      double theta,
                      std::mt19937& rng);

LeidenPartition RefinePartition(const Graph& G,
                                const LeidenGraphStats& stats,
                                const LeidenPartition& partition,
                                const QualityFunction& quality_function,
                                double theta,
                                std::mt19937& rng);

AggregateGraphResult AggregateGraph(const Graph& G,
                                    const LeidenGraphStats& stats,
                                    const LeidenPartition& refined);
