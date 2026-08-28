#pragma once

#include <cstddef>
#include <random>
#include <unordered_map>
#include <vector>

#include "QualityFunction.hpp"

struct MoveNodesFastResult {
    LeidenPartition partition;
    std::size_t num_moves = 0;
    std::size_t num_visits = 0;
};

struct RefinementCommunityEntry {
    int member_count = 0;
    double mass = 0.0;
    double external_weight = 0.0;
};

struct RefinementCommunityStats {
    std::unordered_map<Community, RefinementCommunityEntry> entries;
    std::vector<Community> active_communities;
};

struct AggregateGraphResult {
    Graph graph;
    LeidenGraphStats stats;
    std::vector<Vertex> coarse_of;
};

struct LeidenOptions {
    double theta = 0.01;
    unsigned int seed = 0;
    std::size_t max_levels = 0;
    bool debug = false;
    std::size_t debug_interval = 100000;
};

struct LeidenResult {
    LeidenPartition partition;
    std::size_t num_levels = 0;
    std::size_t total_moves = 0;
};

MoveNodesFastResult MoveNodesFast(const Graph& G,
                                  const LeidenGraphStats& stats,
                                  LeidenPartition partition,
                                  const QualityFunction& quality_function,
                                  std::mt19937& rng,
                                  const LeidenOptions* options = nullptr);

double EdgeWeightFromNodeToSubset(const Graph& G,
                                  Vertex v,
                                  const std::vector<std::size_t>& subset_mark,
                                  std::size_t subset_generation);

bool IsNodeWellConnectedToSubset(const Graph& G,
                                 const LeidenGraphStats& stats,
                                 const QualityFunction& quality_function,
                                 Vertex v,
                                 const std::vector<std::size_t>& subset_mark,
                                 std::size_t subset_generation,
                                 double subset_mass);

bool IsCommunityWellConnectedToSubset(
    const Graph& G,
    const LeidenGraphStats& stats,
    const LeidenPartition& refined,
    const QualityFunction& quality_function,
    Community community,
    const std::vector<Vertex>& subset,
    const std::vector<std::size_t>& subset_mark,
    std::size_t subset_generation);

RefinementCommunityStats BuildRefinementCommunityStats(
    const Graph& G,
    const LeidenGraphStats& stats,
    const LeidenPartition& refined,
    const QualityFunction& quality_function,
    const std::vector<Vertex>& subset,
    const std::vector<std::size_t>& subset_mark,
    std::size_t subset_generation);

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
    const std::vector<std::size_t>& subset_mark,
    std::size_t subset_generation,
    Vertex v,
    Community target,
    RefinementCommunityStats& community_stats);

void MergeNodesSubset(const Graph& G,
                      const LeidenGraphStats& stats,
                      LeidenPartition& refined,
                      const std::vector<Vertex>& subset,
                      const QualityFunction& quality_function,
                      double theta,
                      std::mt19937& rng,
                      const std::vector<std::size_t>& subset_mark,
                      std::size_t subset_generation,
                      const LeidenOptions* options = nullptr);

LeidenPartition RefinePartition(const Graph& G,
                                const LeidenGraphStats& stats,
                                const LeidenPartition& partition,
                                const QualityFunction& quality_function,
                                double theta,
                                std::mt19937& rng,
                                const LeidenOptions* options = nullptr);

AggregateGraphResult AggregateGraph(const Graph& G,
                                    const LeidenGraphStats& stats,
                                    const LeidenPartition& refined);

LeidenPartition BuildCoarsePartition(const AggregateGraphResult& aggregate,
                                     const LeidenPartition& partition,
                                     const LeidenPartition& refined);

LeidenResult Leiden(const Graph& G,
                    const LeidenGraphStats& stats,
                    const QualityFunction& quality_function,
                    const LeidenOptions& options);
