#pragma once

#include <cstddef>
#include <random>
#include <unordered_map>
#include <vector>

#include "QualityFunction.hpp"

// -----------------------------------------------------------------------------
// TEMPORARY MOVENODESFAST PERFORMANCE PROFILING
// Detailed profiling for MoveNodesFast(). Remove this block after performance
// analysis. It is compiled only when ENABLE_MOVENODESFAST_PROFILE is defined.
// -----------------------------------------------------------------------------
#ifdef ENABLE_MOVENODESFAST_PROFILE
struct MoveNodesFastProfile {
    double neighbor_weights = 0.0;
    double candidate_build = 0.0;
    double delta_evaluation = 0.0;
    double move_node = 0.0;
    double neighbor_requeue = 0.0;
    std::size_t num_visits = 0;
    std::size_t total_candidates = 0;
};
#endif
// END TEMPORARY MOVENODESFAST PERFORMANCE PROFILING

struct MoveNodesFastResult {
    LeidenPartition partition;
    std::size_t num_moves = 0;
    std::size_t num_visits = 0;
#ifdef ENABLE_MOVENODESFAST_PROFILE
    // TEMPORARY MOVENODESFAST PERFORMANCE PROFILING
    MoveNodesFastProfile profile;
#endif
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

// -----------------------------------------------------------------------------
// TEMPORARY PERFORMANCE INSTRUMENTATION
// These timers are used only for profiling the current Leiden implementation.
// Keep this block isolated so it can be removed easily after performance study.
// -----------------------------------------------------------------------------
struct LeidenTiming {
    double move_nodes_fast = 0.0;
    double refine_partition = 0.0;
    double aggregate_graph = 0.0;
    double build_coarse_partition = 0.0;
    double total = 0.0;
};
// END TEMPORARY PERFORMANCE INSTRUMENTATION

struct LeidenResult {
    LeidenPartition partition;
    std::size_t num_levels = 0;
    std::size_t total_moves = 0;
    // TEMPORARY PERFORMANCE INSTRUMENTATION
    LeidenTiming timing;
#ifdef ENABLE_MOVENODESFAST_PROFILE
    // TEMPORARY MOVENODESFAST PERFORMANCE PROFILING
    MoveNodesFastProfile move_nodes_fast_profile;
#endif
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
