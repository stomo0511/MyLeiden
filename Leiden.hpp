#pragma once

#include <cstddef>
#include <cstdint>
#include <random>
#include <unordered_map>
#include <vector>

#include "QualityFunction.hpp"

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
    std::size_t active_position = static_cast<std::size_t>(-1);
    std::size_t nonpositive_mass_position = static_cast<std::size_t>(-1);
};

struct RefinementCommunityStats {
    std::unordered_map<Community, RefinementCommunityEntry> entries;
    std::vector<Community> active_communities;
    std::vector<Community> nonpositive_mass_active_communities;
};

std::vector<Community> BuildExactSparseRefinementTargets(
    const LeidenGraphStats& stats,
    const QualityFunction& quality_function,
    const RefinementCommunityStats& community_stats,
    const NeighborCommunityScratch& neighbor_scratch,
    Vertex v,
    Community source,
    std::size_t* exceptional_targets_added = nullptr);

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

// Experimental Stage-4A.1 implementation. Proposal generation is parallel and
// read-only; deterministic commit is serial. The production serial reference
// above intentionally remains available.
MoveNodesFastResult MoveNodesFastParallelStage4A(
    const Graph& G,
    const LeidenGraphStats& stats,
    LeidenPartition partition,
    const QualityFunction& quality_function,
    std::uint64_t global_seed,
    std::uint64_t leiden_level,
    const LeidenOptions* options = nullptr);

// Experimental Stage-4B race-free asynchronous local-moving variant. It is
// intentionally not the sequential Algorithm A.2 execution order.
MoveNodesFastResult MoveNodesFastParallelStage4B(
    const Graph& G,
    const LeidenGraphStats& stats,
    LeidenPartition partition,
    const QualityFunction& quality_function,
    std::uint64_t global_seed,
    std::uint64_t leiden_level,
    const LeidenOptions* options = nullptr);

// Observable result of the Stage-4A.1 serial affected-set update. Kept
// independent of profiling macros so targeted correctness tests can exercise
// the exact production reactivation rule.
struct Stage4A1ReactivationStats {
    std::size_t neighbor_scans = 0;
    std::size_t target_community_exclusions = 0;
    std::size_t newly_activated = 0;
    std::size_t duplicate_attempts = 0;
};

Stage4A1ReactivationStats UpdateStage4A1AffectedNeighbors(
    const Graph& G,
    const LeidenPartition& committed_partition,
    Vertex moved_vertex,
    Community committed_target,
    std::vector<unsigned char>& affected_next);

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

std::uint64_t MakeSubsetSeed(std::uint64_t global_seed,
                             std::uint64_t level,
                             std::uint64_t parent);

LeidenPartition RefinePartition(const Graph& G,
                                const LeidenGraphStats& stats,
                                const LeidenPartition& partition,
                                const QualityFunction& quality_function,
                                double theta,
                                std::uint64_t global_seed,
                                std::uint64_t level,
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
