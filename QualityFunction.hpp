#pragma once

#include "LeidenTypes.hpp"

#include <cstddef>
#include <unordered_map>

#ifdef ENABLE_MOVENODESFAST_PROFILE
// -----------------------------------------------------------------------------
// TEMPORARY MOVENODESFAST PERFORMANCE PROFILING
// Detailed profiling for MoveNodesFast(). Remove this block after performance
// analysis. It is compiled only when ENABLE_MOVENODESFAST_PROFILE is defined.
// -----------------------------------------------------------------------------
struct MoveNodesFastProfile {
    double neighbor_weights = 0.0;
    double candidate_build = 0.0;
    double delta_evaluation = 0.0;
    double move_node = 0.0;

    double move_self_loop_scan = 0.0;
    double move_validation = 0.0;
    double move_ensure_community = 0.0;
    double move_statistics_update = 0.0;
    double move_empty_community = 0.0;

    double neighbor_requeue = 0.0;
    std::size_t num_visits = 0;
    std::size_t total_candidates = 0;

    // TEMPORARY MOVENODESFAST STAGE4A1 PROFILING
    double stage4a_active_scan = 0.0;
    double stage4a_neighbor_weights = 0.0;
    double stage4a_candidate_evaluation = 0.0;
    double stage4a_proposal_generation = 0.0;
    double stage4a_deterministic_commit = 0.0;
    double stage4a_commit_revalidation = 0.0;
    double stage4a_affected_next_update = 0.0;
    double stage4a_buffer_clear = 0.0;
    double stage4a_total = 0.0;
    std::size_t stage4a_rounds = 0;
    std::size_t stage4a_active_vertices = 0;
    std::size_t stage4a_positive_proposals = 0;
    std::size_t stage4a_committed_moves = 0;
    std::size_t stage4a_rejected_proposals = 0;
    std::size_t stage4a_max_active_vertices = 0;
    std::size_t stage4a_reactivation_neighbor_scans = 0;
    std::size_t stage4a_reactivation_target_exclusions = 0;
    std::size_t stage4a_reactivation_new_activations = 0;
    std::size_t stage4a_reactivation_duplicate_attempts = 0;

    // TEMPORARY MOVENODESFAST STAGE4B PROFILING
    double stage4b_neighbor_weights = 0.0;
    double stage4b_candidate_evaluation = 0.0;
    double stage4b_lock_wait = 0.0;
    double stage4b_commit_revalidation = 0.0;
    double stage4b_affected_update = 0.0;
    double stage4b_queue_management = 0.0;
    double stage4b_final_rebuild = 0.0;
    double stage4b_total = 0.0;
    std::size_t stage4b_processed_vertices = 0;
    std::size_t stage4b_successful_moves = 0;
    std::size_t stage4b_failed_validations = 0;
    std::size_t stage4b_neighbor_scans = 0;
    std::size_t stage4b_candidate_evaluations = 0;
    std::size_t stage4b_enqueue_attempts = 0;
    std::size_t stage4b_successful_enqueues = 0;
    std::size_t stage4b_duplicate_suppressions = 0;
    std::size_t stage4b_empty_target_attempts = 0;
    std::size_t stage4b_empty_claim_failures = 0;
    std::size_t stage4b_commit_attempts = 0;
    std::size_t stage4b_commit_retries = 0;
    std::size_t stage4b_lock_attempts = 0;
    std::size_t stage4b_lock_contentions = 0;
    std::size_t stage4b_verification_sweeps = 0;
};
// END TEMPORARY MOVENODESFAST PERFORMANCE PROFILING
#endif

class QualityFunction {
public:
    virtual ~QualityFunction() = default;

    virtual double quality(const Graph& G,
                           const LeidenGraphStats& stats,
                           const LeidenPartition& partition) const = 0;

    virtual double deltaMove(const Graph& G,
                             const LeidenGraphStats& stats,
                             const LeidenPartition& partition,
                             Vertex v,
                             Community target_community) const = 0;

    virtual double deltaMoveFromWeights(const LeidenGraphStats& stats,
                                        const LeidenPartition& partition,
                                        Vertex v,
                                        Community target_community,
                                        double weight_to_source,
                                        double weight_to_target) const = 0;

    virtual double refinementNodeMass(const LeidenGraphStats& stats,
                                      Vertex v) const = 0;

    virtual double refinementResolution(const LeidenGraphStats& stats) const = 0;

    // True only when a positive-refinement-mass vertex cannot have a
    // nonnegative move delta to a non-neighbor community of positive mass.
    // Unknown quality functions retain the exact full-scan fallback.
    virtual bool supportsExactSparseRefinementTargets() const { return false; }

    // Stage-4B asynchronous local moving evaluates from scalar community
    // snapshots protected by the caller. Unknown quality functions fall back
    // to the serial MoveNodesFast implementation.
    virtual bool supportsConcurrentMoveEvaluation() const { return false; }
    virtual double deltaMoveFromCommunitySnapshot(
        const LeidenGraphStats& stats,
        Vertex v,
        double source_mass,
        double target_mass,
        double weight_to_source,
        double weight_to_target) const;
};

struct NeighborCommunityScratch {
    std::vector<double> weights;
    std::vector<std::size_t> marks;
    std::vector<Community> touched;
    std::size_t generation = 0;
};

// Scans G.adj[v] once and sums edge weights by neighbor community.
// Self-loops are excluded because deltaMoveFromWeights() uses the same
// convention as deltaMove(): self-loop contribution is unchanged by moving v.
// Multiple edges to the same community are accumulated.
std::unordered_map<Community, double>
BuildNeighborCommunityWeights(const Graph& G,
                              const LeidenPartition& partition,
                              Vertex v);

// Scratch-buffer overload used by MoveNodesFast(). Storage is reused across
// visits and marks make weights from earlier generations logically absent.
void BuildNeighborCommunityWeights(const Graph& G,
                                   const LeidenPartition& partition,
                                   Vertex v,
                                   NeighborCommunityScratch& scratch);

double LookupNeighborCommunityWeight(
    const NeighborCommunityScratch& scratch,
    Community community);

class CPMQualityFunction final : public QualityFunction {
public:
    explicit CPMQualityFunction(double gamma);

    double quality(const Graph& G,
                   const LeidenGraphStats& stats,
                   const LeidenPartition& partition) const override;

    double deltaMove(const Graph& G,
                     const LeidenGraphStats& stats,
                     const LeidenPartition& partition,
                     Vertex v,
                     Community target_community) const override;

    double deltaMoveFromWeights(const LeidenGraphStats& stats,
                                const LeidenPartition& partition,
                                Vertex v,
                                Community target_community,
                                double weight_to_source,
                                double weight_to_target) const override;

    double refinementNodeMass(const LeidenGraphStats& stats,
                              Vertex v) const override;

    double refinementResolution(const LeidenGraphStats& stats) const override;

    bool supportsExactSparseRefinementTargets() const override { return true; }
    bool supportsConcurrentMoveEvaluation() const override { return true; }
    double deltaMoveFromCommunitySnapshot(
        const LeidenGraphStats& stats, Vertex v, double source_mass,
        double target_mass, double weight_to_source,
        double weight_to_target) const override;

    double gamma() const { return gamma_; }

private:
    double gamma_;
};

class ModularityQualityFunction final : public QualityFunction {
public:
    explicit ModularityQualityFunction(double gamma = 1.0);

    double quality(const Graph& G,
                   const LeidenGraphStats& stats,
                   const LeidenPartition& partition) const override;

    double deltaMove(const Graph& G,
                     const LeidenGraphStats& stats,
                     const LeidenPartition& partition,
                     Vertex v,
                     Community target_community) const override;

    double deltaMoveFromWeights(const LeidenGraphStats& stats,
                                const LeidenPartition& partition,
                                Vertex v,
                                Community target_community,
                                double weight_to_source,
                                double weight_to_target) const override;

    double refinementNodeMass(const LeidenGraphStats& stats,
                              Vertex v) const override;

    double refinementResolution(const LeidenGraphStats& stats) const override;

    bool supportsExactSparseRefinementTargets() const override { return true; }
    bool supportsConcurrentMoveEvaluation() const override { return true; }
    double deltaMoveFromCommunitySnapshot(
        const LeidenGraphStats& stats, Vertex v, double source_mass,
        double target_mass, double weight_to_source,
        double weight_to_target) const override;

    double gamma() const { return gamma_; }

private:
    double gamma_;
};
