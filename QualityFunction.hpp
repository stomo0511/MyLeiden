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

    double gamma() const { return gamma_; }

private:
    double gamma_;
};
