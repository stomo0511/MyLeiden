#pragma once

#include "LeidenTypes.hpp"

#include <unordered_map>

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

// Scans G.adj[v] once and sums edge weights by neighbor community.
// Self-loops are excluded because deltaMoveFromWeights() uses the same
// convention as deltaMove(): self-loop contribution is unchanged by moving v.
// Multiple edges to the same community are accumulated.
std::unordered_map<Community, double>
BuildNeighborCommunityWeights(const Graph& G,
                              const LeidenPartition& partition,
                              Vertex v);

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
