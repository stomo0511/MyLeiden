#pragma once

#include "LeidenTypes.hpp"

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
};

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

    double gamma() const { return gamma_; }

private:
    double gamma_;
};
