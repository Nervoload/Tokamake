#pragma once

#include <vector>

#include "tokamak/config.hpp"
#include "tokamak/particle_system.hpp"

namespace tokamak {

class IElectrostaticSourceModel {
public:
    virtual ~IElectrostaticSourceModel() = default;

    virtual void ApplySourceTerms(
        const RunConfig& runConfig,
        const TokamakConfig& tokamakConfig,
        const ParticleSystem& particles,
        double cellVolume_m3,
        std::vector<double>* chargeDensity_CPerM3) const = 0;
};

class IWallInteractionModel {
public:
    virtual ~IWallInteractionModel() = default;

    virtual ElectrostaticBoundaryCondition ResolveBoundaryCondition(
        const RunConfig& runConfig,
        const TokamakConfig& tokamakConfig) const = 0;

    virtual double DirichletBoundaryPotential_V(
        const RunConfig& runConfig,
        const TokamakConfig& tokamakConfig,
        int gridX,
        int gridY,
        int gridZ) const = 0;
};

class DefaultElectrostaticSourceModel final : public IElectrostaticSourceModel {
public:
    void ApplySourceTerms(
        const RunConfig& runConfig,
        const TokamakConfig& tokamakConfig,
        const ParticleSystem& particles,
        double cellVolume_m3,
        std::vector<double>* chargeDensity_CPerM3) const override;
};

class DefaultWallInteractionModel final : public IWallInteractionModel {
public:
    ElectrostaticBoundaryCondition ResolveBoundaryCondition(
        const RunConfig& runConfig,
        const TokamakConfig& tokamakConfig) const override;

    double DirichletBoundaryPotential_V(
        const RunConfig& runConfig,
        const TokamakConfig& tokamakConfig,
        int gridX,
        int gridY,
        int gridZ) const override;
};

}  // namespace tokamak
