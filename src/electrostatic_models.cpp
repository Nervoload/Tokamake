#include "tokamak/electrostatic_models.hpp"

namespace tokamak {

void DefaultElectrostaticSourceModel::ApplySourceTerms(
    const RunConfig& runConfig,
    const TokamakConfig& tokamakConfig,
    const ParticleSystem& particles,
    double cellVolume_m3,
    std::vector<double>* chargeDensity_CPerM3) const {
    (void)runConfig;
    (void)tokamakConfig;
    (void)particles;
    (void)cellVolume_m3;
    (void)chargeDensity_CPerM3;
}

ElectrostaticBoundaryCondition DefaultWallInteractionModel::ResolveBoundaryCondition(
    const RunConfig& runConfig,
    const TokamakConfig& tokamakConfig) const {
    (void)tokamakConfig;
    return runConfig.electrostaticBoundaryCondition;
}

double DefaultWallInteractionModel::DirichletBoundaryPotential_V(
    const RunConfig& runConfig,
    const TokamakConfig& tokamakConfig,
    int gridX,
    int gridY,
    int gridZ) const {
    (void)runConfig;
    (void)tokamakConfig;
    (void)gridX;
    (void)gridY;
    (void)gridZ;
    return 0.0;
}

}  // namespace tokamak
