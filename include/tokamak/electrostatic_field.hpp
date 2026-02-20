#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "tokamak/config.hpp"
#include "tokamak/diagnostics.hpp"
#include "tokamak/types.hpp"

namespace tokamak {

struct ElectrostaticMeshGeometry {
    int gridWidth = 0;
    double cellSize_m = 0.0;
    double offset_m = 0.0;
};

struct ElectrostaticSolveConfig {
    ElectrostaticBoundaryCondition boundaryCondition = ElectrostaticBoundaryCondition::DirichletZero;
    ChargeAssignmentScheme chargeAssignmentScheme = ChargeAssignmentScheme::CIC;
    double tolerance = 1.0e-6;
    uint32_t maxIterations = 5000;
    double sorOmega = 1.7;
    double neutralizingBackgroundFraction = 1.0;
};

struct ElectrostaticSolveResult {
    SolverResidualSnapshot residual;
    std::vector<double> residualHistoryL2;
};

std::size_t ElectrostaticCellCount(const ElectrostaticMeshGeometry& geometry);

void DepositChargeDensity(
    const std::vector<Vec3>& positions,
    const std::vector<float>& charges_C,
    float macroWeight,
    const ElectrostaticMeshGeometry& geometry,
    ChargeAssignmentScheme assignmentScheme,
    std::vector<double>* chargeDensity_CPerM3);

void ApplyNeutralizingBackground(
    double neutralizingBackgroundFraction,
    std::vector<double>* chargeDensity_CPerM3);

ElectrostaticSolveResult SolvePoissonSor(
    const ElectrostaticMeshGeometry& geometry,
    const ElectrostaticSolveConfig& solveConfig,
    double dirichletBoundaryPotential_V,
    const std::vector<double>& chargeDensity_CPerM3,
    std::vector<double>* potential_V);

void ReconstructElectricField(
    const ElectrostaticMeshGeometry& geometry,
    const std::vector<double>& potential_V,
    std::vector<Vec3>* electricField_VPerM);

Vec3 SampleElectricField(
    const ElectrostaticMeshGeometry& geometry,
    const std::vector<Vec3>& electricField_VPerM,
    const Vec3& position,
    ChargeAssignmentScheme assignmentScheme);

}  // namespace tokamak
