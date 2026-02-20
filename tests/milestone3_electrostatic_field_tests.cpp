#include <cmath>
#include <vector>

#include "gtest/gtest.h"
#include "tokamak/electrostatic_field.hpp"

namespace {

tokamak::ElectrostaticMeshGeometry TestGeometry() {
    tokamak::ElectrostaticMeshGeometry geometry;
    geometry.gridWidth = 8;
    geometry.offset_m = 2.5;
    geometry.cellSize_m = (2.0 * geometry.offset_m) / static_cast<double>(geometry.gridWidth);
    return geometry;
}

}  // namespace

TEST(Milestone3ElectrostaticFieldTest, ElectrostaticCellCountMatchesGridVolume) {
    const tokamak::ElectrostaticMeshGeometry geometry = TestGeometry();
    EXPECT_EQ(tokamak::ElectrostaticCellCount(geometry), static_cast<std::size_t>(512));
}

TEST(Milestone3ElectrostaticFieldTest, DepositChargeDensityProducesFiniteValues) {
    const tokamak::ElectrostaticMeshGeometry geometry = TestGeometry();
    std::vector<double> rho;

    const std::vector<tokamak::Vec3> positions = {
        tokamak::Vec3(0.1f, -0.2f, 0.3f),
        tokamak::Vec3(-0.3f, 0.2f, -0.1f),
    };
    const std::vector<float> charges = {
        tokamak::constants::kElementaryCharge_C,
        tokamak::constants::kElementaryCharge_C,
    };

    tokamak::DepositChargeDensity(
        positions,
        charges,
        1.0e10f,
        geometry,
        tokamak::ChargeAssignmentScheme::CIC,
        &rho);

    EXPECT_EQ(rho.size(), tokamak::ElectrostaticCellCount(geometry));
    for (double value : rho) {
        EXPECT_TRUE(std::isfinite(value));
    }
}

TEST(Milestone3ElectrostaticFieldTest, SolvePoissonSorWithZeroChargeConverges) {
    const tokamak::ElectrostaticMeshGeometry geometry = TestGeometry();
    tokamak::ElectrostaticSolveConfig solveConfig;
    solveConfig.boundaryCondition = tokamak::ElectrostaticBoundaryCondition::DirichletZero;
    solveConfig.tolerance = 1.0e-6;
    solveConfig.maxIterations = 200;
    solveConfig.sorOmega = 1.5;

    std::vector<double> rho(tokamak::ElectrostaticCellCount(geometry), 0.0);
    std::vector<double> phi;

    const tokamak::ElectrostaticSolveResult result =
        tokamak::SolvePoissonSor(geometry, solveConfig, 0.0, rho, &phi);

    EXPECT_EQ(phi.size(), tokamak::ElectrostaticCellCount(geometry));
    EXPECT_TRUE(result.residual.residualAvailable);
    EXPECT_TRUE(std::isfinite(result.residual.residualL2));
    EXPECT_TRUE(result.residual.converged);
}
