#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <random>
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

std::size_t CellIndex(int x, int y, int z, int width) {
    return static_cast<std::size_t>(x + (width * (y + (width * z))));
}

double FieldRoughnessL2(
    const tokamak::ElectrostaticMeshGeometry& geometry,
    const std::vector<tokamak::Vec3>& field) {
    if (geometry.gridWidth < 3 || field.empty()) {
        return 0.0;
    }

    const int width = geometry.gridWidth;
    double sum = 0.0;
    std::size_t samples = 0;
    for (int z = 1; z < width - 1; ++z) {
        for (int y = 1; y < width - 1; ++y) {
            for (int x = 1; x < width - 1; ++x) {
                const tokamak::Vec3 center = field[CellIndex(x, y, z, width)];
                const tokamak::Vec3 neighborAvg =
                    (field[CellIndex(x + 1, y, z, width)] + field[CellIndex(x - 1, y, z, width)] +
                     field[CellIndex(x, y + 1, z, width)] + field[CellIndex(x, y - 1, z, width)] +
                     field[CellIndex(x, y, z + 1, width)] + field[CellIndex(x, y, z - 1, width)]) /
                    6.0f;
                const tokamak::Vec3 delta = center - neighborAvg;
                sum += static_cast<double>(tokamak::Vec3::Dot(delta, delta));
                ++samples;
            }
        }
    }

    if (samples == 0) {
        return 0.0;
    }
    return std::sqrt(sum / static_cast<double>(samples));
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

TEST(Milestone3ElectrostaticFieldTest, CicDepositConservesTotalCharge) {
    const tokamak::ElectrostaticMeshGeometry geometry = TestGeometry();
    std::vector<double> rho;

    const std::vector<tokamak::Vec3> positions = {
        tokamak::Vec3(-0.7f, 0.2f, 0.5f),
        tokamak::Vec3(0.9f, -0.4f, -0.3f),
        tokamak::Vec3(0.1f, 1.0f, -0.8f),
        tokamak::Vec3(-1.1f, -0.6f, 0.7f),
    };
    const std::vector<float> charges = {
        tokamak::constants::kElementaryCharge_C,
        tokamak::constants::kElementaryCharge_C,
        tokamak::constants::kElementaryCharge_C,
        -tokamak::constants::kElementaryCharge_C,
    };
    constexpr float kMacroWeight = 4.2e9f;

    tokamak::DepositChargeDensity(
        positions,
        charges,
        kMacroWeight,
        geometry,
        tokamak::ChargeAssignmentScheme::CIC,
        &rho);

    const double cellVolume_m3 = geometry.cellSize_m * geometry.cellSize_m * geometry.cellSize_m;
    const double depositedCharge_C = std::accumulate(rho.begin(), rho.end(), 0.0) * cellVolume_m3;

    double expectedCharge_C = 0.0;
    for (const float charge : charges) {
        expectedCharge_C += static_cast<double>(charge) * static_cast<double>(kMacroWeight);
    }

    EXPECT_NEAR(depositedCharge_C, expectedCharge_C, std::max(1.0e-18, std::abs(expectedCharge_C) * 1.0e-5));
}

TEST(Milestone3ElectrostaticFieldTest, CicInterpolationReproducesLinearField) {
    const tokamak::ElectrostaticMeshGeometry geometry = TestGeometry();
    const int width = geometry.gridWidth;

    std::vector<tokamak::Vec3> linearField(tokamak::ElectrostaticCellCount(geometry), tokamak::Vec3(0.0f, 0.0f, 0.0f));
    for (int z = 0; z < width; ++z) {
        for (int y = 0; y < width; ++y) {
            for (int x = 0; x < width; ++x) {
                linearField[CellIndex(x, y, z, width)] = tokamak::Vec3(
                    1.0f + (0.2f * static_cast<float>(x)) + (0.3f * static_cast<float>(y)) + (0.4f * static_cast<float>(z)),
                    -2.0f + (0.6f * static_cast<float>(x)) - (0.1f * static_cast<float>(y)) + (0.25f * static_cast<float>(z)),
                    0.5f - (0.15f * static_cast<float>(x)) + (0.45f * static_cast<float>(y)) + (0.05f * static_cast<float>(z)));
            }
        }
    }

    const std::vector<tokamak::Vec3> samplePositions = {
        tokamak::Vec3(-1.3125f, -0.9375f, -0.5625f),
        tokamak::Vec3(-0.8125f, 0.1875f, 0.9375f),
        tokamak::Vec3(0.4375f, -0.4375f, 0.3125f),
        tokamak::Vec3(1.0625f, 0.8125f, -0.1875f),
    };

    for (const tokamak::Vec3& pos : samplePositions) {
        const double gx = (static_cast<double>(pos.x) + geometry.offset_m) / geometry.cellSize_m;
        const double gy = (static_cast<double>(pos.y) + geometry.offset_m) / geometry.cellSize_m;
        const double gz = (static_cast<double>(pos.z) + geometry.offset_m) / geometry.cellSize_m;

        const tokamak::Vec3 expected(
            static_cast<float>(1.0 + (0.2 * gx) + (0.3 * gy) + (0.4 * gz)),
            static_cast<float>(-2.0 + (0.6 * gx) - (0.1 * gy) + (0.25 * gz)),
            static_cast<float>(0.5 - (0.15 * gx) + (0.45 * gy) + (0.05 * gz)));

        const tokamak::Vec3 sampled = tokamak::SampleElectricField(
            geometry,
            linearField,
            pos,
            tokamak::ChargeAssignmentScheme::CIC);

        EXPECT_NEAR(sampled.x, expected.x, 1.0e-4f);
        EXPECT_NEAR(sampled.y, expected.y, 1.0e-4f);
        EXPECT_NEAR(sampled.z, expected.z, 1.0e-4f);
    }
}

TEST(Milestone3ElectrostaticFieldTest, CicProducesSmootherFieldThanNgpForSameChargeCloud) {
    tokamak::ElectrostaticMeshGeometry geometry;
    geometry.gridWidth = 16;
    geometry.offset_m = 2.5;
    geometry.cellSize_m = (2.0 * geometry.offset_m) / static_cast<double>(geometry.gridWidth);

    std::mt19937 rng(4242);
    std::uniform_real_distribution<float> posDist(-2.2f, 2.2f);

    std::vector<tokamak::Vec3> positions;
    std::vector<float> charges;
    positions.reserve(600);
    charges.reserve(600);
    for (int i = 0; i < 600; ++i) {
        positions.push_back(tokamak::Vec3(posDist(rng), posDist(rng), posDist(rng)));
        charges.push_back(tokamak::constants::kElementaryCharge_C);
    }

    tokamak::ElectrostaticSolveConfig solveConfig;
    solveConfig.boundaryCondition = tokamak::ElectrostaticBoundaryCondition::DirichletZero;
    solveConfig.tolerance = 1.0e-5;
    solveConfig.maxIterations = 400;
    solveConfig.sorOmega = 1.6;

    auto roughnessForScheme = [&](tokamak::ChargeAssignmentScheme scheme) {
        std::vector<double> rho;
        tokamak::DepositChargeDensity(positions, charges, 1.0e9f, geometry, scheme, &rho);
        tokamak::ApplyNeutralizingBackground(1.0, &rho);

        std::vector<double> potential;
        const tokamak::ElectrostaticSolveResult solve =
            tokamak::SolvePoissonSor(geometry, solveConfig, 0.0, rho, &potential);
        EXPECT_TRUE(solve.residual.residualAvailable);

        std::vector<tokamak::Vec3> eField;
        tokamak::ReconstructElectricField(geometry, potential, &eField);
        return FieldRoughnessL2(geometry, eField);
    };

    const double ngpRoughness = roughnessForScheme(tokamak::ChargeAssignmentScheme::NGP);
    const double cicRoughness = roughnessForScheme(tokamak::ChargeAssignmentScheme::CIC);

    EXPECT_GT(ngpRoughness, 0.0);
    EXPECT_GT(cicRoughness, 0.0);
    EXPECT_LT(cicRoughness, ngpRoughness);
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
