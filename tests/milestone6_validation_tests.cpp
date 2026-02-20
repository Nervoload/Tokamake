#include <cmath>
#include <cstdint>
#include <vector>

#include "gtest/gtest.h"
#include "tokamak/collision.hpp"
#include "tokamak/electrostatic_field.hpp"
#include "tokamak/particle_push.hpp"
#include "tokamak/spatial_grid.hpp"
#include "validation_harness.hpp"

namespace {

using tokamak_validation::ExpectNearWithTolerance;
using tokamak_validation::ScalarTolerance;

constexpr double kPi = 3.14159265358979323846;
constexpr double kEpsilon0 = 8.8541878128e-12;

// Standalone finite-difference harness for Poisson math sanity checks.
// This does not exercise the engine's runtime electrostatic field path.
struct Poisson1DSolveResult {
    std::vector<double> potential_V;
    double residualL2 = 0.0;
};

double ChargeDensitySinusoid(double x_m, double length_m, double rho0_CPerM3) {
    return rho0_CPerM3 * std::sin((kPi * x_m) / length_m);
}

double AnalyticPotentialSinusoid(double x_m, double length_m, double rho0_CPerM3) {
    const double amplitude = (rho0_CPerM3 / kEpsilon0) * ((length_m * length_m) / (kPi * kPi));
    return amplitude * std::sin((kPi * x_m) / length_m);
}

Poisson1DSolveResult SolvePoisson1DDirichletTridiagonal(
    int nodeCount,
    double length_m,
    double rho0_CPerM3) {
    EXPECT_TRUE(nodeCount >= 5);
    EXPECT_TRUE(length_m > 0.0);

    const int interiorCount = nodeCount - 2;
    const double dx_m = length_m / static_cast<double>(nodeCount - 1);
    const double dxSq_m2 = dx_m * dx_m;

    std::vector<double> rhs(nodeCount, 0.0);
    for (int i = 1; i < nodeCount - 1; ++i) {
        const double x_m = dx_m * static_cast<double>(i);
        rhs[static_cast<std::size_t>(i)] = -ChargeDensitySinusoid(x_m, length_m, rho0_CPerM3) / kEpsilon0;
    }

    std::vector<double> a(static_cast<std::size_t>(interiorCount), 1.0);
    std::vector<double> b(static_cast<std::size_t>(interiorCount), -2.0);
    std::vector<double> c(static_cast<std::size_t>(interiorCount), 1.0);
    std::vector<double> d(static_cast<std::size_t>(interiorCount), 0.0);

    for (int i = 0; i < interiorCount; ++i) {
        const int gridIndex = i + 1;
        d[static_cast<std::size_t>(i)] = dxSq_m2 * rhs[static_cast<std::size_t>(gridIndex)];
    }

    for (int i = 1; i < interiorCount; ++i) {
        const double m = a[static_cast<std::size_t>(i)] / b[static_cast<std::size_t>(i - 1)];
        b[static_cast<std::size_t>(i)] -= m * c[static_cast<std::size_t>(i - 1)];
        d[static_cast<std::size_t>(i)] -= m * d[static_cast<std::size_t>(i - 1)];
    }

    std::vector<double> interiorSolution(static_cast<std::size_t>(interiorCount), 0.0);
    interiorSolution[static_cast<std::size_t>(interiorCount - 1)] =
        d[static_cast<std::size_t>(interiorCount - 1)] / b[static_cast<std::size_t>(interiorCount - 1)];

    for (int i = interiorCount - 2; i >= 0; --i) {
        interiorSolution[static_cast<std::size_t>(i)] =
            (d[static_cast<std::size_t>(i)] -
             (c[static_cast<std::size_t>(i)] * interiorSolution[static_cast<std::size_t>(i + 1)])) /
            b[static_cast<std::size_t>(i)];
    }

    std::vector<double> potential(static_cast<std::size_t>(nodeCount), 0.0);
    for (int i = 0; i < interiorCount; ++i) {
        potential[static_cast<std::size_t>(i + 1)] = interiorSolution[static_cast<std::size_t>(i)];
    }

    double residualSqSum = 0.0;
    for (int i = 1; i < nodeCount - 1; ++i) {
        const double laplacian =
            (potential[static_cast<std::size_t>(i - 1)] -
             (2.0 * potential[static_cast<std::size_t>(i)]) +
             potential[static_cast<std::size_t>(i + 1)]) /
            dxSq_m2;
        const double residual = laplacian - rhs[static_cast<std::size_t>(i)];
        residualSqSum += residual * residual;
    }

    Poisson1DSolveResult out;
    out.residualL2 = std::sqrt(residualSqSum / static_cast<double>(interiorCount));
    out.potential_V = std::move(potential);
    return out;
}

tokamak::SpeciesCounts CountSpecies(const tokamak::ParticleSystem& particles) {
    tokamak::SpeciesCounts counts;
    const auto& species = particles.Species();
    for (const auto type : species) {
        if (type == tokamak::ParticleType::Deuterium) {
            ++counts.deuterium;
        } else if (type == tokamak::ParticleType::Tritium) {
            ++counts.tritium;
        } else if (type == tokamak::ParticleType::Helium) {
            ++counts.helium;
        }
    }
    return counts;
}

int GridIndexFromXYZ(int x, int y, int z, int gridWidth) {
    return x + (gridWidth * (y + (gridWidth * z)));
}

std::vector<uint32_t> AggregateFineCountsToCoarse(const tokamak::SpatialGrid& coarse, const tokamak::SpatialGrid& fine) {
    const int ratio = fine.gridWidth / coarse.gridWidth;
    EXPECT_GT(ratio, 0);
    EXPECT_EQ(fine.gridWidth % coarse.gridWidth, 0);

    std::vector<uint32_t> aggregated(coarse.cellCounts.size(), 0U);
    for (int cz = 0; cz < coarse.gridWidth; ++cz) {
        for (int cy = 0; cy < coarse.gridWidth; ++cy) {
            for (int cx = 0; cx < coarse.gridWidth; ++cx) {
                uint32_t sum = 0U;
                for (int fz = cz * ratio; fz < (cz + 1) * ratio; ++fz) {
                    for (int fy = cy * ratio; fy < (cy + 1) * ratio; ++fy) {
                        for (int fx = cx * ratio; fx < (cx + 1) * ratio; ++fx) {
                            const int fineIndex = GridIndexFromXYZ(fx, fy, fz, fine.gridWidth);
                            sum += fine.cellCounts[static_cast<std::size_t>(fineIndex)];
                        }
                    }
                }

                const int coarseIndex = GridIndexFromXYZ(cx, cy, cz, coarse.gridWidth);
                aggregated[static_cast<std::size_t>(coarseIndex)] = sum;
            }
        }
    }
    return aggregated;
}

}  // namespace

TEST(Milestone6ValidationTest, UniformBGyroFrequencyAndRadiusMatchAnalytic) {
    const float bMag_T = 1.8f;
    const float vPerp_mPerS = 2.5e5f;
    const float qOverM = tokamak::constants::kElementaryCharge_C / tokamak::constants::kMassDeuterium_kg;
    const double omegaExpected_radPerS = std::abs(static_cast<double>(qOverM) * static_cast<double>(bMag_T));
    const double periodExpected_s = (2.0 * kPi) / omegaExpected_radPerS;

    const float dt_s = static_cast<float>(periodExpected_s / 240.0);
    const int totalSteps = 240 * 80;

    tokamak::Vec3 position(0.0f, 0.0f, 0.0f);
    tokamak::Vec3 velocity(vPerp_mPerS, 0.0f, 0.0f);
    const tokamak::Vec3 eField(0.0f, 0.0f, 0.0f);
    const tokamak::Vec3 bField(0.0f, 0.0f, bMag_T);
    const tokamak::Vec3 bHat(0.0f, 0.0f, 1.0f);

    double accumulatedPhase_rad = 0.0;
    double previousPhase_rad = std::atan2(static_cast<double>(velocity.y), static_cast<double>(velocity.x));
    double radiusAccumulator_m = 0.0;
    int radiusSamples = 0;

    for (int i = 0; i < totalSteps; ++i) {
        velocity = tokamak::BorisVelocityStep(velocity, eField, bField, qOverM, dt_s);
        position += velocity * dt_s;

        const double phase_rad = std::atan2(static_cast<double>(velocity.y), static_cast<double>(velocity.x));
        accumulatedPhase_rad += tokamak_validation::WrappedAngleDelta(previousPhase_rad, phase_rad);
        previousPhase_rad = phase_rad;

        if (i > (totalSteps / 4)) {
            const float gyroSigned_radPerS = qOverM * bMag_T;
            const tokamak::Vec3 guidingCenter = position + (tokamak::Vec3::Cross(velocity, bHat) / gyroSigned_radPerS);
            radiusAccumulator_m += static_cast<double>((position - guidingCenter).Magnitude());
            ++radiusSamples;
        }
    }

    ASSERT_TRUE(radiusSamples > 0);

    const double omegaObserved_radPerS =
        std::abs(accumulatedPhase_rad) / (static_cast<double>(totalSteps) * static_cast<double>(dt_s));
    const double gyroradiusObserved_m = radiusAccumulator_m / static_cast<double>(radiusSamples);
    const double gyroradiusExpected_m = static_cast<double>(vPerp_mPerS) / omegaExpected_radPerS;

    const ScalarTolerance gyroFrequencyTolerance{
        0.0,
        0.01,
        "1% relative tolerance covers discretization from finite dt while still catching integrator regressions.",
    };
    const ScalarTolerance gyroRadiusTolerance{
        0.0,
        0.015,
        "1.5% relative tolerance allows orbit-center estimation noise but still enforces correct Larmor scaling.",
    };

    ExpectNearWithTolerance("gyrofrequency(rad/s)", omegaObserved_radPerS, omegaExpected_radPerS, gyroFrequencyTolerance);
    ExpectNearWithTolerance("gyroradius(m)", gyroradiusObserved_m, gyroradiusExpected_m, gyroRadiusTolerance);
}

TEST(Milestone6ValidationTest, ExBDriftMatchesAnalyticVelocity) {
    const tokamak::Vec3 eField(1200.0f, 0.0f, 0.0f);
    const tokamak::Vec3 bField(0.0f, 0.0f, 2.0f);
    const float qOverM = tokamak::constants::kElementaryCharge_C / tokamak::constants::kMassDeuterium_kg;

    const float omega_radPerS = std::abs(qOverM * bField.z);
    const float period_s = static_cast<float>((2.0 * kPi) / static_cast<double>(omega_radPerS));
    const float dt_s = period_s / 240.0f;
    const int totalSteps = 240 * 120;

    tokamak::Vec3 position(0.0f, 0.0f, 0.0f);
    tokamak::Vec3 velocity(0.0f, 0.0f, 0.0f);

    tokamak::Vec3 velocityAccumulator(0.0f, 0.0f, 0.0f);
    int velocitySamples = 0;

    for (int i = 0; i < totalSteps; ++i) {
        velocity = tokamak::BorisVelocityStep(velocity, eField, bField, qOverM, dt_s);
        position += velocity * dt_s;

        if (i > (totalSteps / 3)) {
            velocityAccumulator += velocity;
            ++velocitySamples;
        }
    }

    ASSERT_TRUE(velocitySamples > 0);
    (void)position;

    const tokamak::Vec3 meanVelocity = velocityAccumulator / static_cast<float>(velocitySamples);
    const tokamak::Vec3 expectedDrift = tokamak::Vec3::Cross(eField, bField) / tokamak::Vec3::Dot(bField, bField);

    const ScalarTolerance driftTolerance{
        5.0,
        0.02,
        "2% relative + 5 m/s absolute tolerance covers finite-step gyro averaging while enforcing correct ExB drift.",
    };

    ExpectNearWithTolerance("ExB drift vx(m/s)", meanVelocity.x, expectedDrift.x, driftTolerance);
    ExpectNearWithTolerance("ExB drift vy(m/s)", meanVelocity.y, expectedDrift.y, driftTolerance);
    ExpectNearWithTolerance("ExB drift vz(m/s)", meanVelocity.z, expectedDrift.z, driftTolerance);
}

TEST(Milestone6ValidationTest, WallReflectionConservesSpeed) {
    const float majorRadius_m = 2.0f;
    const float minorRadius_m = 0.5f;
    tokamak::Vec3 position(2.58f, 0.0f, 0.0f);
    tokamak::Vec3 velocity(1.25e5f, -2.0e4f, 3.5e4f);

    const float speedBefore_mPerS = velocity.Magnitude();
    const bool reflected = tokamak::ReflectAtTokamakWall(&position, &velocity, majorRadius_m, minorRadius_m);

    ASSERT_TRUE(reflected);
    const float speedAfter_mPerS = velocity.Magnitude();

    const ScalarTolerance speedTolerance{
        1.0e-3,
        1.0e-6,
        "Reflection is algebraically elastic; tolerance only guards floating-point roundoff.",
    };
    ExpectNearWithTolerance("reflective wall speed conservation", speedAfter_mPerS, speedBefore_mPerS, speedTolerance);

    const float radialPos_m = std::sqrt((position.x * position.x) + (position.y * position.y));
    const tokamak::Vec3 coreCenter(majorRadius_m * (position.x / radialPos_m), majorRadius_m * (position.y / radialPos_m), 0.0f);
    const tokamak::Vec3 normal = (position - coreCenter).Normalized();
    const float radialTube_m = std::sqrt(std::pow(radialPos_m - majorRadius_m, 2.0f) + (position.z * position.z));

    EXPECT_LT(tokamak::Vec3::Dot(velocity, normal), 0.0f);
    EXPECT_NEAR(radialTube_m, minorRadius_m * 0.99f, 1.0e-5f);
}

TEST(Milestone6ValidationTest, CollisionSelectionAndApplicationAreSeedReproducible) {
    constexpr std::size_t kPairsPerSpecies = 800;
    constexpr float kDt_s = 2.0e-12f;
    constexpr uint32_t kSeed = 2026U;

    tokamak::ParticleSystem particlesA = tokamak_validation::MakeHotDtCollisionPool(kPairsPerSpecies);
    tokamak::ParticleSystem particlesB = tokamak_validation::MakeHotDtCollisionPool(kPairsPerSpecies);

    tokamak::SpatialGrid gridA(2.5f, 0.2f);
    tokamak::SpatialGrid gridB(2.5f, 0.2f);
    tokamak::SortParticlesIntoGrid(particlesA, gridA, 2.5f);
    tokamak::SortParticlesIntoGrid(particlesB, gridB, 2.5f);

    EXPECT_TRUE(tokamak_validation::SortedIdsFormCompletePermutation(gridA.sortedParticleIDs, particlesA.Size()));
    EXPECT_TRUE(tokamak_validation::SortedIdsFormCompletePermutation(gridB.sortedParticleIDs, particlesB.Size()));

    std::mt19937 rngA(kSeed);
    std::mt19937 rngB(kSeed);

    std::vector<tokamak::PendingFusionEvent> eventsA;
    std::vector<tokamak::PendingFusionEvent> eventsB;
    const auto summaryA = tokamak::SelectCollisionEvents(particlesA, gridA, rngA, kDt_s, eventsA);
    const auto summaryB = tokamak::SelectCollisionEvents(particlesB, gridB, rngB, kDt_s, eventsB);

    EXPECT_EQ(summaryA.fusionAttempts, summaryB.fusionAttempts);
    EXPECT_EQ(summaryA.selectedEvents, summaryB.selectedEvents);
    EXPECT_EQ(summaryA.maxReactionsInCell, summaryB.maxReactionsInCell);
    ASSERT_EQ(eventsA.size(), eventsB.size());

    for (std::size_t i = 0; i < eventsA.size(); ++i) {
        EXPECT_EQ(eventsA[i].p1, eventsB[i].p1);
        EXPECT_EQ(eventsA[i].p2, eventsB[i].p2);
        EXPECT_NEAR(eventsA[i].centerOfMassPos.x, eventsB[i].centerOfMassPos.x, 0.0f);
        EXPECT_NEAR(eventsA[i].centerOfMassPos.y, eventsB[i].centerOfMassPos.y, 0.0f);
        EXPECT_NEAR(eventsA[i].centerOfMassPos.z, eventsB[i].centerOfMassPos.z, 0.0f);
        EXPECT_NEAR(eventsA[i].heliumVelocity.x, eventsB[i].heliumVelocity.x, 0.0f);
        EXPECT_NEAR(eventsA[i].heliumVelocity.y, eventsB[i].heliumVelocity.y, 0.0f);
        EXPECT_NEAR(eventsA[i].heliumVelocity.z, eventsB[i].heliumVelocity.z, 0.0f);
    }

    tokamak::RuntimeCounters countersA;
    tokamak::RuntimeCounters countersB;
    tokamak::EnergyChargeBudget budgetA;
    tokamak::EnergyChargeBudget budgetB;
    tokamak::ApplyCollisionEvents(particlesA, eventsA, countersA, budgetA);
    tokamak::ApplyCollisionEvents(particlesB, eventsB, countersB, budgetB);
    particlesA.Compact();
    particlesB.Compact();

    const auto speciesA = CountSpecies(particlesA);
    const auto speciesB = CountSpecies(particlesB);

    EXPECT_EQ(countersA.fusionAccepted, countersB.fusionAccepted);
    EXPECT_EQ(countersA.rejectedFusionAsh, countersB.rejectedFusionAsh);
    EXPECT_EQ(particlesA.fusionCountTotal, particlesB.fusionCountTotal);
    EXPECT_EQ(speciesA.deuterium, speciesB.deuterium);
    EXPECT_EQ(speciesA.tritium, speciesB.tritium);
    EXPECT_EQ(speciesA.helium, speciesB.helium);
    EXPECT_NEAR(budgetA.fusionAlphaInjected_J, budgetB.fusionAlphaInjected_J, 0.0);
}

TEST(Milestone6ValidationTest, GridConsistencyConvergesAcrossResolution) {
    constexpr std::size_t kParticleCount = 4000;
    const auto positions = tokamak_validation::MakeDeterministicPositions(kParticleCount, -2.49f, 2.49f, 1337U);

    tokamak::ParticleSystem particles(kParticleCount + 8U);
    for (const auto& position : positions) {
        ASSERT_TRUE(particles.AddParticle(
            position,
            tokamak::Vec3(0.0f, 0.0f, 0.0f),
            tokamak::constants::kMassDeuterium_kg,
            tokamak::constants::kElementaryCharge_C,
            tokamak::ParticleType::Deuterium));
    }

    tokamak::SpatialGrid coarseGrid(2.5f, 0.5f);
    tokamak::SpatialGrid fineGrid(2.5f, 0.25f);
    tokamak::SortParticlesIntoGrid(particles, coarseGrid, 2.5f);
    tokamak::SortParticlesIntoGrid(particles, fineGrid, 2.5f);

    EXPECT_TRUE(tokamak_validation::SortedIdsFormCompletePermutation(coarseGrid.sortedParticleIDs, particles.Size()));
    EXPECT_TRUE(tokamak_validation::SortedIdsFormCompletePermutation(fineGrid.sortedParticleIDs, particles.Size()));

    const auto aggregatedFineCounts = AggregateFineCountsToCoarse(coarseGrid, fineGrid);
    ASSERT_EQ(aggregatedFineCounts.size(), coarseGrid.cellCounts.size());

    uint64_t mismatchCells = 0;
    for (std::size_t i = 0; i < coarseGrid.cellCounts.size(); ++i) {
        if (aggregatedFineCounts[i] != coarseGrid.cellCounts[i]) {
            ++mismatchCells;
        }
    }

    EXPECT_EQ(mismatchCells, static_cast<uint64_t>(0));
}

TEST(Milestone6ValidationTest, SpatialGridGetCellIndexReportsClamp) {
    tokamak::SpatialGrid grid(2.5f, 0.2f);
    bool wasClamped = false;

    const int clampedIndex = grid.GetCellIndex(tokamak::Vec3(1000.0f, -1000.0f, 1000.0f), 2.5f, &wasClamped);
    EXPECT_TRUE(grid.IsValidCellIndex(clampedIndex));
    EXPECT_TRUE(wasClamped);

    wasClamped = false;
    const int insideIndex = grid.GetCellIndex(tokamak::Vec3(0.0f, 0.0f, 0.0f), 2.5f, &wasClamped);
    EXPECT_TRUE(grid.IsValidCellIndex(insideIndex));
    EXPECT_FALSE(wasClamped);
}

TEST(Milestone6ValidationTest, ElectrostaticPoissonStandaloneNumericHarnessMatchesAnalytic) {
    constexpr int kNodeCount = 129;
    constexpr double kLength_m = 1.0;
    constexpr double kRho0_CPerM3 = 2.0e-6;
    constexpr double kResidualTolerance = 1.0e-6;

    const Poisson1DSolveResult solution = SolvePoisson1DDirichletTridiagonal(
        kNodeCount,
        kLength_m,
        kRho0_CPerM3);

    EXPECT_TRUE(solution.residualL2 <= kResidualTolerance);
    ASSERT_EQ(solution.potential_V.size(), static_cast<std::size_t>(kNodeCount));

    const double dx_m = kLength_m / static_cast<double>(kNodeCount - 1);

    double potentialErrorSqSum = 0.0;
    double analyticPotentialSqSum = 0.0;
    for (int i = 1; i < kNodeCount - 1; ++i) {
        const double x_m = static_cast<double>(i) * dx_m;
        const double analytic = AnalyticPotentialSinusoid(x_m, kLength_m, kRho0_CPerM3);
        const double error = solution.potential_V[static_cast<std::size_t>(i)] - analytic;
        potentialErrorSqSum += error * error;
        analyticPotentialSqSum += analytic * analytic;
    }
    const double potentialRelativeL2 = std::sqrt(potentialErrorSqSum / analyticPotentialSqSum);

    const ScalarTolerance potentialRelativeTolerance{
        3.0e-3,
        0.0,
        "Second-order finite differences at 129 nodes should keep potential relative L2 error below 0.3%.",
    };
    ExpectNearWithTolerance(
        "poisson potential relative L2",
        potentialRelativeL2,
        0.0,
        potentialRelativeTolerance);

    double electricErrorSqSum = 0.0;
    double analyticElectricSqSum = 0.0;
    const double amplitude = (kRho0_CPerM3 / kEpsilon0) * ((kLength_m * kLength_m) / (kPi * kPi));
    for (int i = 1; i < kNodeCount - 1; ++i) {
        const double x_m = static_cast<double>(i) * dx_m;
        const double electricNumeric_VPerM =
            -(solution.potential_V[static_cast<std::size_t>(i + 1)] -
              solution.potential_V[static_cast<std::size_t>(i - 1)]) /
            (2.0 * dx_m);
        const double electricAnalytic_VPerM = -amplitude * (kPi / kLength_m) * std::cos((kPi * x_m) / kLength_m);
        const double error = electricNumeric_VPerM - electricAnalytic_VPerM;
        electricErrorSqSum += error * error;
        analyticElectricSqSum += electricAnalytic_VPerM * electricAnalytic_VPerM;
    }
    const double electricRelativeL2 = std::sqrt(electricErrorSqSum / analyticElectricSqSum);

    const ScalarTolerance electricRelativeTolerance{
        0.015,
        0.0,
        "Field reconstruction from central differences is noisier than potential; 1.5% relative L2 is expected.",
    };
    ExpectNearWithTolerance(
        "poisson electric-field relative L2",
        electricRelativeL2,
        0.0,
        electricRelativeTolerance);
}
