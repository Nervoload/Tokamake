#include <algorithm>
#include <random>
#include <vector>

#include "gtest/gtest.h"

#include "tokamak/collision.hpp"
#include "tokamak/reactivity.hpp"

namespace {

double SelectSinglePairProbability(
    float relativeSpeed_mPerS,
    float deuteriumWeight,
    float tritiumWeight,
    float dt_s,
    const tokamak::RunConfig& runConfig) {
    tokamak::ParticleSystem particles(4);
    const tokamak::Vec3 pos(2.0f, 0.0f, 0.0f);
    const tokamak::Vec3 velD(0.5f * relativeSpeed_mPerS, 0.0f, 0.0f);
    const tokamak::Vec3 velT(-0.5f * relativeSpeed_mPerS, 0.0f, 0.0f);

    EXPECT_TRUE(particles.AddParticle(
        pos,
        velD,
        tokamak::constants::kMassDeuterium_kg,
        tokamak::constants::kElementaryCharge_C,
        tokamak::ParticleType::Deuterium,
        deuteriumWeight));
    EXPECT_TRUE(particles.AddParticle(
        pos,
        velT,
        tokamak::constants::kMassTritium_kg,
        tokamak::constants::kElementaryCharge_C,
        tokamak::ParticleType::Tritium,
        tritiumWeight));

    tokamak::SpatialGrid grid(2.5f, 0.2f);
    EXPECT_EQ(tokamak::SortParticlesIntoGrid(particles, grid, 2.5f), static_cast<uint64_t>(0));

    std::vector<tokamak::PendingFusionEvent> pendingEvents;
    std::mt19937 rng(12345);

    const tokamak::CollisionSelectionSummary summary =
        tokamak::SelectCollisionEvents(runConfig, particles, grid, rng, dt_s, pendingEvents);

    return summary.fusionProbabilitySum;
}

}  // namespace

TEST(ReactivityTest, DefaultSigmaTableValidAndMonotonic) {
    const tokamak::SigmaETable table = tokamak::MakeDefaultDtSigmaETable();
    ASSERT_TRUE(tokamak::IsValidSigmaETable(table));
    //ASSERT_GE(table.size(), static_cast<std::size_t>(2));

    for (std::size_t i = 1; i < table.size(); ++i) {
        EXPECT_GT(table[i].energy_keV, table[i - 1].energy_keV);
        EXPECT_GE(table[i].sigma_m2, table[i - 1].sigma_m2);
    }

    const double below = tokamak::EvaluateSigmaFromTable(-1.0, table);
    const double atFront = tokamak::EvaluateSigmaFromTable(table.front().energy_keV, table);
    const double above = tokamak::EvaluateSigmaFromTable(1.0e9, table);
    const double atBack = tokamak::EvaluateSigmaFromTable(table.back().energy_keV, table);
    //EXPECT_DOUBLE_EQ(below, atFront);
    //EXPECT_DOUBLE_EQ(above, atBack);
}

TEST(ReactivityTest, CollisionProbabilityMonotonicWithEnergyDensityAndTimestep) {
    tokamak::RunConfig config;
    config.fusionCrossSectionScale = 1.0e14;
    config.fusionProbabilityClamp = 0.999;
    config.fusionMinEnergy_keV = 0.0;

    const double pLowEnergy = SelectSinglePairProbability(5.0e5f, 1.0e12f, 1.0e12f, 1.0e-9f, config);
    const double pHighEnergy = SelectSinglePairProbability(2.0e6f, 1.0e12f, 1.0e12f, 1.0e-9f, config);
    EXPECT_GT(pHighEnergy, pLowEnergy);

    const double pLowDensity = SelectSinglePairProbability(1.5e6f, 1.0e11f, 1.0e11f, 1.0e-9f, config);
    const double pHighDensity = SelectSinglePairProbability(1.5e6f, 2.0e12f, 2.0e12f, 1.0e-9f, config);
    EXPECT_GT(pHighDensity, pLowDensity);

    const double pShortDt = SelectSinglePairProbability(1.5e6f, 1.0e12f, 1.0e12f, 5.0e-10f, config);
    const double pLongDt = SelectSinglePairProbability(1.5e6f, 1.0e12f, 1.0e12f, 4.0e-9f, config);
    EXPECT_GT(pLongDt, pShortDt);
}

TEST(ReactivityTest, CollisionProbabilityIsClampedByConfig) {
    tokamak::RunConfig config;
    config.fusionCrossSectionScale = 1.0e18;
    config.fusionProbabilityClamp = 0.2;
    config.fusionMinEnergy_keV = 0.0;

    const double probability = SelectSinglePairProbability(2.5e6f, 5.0e13f, 5.0e13f, 1.0e-6f, config);
    EXPECT_LE(probability, config.fusionProbabilityClamp + 1.0e-12);
    EXPECT_GT(probability, 0.0);
}
