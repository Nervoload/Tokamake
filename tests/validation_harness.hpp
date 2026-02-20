#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "tokamak/particle_system.hpp"
#include "tokamak/types.hpp"

namespace tokamak_validation {

struct ScalarTolerance {
    double absolute;
    double relative;
    const char* rationale;
};

inline double WrappedAngleDelta(double previous, double current) {
    constexpr double kTwoPi = 2.0 * 3.14159265358979323846;
    double delta = current - previous;
    while (delta > 3.14159265358979323846) {
        delta -= kTwoPi;
    }
    while (delta < -3.14159265358979323846) {
        delta += kTwoPi;
    }
    return delta;
}

inline void ExpectNearWithTolerance(
    const std::string& label,
    double observed,
    double expected,
    const ScalarTolerance& tolerance) {
    const double absError = std::abs(observed - expected);
    const double allowed =
        tolerance.absolute + (tolerance.relative * std::max(std::abs(expected), std::numeric_limits<double>::epsilon()));
    if (absError > allowed) {
        std::ostringstream oss;
        oss << label << " observed=" << observed << " expected=" << expected << " abs_error=" << absError
            << " allowed=" << allowed << " rationale=" << tolerance.rationale;
        ::testing::AddFailure(__FILE__, __LINE__, oss.str(), false);
    }
}

inline std::vector<tokamak::Vec3> MakeDeterministicPositions(
    std::size_t count,
    float minCoordinate_m,
    float maxCoordinate_m,
    uint32_t seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> dist(minCoordinate_m, maxCoordinate_m);

    std::vector<tokamak::Vec3> positions;
    positions.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        positions.emplace_back(dist(rng), dist(rng), dist(rng));
    }
    return positions;
}

inline tokamak::ParticleSystem MakeHotDtCollisionPool(std::size_t pairsPerSpecies) {
    tokamak::ParticleSystem particles((pairsPerSpecies * 2U) + 8U);
    const tokamak::Vec3 velocityD(2.2e6f, 0.0f, 0.0f);
    const tokamak::Vec3 velocityT(-2.2e6f, 0.0f, 0.0f);

    for (std::size_t i = 0; i < pairsPerSpecies; ++i) {
        const float jitter = 1.0e-4f * static_cast<float>(i % 11U);
        const tokamak::Vec3 position(2.0f + jitter, jitter, -jitter);
        EXPECT_TRUE(particles.AddParticle(
            position,
            velocityD,
            tokamak::constants::kMassDeuterium_kg,
            tokamak::constants::kElementaryCharge_C,
            tokamak::ParticleType::Deuterium));
        EXPECT_TRUE(particles.AddParticle(
            position,
            velocityT,
            tokamak::constants::kMassTritium_kg,
            tokamak::constants::kElementaryCharge_C,
            tokamak::ParticleType::Tritium));
    }

    return particles;
}

inline bool SortedIdsFormCompletePermutation(const std::vector<uint32_t>& sortedIds, std::size_t expectedSize) {
    if (sortedIds.size() < expectedSize) {
        return false;
    }

    std::vector<uint8_t> seen(expectedSize, 0);
    for (std::size_t i = 0; i < expectedSize; ++i) {
        const uint32_t id = sortedIds[i];
        if (id >= expectedSize || seen[id] != 0U) {
            return false;
        }
        seen[id] = 1U;
    }
    return true;
}

}  // namespace tokamak_validation
