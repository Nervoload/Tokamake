#include <cstddef>
#include <limits>

#include "gtest/gtest.h"
#include "tokamak/particle_system.hpp"

namespace {

tokamak::Vec3 ZeroVec() {
    return tokamak::Vec3(0.0f, 0.0f, 0.0f);
}

}  // namespace

TEST(ParticleSystemTest, AddParticleReturnsFalseAtCapacity) {
    tokamak::ParticleSystem particles(2);

    EXPECT_TRUE(particles.AddParticle(
        ZeroVec(),
        ZeroVec(),
        tokamak::constants::kMassDeuterium_kg,
        tokamak::constants::kElementaryCharge_C,
        tokamak::ParticleType::Deuterium));
    EXPECT_TRUE(particles.AddParticle(
        ZeroVec(),
        ZeroVec(),
        tokamak::constants::kMassTritium_kg,
        tokamak::constants::kElementaryCharge_C,
        tokamak::ParticleType::Tritium));

    EXPECT_FALSE(particles.AddParticle(
        ZeroVec(),
        ZeroVec(),
        tokamak::constants::kMassHelium4_kg,
        tokamak::constants::kElementaryCharge_C,
        tokamak::ParticleType::Helium));

    EXPECT_EQ(particles.Size(), static_cast<std::size_t>(2));
    EXPECT_TRUE(particles.IsArrayLengthConsistent());
}

TEST(ParticleSystemTest, AddParticleRejectsNonFiniteVectors) {
    tokamak::ParticleSystem particles(4);
    const float nan = std::numeric_limits<float>::quiet_NaN();
    const float inf = std::numeric_limits<float>::infinity();

    EXPECT_FALSE(particles.AddParticle(
        tokamak::Vec3(nan, 0.0f, 0.0f),
        ZeroVec(),
        tokamak::constants::kMassDeuterium_kg,
        tokamak::constants::kElementaryCharge_C,
        tokamak::ParticleType::Deuterium));

    EXPECT_FALSE(particles.AddParticle(
        ZeroVec(),
        tokamak::Vec3(0.0f, inf, 0.0f),
        tokamak::constants::kMassTritium_kg,
        tokamak::constants::kElementaryCharge_C,
        tokamak::ParticleType::Tritium));

    EXPECT_EQ(particles.Size(), static_cast<std::size_t>(0));
    EXPECT_TRUE(particles.IsArrayLengthConsistent());
}
