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

TEST(ParticleSystemTest, ParticleIdsRemainStableAcrossCompaction) {
    tokamak::ParticleSystem particles(4);

    ASSERT_TRUE(particles.AddParticle(
        ZeroVec(),
        ZeroVec(),
        tokamak::constants::kMassDeuterium_kg,
        tokamak::constants::kElementaryCharge_C,
        tokamak::ParticleType::Deuterium));
    ASSERT_TRUE(particles.AddParticle(
        tokamak::Vec3(1.0f, 0.0f, 0.0f),
        ZeroVec(),
        tokamak::constants::kMassTritium_kg,
        tokamak::constants::kElementaryCharge_C,
        tokamak::ParticleType::Tritium));
    ASSERT_TRUE(particles.AddParticle(
        tokamak::Vec3(2.0f, 0.0f, 0.0f),
        ZeroVec(),
        tokamak::constants::kMassHelium4_kg,
        tokamak::constants::kElementaryCharge_C,
        tokamak::ParticleType::Helium));

    const std::vector<uint64_t> idsBefore = particles.Ids();
    ASSERT_EQ(idsBefore.size(), static_cast<std::size_t>(3));

    particles.MarkDead(1);
    particles.Compact();

    ASSERT_EQ(particles.Size(), static_cast<std::size_t>(2));
    EXPECT_TRUE(particles.IsArrayLengthConsistent());
    EXPECT_EQ(particles.Ids()[0], idsBefore[0]);
    EXPECT_EQ(particles.Ids()[1], idsBefore[2]);
}

TEST(ParticleSystemTest, AddParticleReusesDeadSlotBeforeCompaction) {
    tokamak::ParticleSystem particles(2);

    ASSERT_TRUE(particles.AddParticle(
        ZeroVec(),
        ZeroVec(),
        tokamak::constants::kMassDeuterium_kg,
        tokamak::constants::kElementaryCharge_C,
        tokamak::ParticleType::Deuterium));
    ASSERT_TRUE(particles.AddParticle(
        tokamak::Vec3(1.0f, 0.0f, 0.0f),
        ZeroVec(),
        tokamak::constants::kMassTritium_kg,
        tokamak::constants::kElementaryCharge_C,
        tokamak::ParticleType::Tritium));

    particles.MarkDead(0);
    EXPECT_TRUE(particles.CanInsert(1));
    EXPECT_TRUE(particles.AddParticle(
        tokamak::Vec3(2.0f, 0.0f, 0.0f),
        ZeroVec(),
        tokamak::constants::kMassHelium4_kg,
        tokamak::constants::kElementaryCharge_C,
        tokamak::ParticleType::Helium));
    EXPECT_EQ(particles.Size(), static_cast<std::size_t>(2));
    EXPECT_EQ(particles.Species()[0], tokamak::ParticleType::Helium);
    EXPECT_EQ(particles.Species()[1], tokamak::ParticleType::Tritium);
}
