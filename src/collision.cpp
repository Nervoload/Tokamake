#include "tokamak/collision.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

namespace tokamak {

CollisionSelectionSummary SelectCollisionEvents(
    ParticleSystem& particles,
    SpatialGrid& grid,
    std::mt19937& rng,
    float dt_s,
    std::vector<PendingFusionEvent>& pendingEvents) {
    pendingEvents.clear();
    CollisionSelectionSummary summary;

    std::uniform_real_distribution<float> dist01(0.0f, 1.0f);

    auto& sortedParticleIDs = grid.sortedParticleIDs;
    const auto& cellCounts = grid.cellCounts;
    const auto& cellOffsets = grid.cellOffsets;

    for (int c = 0; c < grid.totalCells; ++c) {
        const uint32_t count = cellCounts[static_cast<std::size_t>(c)];
        if (count < 2) {
            continue;
        }

        const uint32_t start = cellOffsets[static_cast<std::size_t>(c)];
        const std::size_t rangeEnd = static_cast<std::size_t>(start + count);
        assert(rangeEnd <= sortedParticleIDs.size());

        std::shuffle(
            sortedParticleIDs.begin() + static_cast<std::size_t>(start),
            sortedParticleIDs.begin() + rangeEnd,
            rng);

        uint32_t reactionsInCell = 0;

        for (uint32_t i = 0; i + 1 < count; i += 2) {
            const std::size_t p1 = sortedParticleIDs[static_cast<std::size_t>(start + i)];
            const std::size_t p2 = sortedParticleIDs[static_cast<std::size_t>(start + i + 1)];
            assert(p1 < particles.Size());
            assert(p2 < particles.Size());

            const auto species1 = particles.Species()[p1];
            const auto species2 = particles.Species()[p2];
            const bool isDT =
                (species1 == ParticleType::Deuterium && species2 == ParticleType::Tritium) ||
                (species1 == ParticleType::Tritium && species2 == ParticleType::Deuterium);
            if (!isDT) {
                continue;
            }

            const Vec3 vRel = particles.Velocities()[p1] - particles.Velocities()[p2];
            const float reducedMass =
                (constants::kMassDeuterium_kg * constants::kMassTritium_kg) /
                (constants::kMassDeuterium_kg + constants::kMassTritium_kg);
            const float eKinetic_J = 0.5f * reducedMass * Vec3::Dot(vRel, vRel);
            const float eKinetic_keV = eKinetic_J / (1000.0f * constants::kElementaryCharge_C);

            if (eKinetic_keV <= 15.0f) {
                continue;
            }

            ++summary.fusionAttempts;
            const float lambda = std::max(0.0f, 0.05f * particles.macroWeight);
            const float probability = 1.0f - std::exp(-(lambda * std::max(0.0f, dt_s)));
            if (dist01(rng) >= probability) {
                continue;
            }

            Vec3 randomDir(dist01(rng) - 0.5f, dist01(rng) - 0.5f, dist01(rng) - 0.5f);
            if (randomDir.Magnitude() < 1.0e-6f) {
                randomDir = Vec3(1.0f, 0.0f, 0.0f);
            }

            const Vec3 centerOfMassPos = (particles.Positions()[p1] + particles.Positions()[p2]) * 0.5f;
            const float mass1 =
                (species1 == ParticleType::Deuterium) ? constants::kMassDeuterium_kg : constants::kMassTritium_kg;
            const float mass2 =
                (species2 == ParticleType::Deuterium) ? constants::kMassDeuterium_kg : constants::kMassTritium_kg;
            const Vec3 centerOfMassVel =
                (particles.Velocities()[p1] * mass1 + particles.Velocities()[p2] * mass2) /
                (mass1 + mass2);
            const float alphaEnergy_J = 3.5e6f * constants::kElementaryCharge_C;
            const float alphaSpeed = std::sqrt((2.0f * alphaEnergy_J) / constants::kMassHelium4_kg);
            const Vec3 heliumVelocity = centerOfMassVel + (randomDir.Normalized() * alphaSpeed);

            pendingEvents.push_back(PendingFusionEvent{p1, p2, centerOfMassPos, heliumVelocity});
            ++reactionsInCell;
        }

        if (reactionsInCell > summary.maxReactionsInCell) {
            summary.maxReactionsInCell = reactionsInCell;
        }
    }

    summary.selectedEvents = static_cast<uint64_t>(pendingEvents.size());
    return summary;
}

void ApplyCollisionEvents(
    ParticleSystem& particles,
    const std::vector<PendingFusionEvent>& pendingEvents,
    RuntimeCounters& counters,
    EnergyChargeBudget& budget,
    std::vector<Vec3>* acceptedFusionPositions) {
    if (acceptedFusionPositions != nullptr) {
        acceptedFusionPositions->clear();
    }

    for (const auto& event : pendingEvents) {
        if (event.p1 >= particles.Size() || event.p2 >= particles.Size()) {
            continue;
        }

        const auto species1 = particles.Species()[event.p1];
        const auto species2 = particles.Species()[event.p2];
        const bool isDT =
            (species1 == ParticleType::Deuterium && species2 == ParticleType::Tritium) ||
            (species1 == ParticleType::Tritium && species2 == ParticleType::Deuterium);
        if (!isDT) {
            continue;
        }

        if (!particles.CanInsert(1)) {
            ++counters.rejectedFusionAsh;
            ++counters.particleCapHitEvents;
            continue;
        }

        const bool added = particles.AddParticle(
            event.centerOfMassPos,
            event.heliumVelocity,
            constants::kMassHelium4_kg,
            constants::kElementaryCharge_C * 2.0f,
            ParticleType::Helium);

        if (!added) {
            ++counters.rejectedFusionAsh;
            ++counters.particleCapHitEvents;
            continue;
        }

        particles.MarkDead(event.p1);
        particles.MarkDead(event.p2);
        ++particles.fusionCountTotal;
        ++counters.fusionAccepted;
        budget.fusionAlphaInjected_J += static_cast<double>(3.5e6f * constants::kElementaryCharge_C);
        if (acceptedFusionPositions != nullptr) {
            acceptedFusionPositions->push_back(event.centerOfMassPos);
        }
    }
}

}  // namespace tokamak
