#include "tokamak/collision.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

#include "tokamak/reactivity.hpp"

namespace tokamak {

CollisionSelectionSummary SelectCollisionEvents(
    const RunConfig& runConfig,
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

    const auto& velocities = particles.Velocities();
    const auto& positions = particles.Positions();
    const auto& species = particles.Species();
    const auto& weights = particles.Weights();

    const double clampedDt_s = std::max(0.0, static_cast<double>(dt_s));
    const double cellVolume_m3 =
        std::max(1.0e-18, static_cast<double>(grid.cellSize_m) * static_cast<double>(grid.cellSize_m) *
                              static_cast<double>(grid.cellSize_m));
    const double sigmaScale = std::max(0.0, runConfig.fusionCrossSectionScale);
    const double energyThreshold_keV = std::max(0.0, runConfig.fusionMinEnergy_keV);
    const double probabilityClamp = std::max(0.0, std::min(1.0, runConfig.fusionProbabilityClamp));

    for (int c = 0; c < grid.totalCells; ++c) {
        const uint32_t count = cellCounts[static_cast<std::size_t>(c)];
        if (count < 2) {
            continue;
        }

        const uint32_t start = cellOffsets[static_cast<std::size_t>(c)];
        const std::size_t rangeEnd = static_cast<std::size_t>(start + count);
        assert(rangeEnd <= sortedParticleIDs.size());
        const auto beginOffset = static_cast<std::vector<uint32_t>::difference_type>(start);
        const auto endOffset = static_cast<std::vector<uint32_t>::difference_type>(rangeEnd);

        std::shuffle(
            sortedParticleIDs.begin() + beginOffset,
            sortedParticleIDs.begin() + endOffset,
            rng);

        uint32_t reactionsInCell = 0;

        for (uint32_t i = 0; i + 1 < count; i += 2) {
            const std::size_t p1 = sortedParticleIDs[static_cast<std::size_t>(start + i)];
            const std::size_t p2 = sortedParticleIDs[static_cast<std::size_t>(start + i + 1)];
            assert(p1 < particles.Size());
            assert(p2 < particles.Size());

            const auto species1 = species[p1];
            const auto species2 = species[p2];
            const bool isDT =
                (species1 == ParticleType::Deuterium && species2 == ParticleType::Tritium) ||
                (species1 == ParticleType::Tritium && species2 == ParticleType::Deuterium);
            if (!isDT) {
                continue;
            }

            const double reactedWeight = std::max(0.0, std::min(weights[p1], weights[p2]));
            if (reactedWeight <= 0.0) {
                continue;
            }

            const Vec3 vRel = velocities[p1] - velocities[p2];
            const double relativeSpeed_mPerS = static_cast<double>(vRel.Magnitude());

            const double reducedMass =
                (constants::kMassDeuterium_kg * constants::kMassTritium_kg) /
                (constants::kMassDeuterium_kg + constants::kMassTritium_kg);
            const double eKinetic_J = 0.5 * reducedMass * static_cast<double>(Vec3::Dot(vRel, vRel));
            const double eKinetic_keV = eKinetic_J / (1000.0 * constants::kElementaryCharge_C);

            const double sigma_m2 =
                (eKinetic_keV < energyThreshold_keV)
                    ? 0.0
                    : (sigmaScale * EvaluateDtSigma_m2(eKinetic_keV, runConfig.fusionReactivityModelKind));
            if (sigma_m2 <= 0.0) {
                continue;
            }

            ++summary.fusionAttempts;
            ++summary.fusionKineticsSamples;
            summary.fusionWeightAttempted += reactedWeight;

            const double n_eff = reactedWeight / cellVolume_m3;
            const double lambda = sigma_m2 * relativeSpeed_mPerS * n_eff;
            const double rawProbability = 1.0 - std::exp(-(lambda * clampedDt_s));
            const double clampedProbability = std::max(0.0, std::min(probabilityClamp, rawProbability));

            summary.fusionSigmaSum_m2 += sigma_m2;
            summary.fusionProbabilitySum += clampedProbability;
            summary.fusionRelativeSpeedSum_mPerS += relativeSpeed_mPerS;

            if (dist01(rng) >= clampedProbability) {
                continue;
            }

            Vec3 randomDir(dist01(rng) - 0.5f, dist01(rng) - 0.5f, dist01(rng) - 0.5f);
            if (randomDir.Magnitude() < 1.0e-6f) {
                randomDir = Vec3(1.0f, 0.0f, 0.0f);
            }

            const Vec3 centerOfMassPos = (positions[p1] + positions[p2]) * 0.5f;
            const double mass1 =
                (species1 == ParticleType::Deuterium) ? constants::kMassDeuterium_kg : constants::kMassTritium_kg;
            const double mass2 =
                (species2 == ParticleType::Deuterium) ? constants::kMassDeuterium_kg : constants::kMassTritium_kg;
            const Vec3 centerOfMassVel =
                (velocities[p1] * static_cast<float>(mass1) + velocities[p2] * static_cast<float>(mass2)) /
                static_cast<float>(mass1 + mass2);
            const double alphaEnergy_J = 3.5e6 * constants::kElementaryCharge_C;
            const float alphaSpeed = static_cast<float>(std::sqrt((2.0 * alphaEnergy_J) / constants::kMassHelium4_kg));
            const Vec3 heliumVelocity = centerOfMassVel + (randomDir.Normalized() * alphaSpeed);

            pendingEvents.push_back(PendingFusionEvent{
                p1,
                p2,
                centerOfMassPos,
                heliumVelocity,
                reactedWeight,
                relativeSpeed_mPerS,
                sigma_m2,
                clampedProbability,
            });
            ++reactionsInCell;
            summary.fusionWeightAccepted += reactedWeight;
        }

        if (reactionsInCell > summary.maxReactionsInCell) {
            summary.maxReactionsInCell = reactionsInCell;
        }
    }

    summary.selectedEvents = static_cast<uint64_t>(pendingEvents.size());
    return summary;
}

CollisionSelectionSummary SelectCollisionEvents(
    ParticleSystem& particles,
    SpatialGrid& grid,
    std::mt19937& rng,
    float dt_s,
    std::vector<PendingFusionEvent>& pendingEvents) {
    RunConfig defaults;
    return SelectCollisionEvents(defaults, particles, grid, rng, dt_s, pendingEvents);
}

void ApplyCollisionEvents(
    ParticleSystem& particles,
    const std::vector<PendingFusionEvent>& pendingEvents,
    RuntimeCounters& counters,
    EnergyChargeBudget& budget,
    std::vector<Vec3>* acceptedFusionPositions,
    std::vector<double>* acceptedFusionWeights,
    CollisionKineticsSnapshot* outKinetics) {
    if (acceptedFusionPositions != nullptr) {
        acceptedFusionPositions->clear();
    }
    if (acceptedFusionWeights != nullptr) {
        acceptedFusionWeights->clear();
    }
    if (outKinetics != nullptr) {
        *outKinetics = CollisionKineticsSnapshot{};
    }

    auto& weights = particles.MutableWeights();
    auto& species = particles.MutableSpecies();

    for (const auto& event : pendingEvents) {
        if (event.p1 >= particles.Size() || event.p2 >= particles.Size()) {
            continue;
        }

        const auto species1 = species[event.p1];
        const auto species2 = species[event.p2];
        const bool isDT =
            (species1 == ParticleType::Deuterium && species2 == ParticleType::Tritium) ||
            (species1 == ParticleType::Tritium && species2 == ParticleType::Deuterium);
        if (!isDT) {
            continue;
        }

        double requestedReactedWeight = event.reactedWeight;
        if (requestedReactedWeight <= 0.0) {
            // Backward-compatible default for callers that do not set weighted-event fields.
            requestedReactedWeight = std::min(weights[event.p1], weights[event.p2]);
        }
        const double reactedWeight =
            std::max(0.0, std::min(requestedReactedWeight, std::min(weights[event.p1], weights[event.p2])));
        if (reactedWeight <= 0.0) {
            continue;
        }

        constexpr double kWeightEpsilon = 1.0e-9;
        const double remainingWeight1 = weights[event.p1] - reactedWeight;
        const double remainingWeight2 = weights[event.p2] - reactedWeight;
        const bool consumeP1 = remainingWeight1 <= kWeightEpsilon;
        const bool consumeP2 = remainingWeight2 <= kWeightEpsilon;

        if (!particles.CanInsert(1) && !consumeP1 && !consumeP2) {
            ++counters.rejectedFusionAsh;
            ++counters.particleCapHitEvents;
            continue;
        }

        if (!consumeP1) {
            weights[event.p1] = remainingWeight1;
        } else {
            particles.MarkDead(event.p1);
        }
        if (!consumeP2) {
            weights[event.p2] = remainingWeight2;
        } else {
            particles.MarkDead(event.p2);
        }

        const bool added = particles.AddParticle(
            event.centerOfMassPos,
            event.heliumVelocity,
            constants::kMassHelium4_kg,
            constants::kElementaryCharge_C * 2.0,
            ParticleType::Helium,
            reactedWeight);

        if (!added) {
            ++counters.rejectedFusionAsh;
            ++counters.particleCapHitEvents;
            if (!consumeP1) {
                weights[event.p1] += reactedWeight;
            }
            if (!consumeP2) {
                weights[event.p2] += reactedWeight;
            }
            if (consumeP1) {
                particles.RestoreSpecies(event.p1, species1);
            }
            if (consumeP2) {
                particles.RestoreSpecies(event.p2, species2);
            }
            continue;
        }

        ++particles.fusionCountTotal;
        ++counters.fusionAccepted;
        counters.fusionWeightAccepted += static_cast<double>(reactedWeight);
        counters.ashWeightProducedHe += static_cast<double>(reactedWeight);
        counters.fuelWeightConsumedD += static_cast<double>(reactedWeight);
        counters.fuelWeightConsumedT += static_cast<double>(reactedWeight);

        budget.fusionAlphaInjected_J += reactedWeight * (3.5e6 * constants::kElementaryCharge_C);
        if (acceptedFusionPositions != nullptr) {
            acceptedFusionPositions->push_back(event.centerOfMassPos);
        }
        if (acceptedFusionWeights != nullptr) {
            acceptedFusionWeights->push_back(reactedWeight);
        }

        if (outKinetics != nullptr) {
            ++outKinetics->samples;
            outKinetics->acceptedWeight += static_cast<double>(reactedWeight);
            outKinetics->avgSigma_m2 += static_cast<double>(event.sigma_m2);
            outKinetics->avgProbability += static_cast<double>(event.probability);
            outKinetics->avgRelativeSpeed_mPerS += static_cast<double>(event.relativeSpeed_mPerS);
        }
    }

    if (outKinetics != nullptr && outKinetics->samples > 0) {
        const double invSamples = 1.0 / static_cast<double>(outKinetics->samples);
        outKinetics->avgSigma_m2 *= invSamples;
        outKinetics->avgProbability *= invSamples;
        outKinetics->avgRelativeSpeed_mPerS *= invSamples;
    }
}

}  // namespace tokamak
