#pragma once

#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "tokamak/config.hpp"
#include "tokamak/diagnostics.hpp"
#include "tokamak/spatial_grid.hpp"

namespace tokamak {

struct PendingFusionEvent {
    std::size_t p1 = 0;
    std::size_t p2 = 0;
    Vec3 centerOfMassPos;
    Vec3 heliumVelocity;
    float reactedWeight = 0.0f;
    float relativeSpeed_mPerS = 0.0f;
    float sigma_m2 = 0.0f;
    float probability = 0.0f;
};

struct CollisionSelectionSummary {
    uint64_t fusionAttempts = 0;
    uint64_t selectedEvents = 0;
    uint32_t maxReactionsInCell = 0;
    double fusionWeightAttempted = 0.0;
    double fusionWeightAccepted = 0.0;
    double fusionSigmaSum_m2 = 0.0;
    double fusionProbabilitySum = 0.0;
    double fusionRelativeSpeedSum_mPerS = 0.0;
    uint64_t fusionKineticsSamples = 0;
};

struct CollisionKineticsSnapshot {
    uint64_t samples = 0;
    double avgSigma_m2 = 0.0;
    double avgProbability = 0.0;
    double avgRelativeSpeed_mPerS = 0.0;
    double attemptedWeight = 0.0;
    double acceptedWeight = 0.0;
};

CollisionSelectionSummary SelectCollisionEvents(
    const RunConfig& runConfig,
    ParticleSystem& particles,
    SpatialGrid& grid,
    std::mt19937& rng,
    float dt_s,
    std::vector<PendingFusionEvent>& pendingEvents);

// Compatibility overload for existing call sites that do not need run-config inputs.
CollisionSelectionSummary SelectCollisionEvents(
    ParticleSystem& particles,
    SpatialGrid& grid,
    std::mt19937& rng,
    float dt_s,
    std::vector<PendingFusionEvent>& pendingEvents);

void ApplyCollisionEvents(
    ParticleSystem& particles,
    const std::vector<PendingFusionEvent>& pendingEvents,
    RuntimeCounters& counters,
    EnergyChargeBudget& budget,
    std::vector<Vec3>* acceptedFusionPositions = nullptr,
    CollisionKineticsSnapshot* outKinetics = nullptr);

}  // namespace tokamak
