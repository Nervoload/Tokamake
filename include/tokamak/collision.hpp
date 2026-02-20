#pragma once

#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "tokamak/diagnostics.hpp"
#include "tokamak/spatial_grid.hpp"

namespace tokamak {

struct PendingFusionEvent {
    std::size_t p1 = 0;
    std::size_t p2 = 0;
    Vec3 centerOfMassPos;
    Vec3 heliumVelocity;
};

struct CollisionSelectionSummary {
    uint64_t fusionAttempts = 0;
    uint64_t selectedEvents = 0;
    uint32_t maxReactionsInCell = 0;
};

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
    std::vector<Vec3>* acceptedFusionPositions = nullptr);

}  // namespace tokamak
