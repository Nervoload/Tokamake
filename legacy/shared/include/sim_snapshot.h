#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "tokamak_types.h"

struct SnapshotRequest {
    uint32_t maxSampledParticles = 100000;
    bool includeVelocities = true;
    bool includeMassAndCharge = true;
};

struct ParticleSample {
    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;
    float vx = 0.0f;
    float vy = 0.0f;
    float vz = 0.0f;
    float mass = 0.0f;
    float charge = 0.0f;
    float qOverM = 0.0f;
    uint8_t species = static_cast<uint8_t>(ParticleType::Dead);
};

struct TelemetrySummary {
    uint32_t totalIons = 0;
    uint32_t deuteriumCount = 0;
    uint32_t tritiumCount = 0;
    uint32_t heliumCount = 0;
    double avgTempKeV = 0.0;
    uint64_t fusionEvents = 0;
};

struct SnapshotPerformance {
    double simStepMs = 0.0;
    double exportMs = 0.0;
    double frameSubmitMs = 0.0;
};

struct SimulationSnapshot {
    uint64_t stepIndex = 0;
    double simTimeSeconds = 0.0;
    TokamakConfig config;
    NBIConfig nbi;

    uint32_t totalParticleCount = 0;
    uint32_t sampleStride = 1;
    TelemetrySummary telemetry;
    SnapshotPerformance performance;

    std::vector<ParticleSample> sampledParticles;

    void Reserve(size_t count) { sampledParticles.reserve(count); }
};
