#pragma once

#include <cstdint>

#include "math_types.h"

enum class Scenario : uint8_t {
    COLD_VACUUM = 0,
    NBI_IGNITION = 1,
    MAGNETIC_FAILURE = 2,
};

enum class ParticleType : uint8_t {
    Deuterium = 0,
    Tritium = 1,
    Helium = 2,
    Dead = 3,
};

struct TokamakConfig {
    float majorRadius = 2.0f;
    float minorRadius = 0.5f;
    float toroidalCurrent = 15.0e6f;
    int32_t toroidalCoilTurns = 18;
    float plasmaCurrent = 2.0e6f;
};

struct NBIConfig {
    bool isActive = true;
    float beamEnergyKeV = 100.0f;
    int32_t particlesPerStep = 50;
    Vec3 injectorPos = Vec3(2.5f, 0.0f, 0.0f);
    Vec3 injectionNormal = Vec3(-1.0f, 0.2f, 0.0f).normalize();
};
