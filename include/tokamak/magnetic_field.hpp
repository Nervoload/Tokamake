#pragma once

#include "tokamak/config.hpp"

namespace tokamak {

struct PlasmaCurrentProfileConfig {
    PlasmaCurrentProfileKind kind = PlasmaCurrentProfileKind::Uniform;
    float axisEpsilon_m = 1.0e-4f;
    float customAxisBlendRadius_m = 1.0e-3f;
    std::vector<CurrentProfilePoint> customTable;
};

struct MagneticFieldSample {
    Vec3 totalField_T;
    Vec3 toroidalField_T;
    Vec3 poloidalField_T;

    float totalMagnitude_T = 0.0f;
    float majorRadius_m = 0.0f;
    float minorRadius_m = 0.0f;
    float normalizedMinorRadius = 0.0f;
    float enclosedCurrentFraction = 0.0f;
};

float ComputeEnclosedCurrentFraction(
    float normalizedMinorRadius,
    const PlasmaCurrentProfileConfig& profileConfig);

MagneticFieldSample EvaluateMagneticFieldSample(
    const TokamakConfig& config,
    const PlasmaCurrentProfileConfig& profileConfig,
    const Vec3& position_m);

double RecommendDtFromMaxField(double maxField_T, double safetyFraction = 0.05);

}  // namespace tokamak
