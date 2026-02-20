#include "tokamak/magnetic_field.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace tokamak {

namespace {

constexpr float kTwoPi = 2.0f * constants::kPi;

float Clamp01(float value) {
    return std::max(0.0f, std::min(value, 1.0f));
}

float SmoothStep01(float t) {
    const float clamped = Clamp01(t);
    return clamped * clamped * (3.0f - (2.0f * clamped));
}

std::vector<CurrentProfilePoint> SanitizeCustomTable(const std::vector<CurrentProfilePoint>& table) {
    std::vector<CurrentProfilePoint> points;
    points.reserve(table.size() + 2);

    for (const auto& point : table) {
        points.push_back(CurrentProfilePoint{
            Clamp01(point.normalizedMinorRadius),
            Clamp01(point.enclosedCurrentFraction),
        });
    }

    if (points.empty()) {
        points.push_back(CurrentProfilePoint{0.0f, 0.0f});
        points.push_back(CurrentProfilePoint{1.0f, 1.0f});
        return points;
    }

    std::sort(
        points.begin(),
        points.end(),
        [](const CurrentProfilePoint& a, const CurrentProfilePoint& b) {
            return a.normalizedMinorRadius < b.normalizedMinorRadius;
        });

    std::vector<CurrentProfilePoint> deduped;
    deduped.reserve(points.size() + 2);

    constexpr float kMergeEpsilon = 1.0e-6f;
    for (const auto& point : points) {
        if (!deduped.empty() &&
            std::fabs(point.normalizedMinorRadius - deduped.back().normalizedMinorRadius) <= kMergeEpsilon) {
            deduped.back().enclosedCurrentFraction =
                std::max(deduped.back().enclosedCurrentFraction, point.enclosedCurrentFraction);
            continue;
        }
        deduped.push_back(point);
    }

    if (deduped.front().normalizedMinorRadius > 0.0f) {
        deduped.insert(deduped.begin(), CurrentProfilePoint{0.0f, 0.0f});
    } else {
        deduped.front().normalizedMinorRadius = 0.0f;
        deduped.front().enclosedCurrentFraction = 0.0f;
    }

    if (deduped.back().normalizedMinorRadius < 1.0f) {
        deduped.push_back(CurrentProfilePoint{1.0f, 1.0f});
    } else {
        deduped.back().normalizedMinorRadius = 1.0f;
        deduped.back().enclosedCurrentFraction = 1.0f;
    }

    float running = 0.0f;
    for (auto& point : deduped) {
        running = std::max(running, point.enclosedCurrentFraction);
        point.enclosedCurrentFraction = running;
    }
    deduped.back().enclosedCurrentFraction = 1.0f;
    return deduped;
}

float InterpolateCustomFraction(
    float rho,
    const std::vector<CurrentProfilePoint>& sanitizedTable) {
    if (rho <= 0.0f) {
        return 0.0f;
    }
    if (rho >= 1.0f) {
        return 1.0f;
    }
    if (sanitizedTable.size() < 2) {
        return rho * rho;
    }

    for (std::size_t i = 1; i < sanitizedTable.size(); ++i) {
        const auto& left = sanitizedTable[i - 1];
        const auto& right = sanitizedTable[i];
        if (rho <= right.normalizedMinorRadius) {
            const float width = std::max(1.0e-6f, right.normalizedMinorRadius - left.normalizedMinorRadius);
            const float t = (rho - left.normalizedMinorRadius) / width;
            return Clamp01(left.enclosedCurrentFraction +
                           ((right.enclosedCurrentFraction - left.enclosedCurrentFraction) * t));
        }
    }

    return 1.0f;
}

float ComputeEnclosedCurrentFractionCustom(
    float rho,
    const PlasmaCurrentProfileConfig& profileConfig) {
    const std::vector<CurrentProfilePoint> sanitized = SanitizeCustomTable(profileConfig.customTable);
    return InterpolateCustomFraction(rho, sanitized);
}

}  // namespace

float ComputeEnclosedCurrentFraction(
    float normalizedMinorRadius,
    const PlasmaCurrentProfileConfig& profileConfig) {
    const float rho = Clamp01(normalizedMinorRadius);
    if (rho >= 1.0f) {
        return 1.0f;
    }

    switch (profileConfig.kind) {
        case PlasmaCurrentProfileKind::Uniform:
            return rho * rho;
        case PlasmaCurrentProfileKind::Parabolic: {
            const float rhoSq = rho * rho;
            return Clamp01((2.0f * rhoSq) - (rhoSq * rhoSq));
        }
        case PlasmaCurrentProfileKind::CustomTable:
            return ComputeEnclosedCurrentFractionCustom(rho, profileConfig);
    }

    return rho * rho;
}

MagneticFieldSample EvaluateMagneticFieldSample(
    const TokamakConfig& config,
    const PlasmaCurrentProfileConfig& profileConfig,
    const Vec3& position_m) {
    MagneticFieldSample sample;

    const float axisEpsilon_m = std::max(1.0e-7f, profileConfig.axisEpsilon_m);
    const float R = std::sqrt((position_m.x * position_m.x) + (position_m.y * position_m.y));
    const float RSafe = std::max(R, axisEpsilon_m);
    sample.majorRadius_m = R;

    float xHat = 0.0f;
    float yHat = 0.0f;
    if (R > axisEpsilon_m) {
        xHat = position_m.x / R;
        yHat = position_m.y / R;
    }

    const float bTorMagnitude_T =
        (constants::kMu0 * static_cast<float>(config.toroidalCoilTurns) * config.toroidalCurrent_A) /
        (kTwoPi * RSafe);
    sample.toroidalField_T = Vec3(-bTorMagnitude_T * yHat, bTorMagnitude_T * xHat, 0.0f);

    const float dR = R - config.majorRadius_m;
    const float rMinor = std::sqrt((dR * dR) + (position_m.z * position_m.z));
    sample.minorRadius_m = rMinor;

    const float minorRadiusSafe = std::max(config.minorRadius_m, axisEpsilon_m);
    const float rho = Clamp01(rMinor / minorRadiusSafe);
    sample.normalizedMinorRadius = rho;

    float enclosedFraction = ComputeEnclosedCurrentFraction(rho, profileConfig);
    if (profileConfig.kind == PlasmaCurrentProfileKind::CustomTable && profileConfig.customAxisBlendRadius_m > 0.0f) {
        const float blendRho = Clamp01(profileConfig.customAxisBlendRadius_m / minorRadiusSafe);
        if (blendRho > 0.0f && rho < blendRho) {
            const float nearAxisQuadratic = rho * rho;
            const float t = SmoothStep01(rho / blendRho);
            enclosedFraction = Clamp01(nearAxisQuadratic + ((enclosedFraction - nearAxisQuadratic) * t));
        }
    }
    sample.enclosedCurrentFraction = enclosedFraction;
    const float enclosedCurrent_A = config.plasmaCurrent_A * enclosedFraction;

    const float rSafe = std::max(rMinor, axisEpsilon_m);
    const float bPolMagnitude_T = (constants::kMu0 * enclosedCurrent_A) / (kTwoPi * rSafe);
    const float bR = -bPolMagnitude_T * (position_m.z / rSafe);
    const float bZ = bPolMagnitude_T * (dR / rSafe);
    sample.poloidalField_T = Vec3(bR * xHat, bR * yHat, bZ);

    sample.totalField_T = sample.toroidalField_T + sample.poloidalField_T;
    sample.totalMagnitude_T = sample.totalField_T.Magnitude();
    return sample;
}

double RecommendDtFromMaxField(double maxField_T, double safetyFraction) {
    if (!std::isfinite(maxField_T) || maxField_T <= 0.0) {
        return std::numeric_limits<double>::infinity();
    }

    const double clampedSafety = std::max(0.0, safetyFraction);
    const double qOverMDeuterium =
        static_cast<double>(constants::kElementaryCharge_C) / static_cast<double>(constants::kMassDeuterium_kg);
    const double qOverMTritium =
        static_cast<double>(constants::kElementaryCharge_C) / static_cast<double>(constants::kMassTritium_kg);
    const double qOverMHelium =
        static_cast<double>(2.0f * constants::kElementaryCharge_C) /
        static_cast<double>(constants::kMassHelium4_kg);
    const double maxQOverM = std::max(qOverMDeuterium, std::max(qOverMTritium, qOverMHelium));

    const double maxGyroOmega = maxQOverM * maxField_T;
    if (!std::isfinite(maxGyroOmega) || maxGyroOmega <= 0.0) {
        return std::numeric_limits<double>::infinity();
    }

    return clampedSafety / maxGyroOmega;
}

}  // namespace tokamak
