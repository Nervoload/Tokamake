#include "tokamak/viewer/viewer_app.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <utility>

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

#include "tokamak/magnetic_field.hpp"
#include "tokamak/viewer/camera.hpp"
#include "tokamak/viewer/gl_renderer.hpp"

namespace tokamak::viewer {
namespace {

struct WindowInputState {
    double scrollDeltaY = 0.0;
};

void ScrollCallback(GLFWwindow* window, double /*xoffset*/, double yoffset) {
    WindowInputState* input = static_cast<WindowInputState*>(glfwGetWindowUserPointer(window));
    if (input != nullptr) {
        input->scrollDeltaY += yoffset;
    }
}

bool ParseFloatValue(const char* text, float* outValue) {
    char* end = nullptr;
    const float value = std::strtof(text, &end);
    if (end == text || *end != '\0' || !std::isfinite(value)) {
        return false;
    }
    *outValue = value;
    return true;
}

bool ParseIntValue(const char* text, int* outValue) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value < std::numeric_limits<int>::min() ||
        value > std::numeric_limits<int>::max()) {
        return false;
    }
    *outValue = static_cast<int>(value);
    return true;
}

std::size_t ClampToSizeT(uint64_t value) {
    const uint64_t maxValue = static_cast<uint64_t>(std::numeric_limits<std::size_t>::max());
    return static_cast<std::size_t>(std::min(value, maxValue));
}

constexpr double kMinDisplayFrameDuration_s = 1.0 / 24.0;

float Clamp01(float value) {
    return std::max(0.0f, std::min(value, 1.0f));
}

double Clamp01Double(double value) {
    return std::max(0.0, std::min(value, 1.0));
}

double SmoothStep01(double value) {
    const double clamped = Clamp01Double(value);
    return clamped * clamped * (3.0 - (2.0 * clamped));
}

constexpr float kPi = 3.14159265359f;

float RadiansToDegrees(float radians) {
    return radians * (180.0f / kPi);
}

float WrapDegrees(float degrees) {
    float wrapped = std::fmod(degrees, 360.0f);
    if (wrapped < 0.0f) {
        wrapped += 360.0f;
    }
    return wrapped;
}

float AngularDistanceDegrees(float a_deg, float b_deg) {
    float diff = std::fabs(WrapDegrees(a_deg) - WrapDegrees(b_deg));
    if (diff > 180.0f) {
        diff = 360.0f - diff;
    }
    return diff;
}

float MinorRadiusMeters(const Vec3& position, float majorRadius) {
    const float radialXY = std::sqrt((position.x * position.x) + (position.y * position.y));
    return std::sqrt(((radialXY - majorRadius) * (radialXY - majorRadius)) + (position.z * position.z));
}

float ToroidalAngleDegrees(const Vec3& position) {
    return WrapDegrees(RadiansToDegrees(std::atan2(position.y, position.x)));
}

float PoloidalAngleDegrees(const Vec3& position, float majorRadius) {
    const float radialXY = std::sqrt((position.x * position.x) + (position.y * position.y));
    return WrapDegrees(RadiansToDegrees(std::atan2(position.z, radialXY - majorRadius)));
}

bool ParticleMatchesSlice(const ReplayParticle& particle, float majorRadius, const ParticleSliceSettings& slice) {
    if (slice.enableToroidalSlice) {
        const float delta = AngularDistanceDegrees(ToroidalAngleDegrees(particle.position_m), slice.toroidalCenter_deg);
        if (delta > std::max(0.0f, slice.toroidalHalfWidth_deg)) {
            return false;
        }
    }
    if (slice.enablePoloidalSlice) {
        const float delta = AngularDistanceDegrees(PoloidalAngleDegrees(particle.position_m, majorRadius), slice.poloidalCenter_deg);
        if (delta > std::max(0.0f, slice.poloidalHalfWidth_deg)) {
            return false;
        }
    }
    return true;
}

bool ProjectWorldToScreen(
    const Vec3& position,
    const Mat4& viewProjection,
    int viewportWidth,
    int viewportHeight,
    ImVec2* outScreen,
    float* outDepth01) {
    const auto& m = viewProjection.elements;
    const float clipX = (m[0] * position.x) + (m[4] * position.y) + (m[8] * position.z) + m[12];
    const float clipY = (m[1] * position.x) + (m[5] * position.y) + (m[9] * position.z) + m[13];
    const float clipZ = (m[2] * position.x) + (m[6] * position.y) + (m[10] * position.z) + m[14];
    const float clipW = (m[3] * position.x) + (m[7] * position.y) + (m[11] * position.z) + m[15];
    if (clipW <= 1.0e-5f) {
        return false;
    }

    const float invW = 1.0f / clipW;
    const float ndcX = clipX * invW;
    const float ndcY = clipY * invW;
    const float ndcZ = clipZ * invW;
    if (!std::isfinite(ndcX) || !std::isfinite(ndcY) || !std::isfinite(ndcZ)) {
        return false;
    }

    if (outScreen != nullptr) {
        outScreen->x = ((ndcX * 0.5f) + 0.5f) * static_cast<float>(viewportWidth);
        outScreen->y = ((-ndcY * 0.5f) + 0.5f) * static_cast<float>(viewportHeight);
    }
    if (outDepth01 != nullptr) {
        *outDepth01 = Clamp01((ndcZ * 0.5f) + 0.5f);
    }
    return ndcZ >= -1.2f && ndcZ <= 1.2f;
}

void ContinuousPalette(float t, float* r, float* g, float* b) {
    const float clamped = Clamp01(t);

    struct Stop {
        float t;
        float r;
        float g;
        float b;
    };
    constexpr Stop kStops[] = {
        {0.00f, 0.11f, 0.28f, 0.78f},
        {0.35f, 0.12f, 0.80f, 0.86f},
        {0.70f, 0.96f, 0.88f, 0.25f},
        {1.00f, 1.00f, 0.40f, 0.12f},
    };

    for (std::size_t i = 1; i < (sizeof(kStops) / sizeof(kStops[0])); ++i) {
        if (clamped <= kStops[i].t) {
            const float localT = (clamped - kStops[i - 1].t) / (kStops[i].t - kStops[i - 1].t);
            *r = kStops[i - 1].r + ((kStops[i].r - kStops[i - 1].r) * localT);
            *g = kStops[i - 1].g + ((kStops[i].g - kStops[i - 1].g) * localT);
            *b = kStops[i - 1].b + ((kStops[i].b - kStops[i - 1].b) * localT);
            return;
        }
    }

    *r = kStops[3].r;
    *g = kStops[3].g;
    *b = kStops[3].b;
}

void DivergingPalette(float t, float* r, float* g, float* b) {
    const float clamped = Clamp01(t);
    if (clamped <= 0.5f) {
        const float local = clamped / 0.5f;
        *r = 0.12f + ((0.70f - 0.12f) * local);
        *g = 0.32f + ((0.72f - 0.32f) * local);
        *b = 0.84f + ((0.74f - 0.84f) * local);
        return;
    }
    const float local = (clamped - 0.5f) / 0.5f;
    *r = 0.70f + ((0.96f - 0.70f) * local);
    *g = 0.72f + ((0.28f - 0.72f) * local);
    *b = 0.74f + ((0.14f - 0.74f) * local);
}

ImU32 ToImColor(float r, float g, float b, float alpha = 1.0f) {
    return IM_COL32(
        static_cast<int>(Clamp01(r) * 255.0f),
        static_cast<int>(Clamp01(g) * 255.0f),
        static_cast<int>(Clamp01(b) * 255.0f),
        static_cast<int>(Clamp01(alpha) * 255.0f));
}

template <typename T, typename PositionFn>
const T* NearestByScreenDistance(
    const std::vector<T>* points,
    const Mat4& viewProjection,
    int viewportWidth,
    int viewportHeight,
    const ImVec2& mouse,
    float maxDistancePx,
    PositionFn positionFn,
    float* outDistancePx = nullptr) {
    if (points == nullptr) {
        return nullptr;
    }

    const T* nearest = nullptr;
    float nearestDistanceSq = maxDistancePx * maxDistancePx;
    float nearestDepth = std::numeric_limits<float>::infinity();
    for (const T& point : *points) {
        ImVec2 screen;
        float depth01 = 1.0f;
        if (!ProjectWorldToScreen(positionFn(point), viewProjection, viewportWidth, viewportHeight, &screen, &depth01)) {
            continue;
        }
        const float dx = screen.x - mouse.x;
        const float dy = screen.y - mouse.y;
        const float distanceSq = (dx * dx) + (dy * dy);
        if (distanceSq < nearestDistanceSq || (distanceSq <= nearestDistanceSq * 1.05f && depth01 < nearestDepth)) {
            nearest = &point;
            nearestDistanceSq = distanceSq;
            nearestDepth = depth01;
        }
    }

    if (nearest != nullptr && outDistancePx != nullptr) {
        *outDistancePx = std::sqrt(nearestDistanceSq);
    }
    return nearest;
}

struct PowerOnVisualState {
    float startupProgress = 1.0f;
    float fusionGateProgress = 1.0f;
    float ignitionIntensity = 0.0f;
    uint64_t recentFusionEvents = 0;
    std::string label = "steady-state";
};

struct SceneSelectionState {
    std::optional<ReplayParticle> hoveredParticle;
    std::optional<ReplayParticle> pinnedParticle;
    std::optional<ReplayFieldProbe> hoveredProbe;
    std::optional<ReplayFieldProbe> pinnedProbe;
};

PowerOnVisualState ComputePowerOnVisualState(
    double displayedTime_s,
    const ReplayRunConfig& runConfig,
    const ReplaySummaryPoint* summary,
    const ReplaySummaryPoint* windowStartSummary) {
    PowerOnVisualState state;

    if (runConfig.startupRampDuration_s > 0.0) {
        state.startupProgress = static_cast<float>(
            SmoothStep01(displayedTime_s / runConfig.startupRampDuration_s));
    }

    const double fusionGateStart_s =
        std::max(0.0, runConfig.startupRampDuration_s) + std::max(0.0, runConfig.fusionStartDelay_s);
    if (displayedTime_s < fusionGateStart_s) {
        state.fusionGateProgress = 0.0f;
    } else if (runConfig.fusionRampDuration_s > 0.0) {
        state.fusionGateProgress = static_cast<float>(
            SmoothStep01((displayedTime_s - fusionGateStart_s) / runConfig.fusionRampDuration_s));
    } else {
        state.fusionGateProgress = 1.0f;
    }

    if (summary != nullptr) {
        uint64_t windowStartTotal = 0;
        if (windowStartSummary != nullptr &&
            windowStartSummary->fusionEventsTotal <= summary->fusionEventsTotal) {
            windowStartTotal = windowStartSummary->fusionEventsTotal;
        }
        state.recentFusionEvents = summary->fusionEventsTotal - windowStartTotal;

        const float recentFusionSignal = Clamp01(
            static_cast<float>(std::log10(1.0 + static_cast<double>(state.recentFusionEvents)) / 2.0));
        const float cumulativeFusionSignal = Clamp01(
            static_cast<float>(std::log10(1.0 + static_cast<double>(summary->fusionEventsTotal)) / 2.4));
        const float energySignal = Clamp01(static_cast<float>((summary->avgEnergy_keV - 90.0) / 140.0));
        const float fusionEvidence =
            std::max((0.65f * recentFusionSignal) + (0.35f * cumulativeFusionSignal), 0.20f * energySignal);
        state.ignitionIntensity = state.fusionGateProgress * fusionEvidence;
    }

    if (state.startupProgress < 0.995f) {
        state.label = "powering magnetics + beam";
    } else if (state.fusionGateProgress < 0.05f) {
        state.label = "thermalizing core";
    } else if (state.recentFusionEvents == 0) {
        state.label = (state.fusionGateProgress < 0.995f) ? "fusion conditions forming" : "threshold not yet sustained";
    } else if (state.ignitionIntensity < 0.35f) {
        state.label = "fusion onset";
    } else if (state.ignitionIntensity < 0.75f) {
        state.label = "burn building";
    } else {
        state.label = "ignition established";
    }

    return state;
}

struct PhysicsOverlaySettings {
    bool showMagneticArrows = true;
    bool showElectricArrows = false;
    bool showMagneticGuideLines = true;
    bool showParticleForceArrows = false;
    bool showSeededStreamlines = true;
    bool showSliceGuides = true;
    float magneticArrowScale = 0.18f;
    float electricArrowScale = 0.18f;
    float forceArrowScale = 0.28f;
    int maxParticleForceArrows = 48;
};

struct PhysicsOverlayStats {
    double maxMagneticMagnitude_T = 0.0;
    double maxElectricMagnitude_VPerM = 0.0;
    double maxLorentzAcceleration_mPerS2 = 0.0;
    std::size_t fieldProbeCount = 0;
    std::size_t forceArrowCount = 0;
    std::size_t guideLineCount = 0;
    std::size_t seededGuideLineCount = 0;
};

float DegreesToRadians(float degrees) {
    return degrees * (3.14159265359f / 180.0f);
}

Vec3 TorusPoint(float majorRadius, float minorRadius, float phiRadians, float thetaRadians, float rho) {
    const float ringRadius = majorRadius + (minorRadius * rho * std::cos(thetaRadians));
    return Vec3(
        ringRadius * std::cos(phiRadians),
        ringRadius * std::sin(phiRadians),
        minorRadius * rho * std::sin(thetaRadians));
}

void PushColoredVertex(std::vector<float>* out, const Vec3& position, float r, float g, float b) {
    out->push_back(position.x);
    out->push_back(position.y);
    out->push_back(position.z);
    out->push_back(r);
    out->push_back(g);
    out->push_back(b);
}

void PushColoredLine(
    std::vector<float>* out,
    const Vec3& a,
    const Vec3& b,
    float r,
    float g,
    float bl) {
    PushColoredVertex(out, a, r, g, bl);
    PushColoredVertex(out, b, r, g, bl);
}

bool IsInsideTorus(const Vec3& position, float majorRadius, float minorRadius) {
    const float radialXY = std::sqrt((position.x * position.x) + (position.y * position.y));
    const float tubeRadius =
        std::sqrt(((radialXY - majorRadius) * (radialXY - majorRadius)) + (position.z * position.z));
    return tubeRadius <= minorRadius;
}

float MagnitudeScale(double magnitude, double maxMagnitude) {
    if (!std::isfinite(magnitude) || magnitude <= 0.0 || !std::isfinite(maxMagnitude) || maxMagnitude <= 0.0) {
        return 0.0f;
    }
    return Clamp01(static_cast<float>(std::sqrt(magnitude / maxMagnitude)));
}

void PushArrow(
    std::vector<float>* out,
    const Vec3& start,
    const Vec3& direction,
    float length,
    float r,
    float g,
    float b) {
    const Vec3 dir = direction.Normalized();
    if (dir.Magnitude() <= 1.0e-6f || length <= 1.0e-6f) {
        return;
    }

    const Vec3 end = start + (dir * length);
    PushColoredLine(out, start, end, r, g, b);

    Vec3 up(0.0f, 0.0f, 1.0f);
    Vec3 side = Vec3::Cross(dir, up);
    if (side.Magnitude() <= 1.0e-6f) {
        side = Vec3::Cross(dir, Vec3(0.0f, 1.0f, 0.0f));
    }
    side = side.Normalized();

    const float headLength = std::max(0.02f, length * 0.22f);
    const float headWidth = std::max(0.01f, length * 0.12f);
    const Vec3 headBase = end - (dir * headLength);
    PushColoredLine(out, end, headBase + (side * headWidth), r, g, b);
    PushColoredLine(out, end, headBase - (side * headWidth), r, g, b);
}

TokamakConfig EffectiveTokamakConfigForViewer(const ReplayRunConfig& runConfig, double displayedTime_s) {
    TokamakConfig effective = runConfig.hasTokamakConfig ? runConfig.tokamakConfig : TokamakConfig{};
    double startupFraction = 1.0;
    if (runConfig.startupRampDuration_s > 0.0) {
        startupFraction = SmoothStep01(displayedTime_s / runConfig.startupRampDuration_s);
    }
    effective.toroidalCurrent_A *= static_cast<float>(startupFraction);
    effective.plasmaCurrent_A *= static_cast<float>(startupFraction);
    return effective;
}

bool CanDrawAnalyticMagneticGuideLines(const ReplayRunConfig& runConfig) {
    return runConfig.hasTokamakConfig &&
           runConfig.hasPlasmaCurrentProfile &&
           runConfig.plasmaCurrentProfile.kind != PlasmaCurrentProfileKind::CustomTable;
}

void AppendFieldProbeArrows(
    std::vector<float>* out,
    const std::vector<ReplayFieldProbe>* fieldProbes,
    const PhysicsOverlaySettings& settings,
    PhysicsOverlayStats* stats) {
    if (fieldProbes == nullptr || stats == nullptr) {
        return;
    }

    for (const ReplayFieldProbe& probe : *fieldProbes) {
        stats->maxMagneticMagnitude_T = std::max(stats->maxMagneticMagnitude_T, probe.magneticMagnitude_T);
        stats->maxElectricMagnitude_VPerM = std::max(stats->maxElectricMagnitude_VPerM, probe.electricMagnitude_VPerM);
    }
    stats->fieldProbeCount = fieldProbes->size();

    for (const ReplayFieldProbe& probe : *fieldProbes) {
        if (settings.showMagneticArrows) {
            const float normalized = MagnitudeScale(probe.magneticMagnitude_T, stats->maxMagneticMagnitude_T);
            PushArrow(
                out,
                probe.position_m,
                probe.magneticField_T,
                settings.magneticArrowScale * (0.04f + (0.22f * normalized)),
                0.18f,
                0.96f,
                0.84f);
        }
        if (settings.showElectricArrows) {
            const float normalized = MagnitudeScale(probe.electricMagnitude_VPerM, stats->maxElectricMagnitude_VPerM);
            PushArrow(
                out,
                probe.position_m,
                probe.electricField_VPerM,
                settings.electricArrowScale * (0.04f + (0.22f * normalized)),
                1.00f,
                0.44f,
                0.18f);
        }
    }
}

void AppendParticleForceArrows(
    std::vector<float>* out,
    const ReplayFrame& frame,
    const PhysicsOverlaySettings& settings,
    PhysicsOverlayStats* stats) {
    if (stats == nullptr || !settings.showParticleForceArrows || frame.particles.empty()) {
        return;
    }

    for (const ReplayParticle& particle : frame.particles) {
        stats->maxLorentzAcceleration_mPerS2 =
            std::max(stats->maxLorentzAcceleration_mPerS2, particle.lorentzAccelerationMagnitude_mPerS2);
    }

    const std::size_t budget = static_cast<std::size_t>(std::max(1, settings.maxParticleForceArrows));
    const std::size_t stride = std::max<std::size_t>(1, frame.particles.size() / budget);
    for (std::size_t i = 0; i < frame.particles.size(); i += stride) {
        const ReplayParticle& particle = frame.particles[i];
        if (particle.lorentzAccelerationMagnitude_mPerS2 <= 0.0) {
            continue;
        }

        float r = 1.0f;
        float g = 0.92f;
        float b = 0.24f;
        switch (particle.species) {
            case ReplaySpecies::Deuterium:
                r = 0.58f;
                g = 0.86f;
                b = 1.00f;
                break;
            case ReplaySpecies::Tritium:
                r = 1.00f;
                g = 0.72f;
                b = 0.34f;
                break;
            case ReplaySpecies::Helium:
                r = 1.00f;
                g = 0.98f;
                b = 0.50f;
                break;
            case ReplaySpecies::Unknown:
                break;
        }

        const float normalized =
            MagnitudeScale(particle.lorentzAccelerationMagnitude_mPerS2, stats->maxLorentzAcceleration_mPerS2);
        PushArrow(
            out,
            particle.position_m,
            particle.lorentzAcceleration_mPerS2,
            settings.forceArrowScale * (0.03f + (0.18f * normalized)),
            r,
            g,
            b);
        ++stats->forceArrowCount;
    }
}

void AppendMagneticGuideLines(
    std::vector<float>* out,
    const ReplayRunConfig& runConfig,
    double displayedTime_s,
    const PhysicsOverlaySettings& settings,
    PhysicsOverlayStats* stats) {
    if (stats == nullptr || !settings.showMagneticGuideLines || !CanDrawAnalyticMagneticGuideLines(runConfig)) {
        return;
    }

    const TokamakConfig effectiveConfig = EffectiveTokamakConfigForViewer(runConfig, displayedTime_s);
    const PlasmaCurrentProfileConfig& profile = runConfig.plasmaCurrentProfile;
    constexpr float kSeedRho[] = {0.28f, 0.52f, 0.76f};
    constexpr float kSeedTheta_deg[] = {20.0f, 180.0f};
    const float seedPhi = DegreesToRadians(15.0f);

    for (const float rho : kSeedRho) {
        for (const float thetaDeg : kSeedTheta_deg) {
            Vec3 position = TorusPoint(
                effectiveConfig.majorRadius_m,
                effectiveConfig.minorRadius_m,
                seedPhi,
                DegreesToRadians(thetaDeg),
                rho);
            constexpr int kSteps = 120;
            constexpr float kStepDistance_m = 0.055f;
            for (int step = 0; step < kSteps; ++step) {
                const MagneticFieldSample magnetic = EvaluateMagneticFieldSample(effectiveConfig, profile, position);
                const Vec3 direction = magnetic.totalField_T.Normalized();
                if (direction.Magnitude() <= 1.0e-6f) {
                    break;
                }
                const Vec3 next = position + (direction * kStepDistance_m);
                if (!IsInsideTorus(next, effectiveConfig.majorRadius_m, effectiveConfig.minorRadius_m)) {
                    break;
                }
                PushColoredLine(out, position, next, 0.22f, 0.82f, 1.00f);
                position = next;
                ++stats->guideLineCount;
            }
        }
    }
}

void AppendSliceGuides(
    std::vector<float>* out,
    const ReplayRunConfig& runConfig,
    const PhysicsOverlaySettings& settings,
    const ParticleSliceSettings& slice) {
    if (!settings.showSliceGuides || !runConfig.hasTokamakGeometry) {
        return;
    }

    constexpr int kSegments = 72;
    constexpr float kTwoPi = 2.0f * kPi;

    auto appendPoloidalLoop = [&](float phi_deg, float r, float g, float b) {
        const float phi = DegreesToRadians(phi_deg);
        Vec3 previous = TorusPoint(runConfig.majorRadius_m, runConfig.minorRadius_m, phi, 0.0f, 1.0f);
        for (int segment = 1; segment <= kSegments; ++segment) {
            const float theta = kTwoPi * static_cast<float>(segment) / static_cast<float>(kSegments);
            const Vec3 current = TorusPoint(runConfig.majorRadius_m, runConfig.minorRadius_m, phi, theta, 1.0f);
            PushColoredLine(out, previous, current, r, g, b);
            previous = current;
        }
    };

    auto appendToroidalLoop = [&](float theta_deg, float r, float g, float b) {
        const float theta = DegreesToRadians(theta_deg);
        Vec3 previous = TorusPoint(runConfig.majorRadius_m, runConfig.minorRadius_m, 0.0f, theta, 1.0f);
        for (int segment = 1; segment <= kSegments; ++segment) {
            const float phi = kTwoPi * static_cast<float>(segment) / static_cast<float>(kSegments);
            const Vec3 current = TorusPoint(runConfig.majorRadius_m, runConfig.minorRadius_m, phi, theta, 1.0f);
            PushColoredLine(out, previous, current, r, g, b);
            previous = current;
        }
    };

    if (slice.enableToroidalSlice) {
        appendPoloidalLoop(slice.toroidalCenter_deg, 1.00f, 0.82f, 0.32f);
        appendPoloidalLoop(slice.toroidalCenter_deg - slice.toroidalHalfWidth_deg, 0.32f, 0.86f, 1.00f);
        appendPoloidalLoop(slice.toroidalCenter_deg + slice.toroidalHalfWidth_deg, 0.32f, 0.86f, 1.00f);
    }
    if (slice.enablePoloidalSlice) {
        appendToroidalLoop(slice.poloidalCenter_deg, 1.00f, 0.62f, 0.22f);
        appendToroidalLoop(slice.poloidalCenter_deg - slice.poloidalHalfWidth_deg, 0.92f, 0.28f, 0.22f);
        appendToroidalLoop(slice.poloidalCenter_deg + slice.poloidalHalfWidth_deg, 0.92f, 0.28f, 0.22f);
    }
}

void AppendMagneticStreamlineFromSeed(
    std::vector<float>* out,
    const TokamakConfig& config,
    const PlasmaCurrentProfileConfig& profile,
    const Vec3& seed,
    float r,
    float g,
    float b,
    PhysicsOverlayStats* stats) {
    if (out == nullptr || stats == nullptr || !seed.IsFinite()) {
        return;
    }

    auto traceDirection = [&](float sign) {
        Vec3 position = seed;
        constexpr int kSteps = 96;
        constexpr float kStepDistance_m = 0.045f;
        for (int step = 0; step < kSteps; ++step) {
            const MagneticFieldSample magnetic = EvaluateMagneticFieldSample(config, profile, position);
            const Vec3 direction = magnetic.totalField_T.Normalized() * sign;
            if (direction.Magnitude() <= 1.0e-6f) {
                break;
            }
            const Vec3 next = position + (direction * kStepDistance_m);
            if (!IsInsideTorus(next, config.majorRadius_m, config.minorRadius_m)) {
                break;
            }
            PushColoredLine(out, position, next, r, g, b);
            position = next;
            ++stats->seededGuideLineCount;
        }
    };

    traceDirection(1.0f);
    traceDirection(-1.0f);
}

void AppendSeededMagneticStreamlines(
    std::vector<float>* out,
    const ReplayRunConfig& runConfig,
    double displayedTime_s,
    const PhysicsOverlaySettings& settings,
    const SceneSelectionState& selection,
    PhysicsOverlayStats* stats) {
    if (!settings.showSeededStreamlines || stats == nullptr || !CanDrawAnalyticMagneticGuideLines(runConfig)) {
        return;
    }

    std::optional<Vec3> seedPosition;
    if (selection.pinnedProbe.has_value()) {
        seedPosition = selection.pinnedProbe->position_m;
    } else if (selection.pinnedParticle.has_value()) {
        seedPosition = selection.pinnedParticle->position_m;
    } else if (selection.hoveredProbe.has_value()) {
        seedPosition = selection.hoveredProbe->position_m;
    } else if (selection.hoveredParticle.has_value()) {
        seedPosition = selection.hoveredParticle->position_m;
    }
    if (!seedPosition.has_value()) {
        return;
    }

    const TokamakConfig effectiveConfig = EffectiveTokamakConfigForViewer(runConfig, displayedTime_s);
    const Vec3 radial = Vec3(seedPosition->x, seedPosition->y, 0.0f).Normalized();
    const Vec3 vertical(0.0f, 0.0f, 1.0f);
    Vec3 binormal = Vec3::Cross(radial, vertical).Normalized();
    if (binormal.Magnitude() <= 1.0e-6f) {
        binormal = Vec3(0.0f, 1.0f, 0.0f);
    }
    const float offset = std::max(0.015f, effectiveConfig.minorRadius_m * 0.06f);
    constexpr float kColors[][3] = {
        {0.95f, 0.92f, 0.30f},
        {0.56f, 0.88f, 1.00f},
        {1.00f, 0.54f, 0.22f},
    };
    const Vec3 seeds[] = {
        *seedPosition,
        *seedPosition + (vertical * offset),
        *seedPosition + (binormal * offset),
    };
    for (std::size_t i = 0; i < (sizeof(seeds) / sizeof(seeds[0])); ++i) {
        const Vec3 seed = IsInsideTorus(seeds[i], effectiveConfig.majorRadius_m, effectiveConfig.minorRadius_m)
            ? seeds[i]
            : *seedPosition;
        AppendMagneticStreamlineFromSeed(
            out,
            effectiveConfig,
            runConfig.plasmaCurrentProfile,
            seed,
            kColors[i][0],
            kColors[i][1],
            kColors[i][2],
            stats);
    }
}

enum class RadialChartMetric : int {
    Density = 0,
    Temperature = 1,
    FusionRate = 2,
};

enum class ToroidalUnwrapMetric : int {
    MagneticMagnitude = 0,
    ElectricMagnitude = 1,
};

struct AnalyticsUiState {
    ParticleViewMode particleViewMode = ParticleViewMode::Species;
    RadialChartMetric radialMetric = RadialChartMetric::Density;
    int focusStartOrderedIndex = 0;
    int focusEndOrderedIndex = -1;
    int selectedRadialBin = -1;
    ParticleSliceSettings slice;
    bool compareEnabled = false;
    bool diffAgainstCompare = false;
    bool heatmapShowsDelta = false;
    ToroidalUnwrapMetric unwrapMetric = ToroidalUnwrapMetric::MagneticMagnitude;
    int unwrapRhoBucket = 0;
};


struct PlotInteraction {
    int hoveredIndex = -1;
    int clickedIndex = -1;
};

const char* RadialChartMetricName(RadialChartMetric metric) {
    switch (metric) {
        case RadialChartMetric::Density:
            return "Density";
        case RadialChartMetric::Temperature:
            return "Avg Ion Energy";
        case RadialChartMetric::FusionRate:
            return "Fusion Rate";
    }
    return "Unknown";
}

std::vector<float> ToPlotValues(const std::vector<double>& values) {
    std::vector<float> out;
    out.reserve(values.size());
    for (double value : values) {
        out.push_back(std::isfinite(value) ? static_cast<float>(value) : 0.0f);
    }
    return out;
}

void ComputeFiniteRange(const std::vector<double>& values, float* outMin, float* outMax) {
    double minValue = std::numeric_limits<double>::infinity();
    double maxValue = -std::numeric_limits<double>::infinity();
    for (double value : values) {
        if (!std::isfinite(value)) {
            continue;
        }
        minValue = std::min(minValue, value);
        maxValue = std::max(maxValue, value);
    }

    if (!std::isfinite(minValue) || !std::isfinite(maxValue)) {
        *outMin = 0.0f;
        *outMax = 1.0f;
        return;
    }
    if (minValue == maxValue) {
        const double padding = (minValue == 0.0) ? 1.0 : std::fabs(minValue) * 0.1;
        minValue -= padding;
        maxValue += padding;
    }

    *outMin = static_cast<float>(minValue);
    *outMax = static_cast<float>(maxValue);
}

int FindNearestSeriesIndexForStep(const std::vector<int>& steps, int step) {
    if (steps.empty()) {
        return -1;
    }
    const auto it = std::lower_bound(steps.begin(), steps.end(), step);
    if (it == steps.end()) {
        return static_cast<int>(steps.size() - 1);
    }
    if (it == steps.begin()) {
        return 0;
    }
    const int upperIndex = static_cast<int>(std::distance(steps.begin(), it));
    const int lowerIndex = upperIndex - 1;
    return (std::abs(steps[upperIndex] - step) < std::abs(steps[lowerIndex] - step))
        ? upperIndex
        : lowerIndex;
}

PlotInteraction PlotSeriesWithOverlay(
    const char* label,
    const std::vector<double>& values,
    int currentIndex,
    int focusStartIndex,
    int focusEndIndex,
    const char* overlayText,
    const ImVec2& size) {
    PlotInteraction interaction;
    const std::vector<float> plotValues = ToPlotValues(values);
    if (plotValues.empty()) {
        ImGui::Text("%s: no data", label);
        return interaction;
    }
    float minValue = 0.0f;
    float maxValue = 1.0f;
    ComputeFiniteRange(values, &minValue, &maxValue);

    ImGui::PlotLines(label, plotValues.data(), static_cast<int>(plotValues.size()), 0, overlayText, minValue, maxValue, size);

    const ImVec2 rectMin = ImGui::GetItemRectMin();
    const ImVec2 rectMax = ImGui::GetItemRectMax();
    const ImVec2 rectSize(rectMax.x - rectMin.x, rectMax.y - rectMin.y);
    if (rectSize.x <= 1.0f || rectSize.y <= 1.0f) {
        return interaction;
    }

    auto sampleX = [&](int index) -> float {
        if (plotValues.size() <= 1) {
            return rectMin.x + (rectSize.x * 0.5f);
        }
        const float t = static_cast<float>(index) / static_cast<float>(plotValues.size() - 1);
        return rectMin.x + (t * rectSize.x);
    };

    ImDrawList* drawList = ImGui::GetWindowDrawList();
    if (focusStartIndex >= 0 && focusEndIndex >= focusStartIndex &&
        focusStartIndex < static_cast<int>(plotValues.size())) {
        const int clampedEnd = std::min(focusEndIndex, static_cast<int>(plotValues.size() - 1));
        drawList->AddRectFilled(
            ImVec2(sampleX(focusStartIndex), rectMin.y),
            ImVec2(sampleX(clampedEnd), rectMax.y),
            IM_COL32(44, 120, 140, 32));
    }

    if (currentIndex >= 0 && currentIndex < static_cast<int>(plotValues.size())) {
        const float x = sampleX(currentIndex);
        drawList->AddLine(ImVec2(x, rectMin.y), ImVec2(x, rectMax.y), IM_COL32(255, 220, 120, 220), 2.0f);
    }

    if (ImGui::IsItemHovered()) {
        const ImVec2 mouse = ImGui::GetIO().MousePos;
        const float localX = Clamp01((mouse.x - rectMin.x) / std::max(rectSize.x, 1.0f));
        const int hovered = static_cast<int>(std::round(localX * static_cast<float>(plotValues.size() - 1)));
        interaction.hoveredIndex = std::max(0, std::min(hovered, static_cast<int>(plotValues.size() - 1)));
        if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
            interaction.clickedIndex = interaction.hoveredIndex;
        }
    }

    return interaction;
}

PlotInteraction PlotHistogramWithOverlay(
    const char* label,
    const std::vector<double>& values,
    int selectedIndex,
    const char* overlayText,
    const ImVec2& size) {
    PlotInteraction interaction;
    const std::vector<float> plotValues = ToPlotValues(values);
    if (plotValues.empty()) {
        ImGui::Text("%s: no data", label);
        return interaction;
    }
    float minValue = 0.0f;
    float maxValue = 1.0f;
    ComputeFiniteRange(values, &minValue, &maxValue);
    minValue = std::min(0.0f, minValue);

    ImGui::PlotHistogram(label, plotValues.data(), static_cast<int>(plotValues.size()), 0, overlayText, minValue, maxValue, size);

    const ImVec2 rectMin = ImGui::GetItemRectMin();
    const ImVec2 rectMax = ImGui::GetItemRectMax();
    const ImVec2 rectSize(rectMax.x - rectMin.x, rectMax.y - rectMin.y);
    if (rectSize.x <= 1.0f || rectSize.y <= 1.0f) {
        return interaction;
    }

    const float barWidth = rectSize.x / static_cast<float>(plotValues.size());
    ImDrawList* drawList = ImGui::GetWindowDrawList();
    if (selectedIndex >= 0 && selectedIndex < static_cast<int>(plotValues.size())) {
        const float x0 = rectMin.x + (barWidth * static_cast<float>(selectedIndex));
        const float x1 = x0 + barWidth;
        drawList->AddRect(ImVec2(x0, rectMin.y), ImVec2(x1, rectMax.y), IM_COL32(255, 220, 120, 220), 0.0f, 0, 2.0f);
    }

    if (ImGui::IsItemHovered()) {
        const ImVec2 mouse = ImGui::GetIO().MousePos;
        const int hovered = static_cast<int>((mouse.x - rectMin.x) / std::max(barWidth, 1.0f));
        interaction.hoveredIndex = std::max(0, std::min(hovered, static_cast<int>(plotValues.size() - 1)));
        if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
            interaction.clickedIndex = interaction.hoveredIndex;
        }
    }

    return interaction;
}

template <typename T, typename ValueFn>
void BuildTimelineSeries(const std::vector<T>& rows, ValueFn valueFn, std::vector<int>* outSteps, std::vector<double>* outValues) {
    outSteps->clear();
    outValues->clear();
    outSteps->reserve(rows.size());
    outValues->reserve(rows.size());
    for (const T& row : rows) {
        outSteps->push_back(row.step);
        outValues->push_back(valueFn(row));
    }
}

double RadialMetricValue(const ReplayRadialProfileBin& bin, RadialChartMetric metric) {
    switch (metric) {
        case RadialChartMetric::Density:
            return bin.density_m3;
        case RadialChartMetric::Temperature:
            return bin.avgIonEnergy_keV;
        case RadialChartMetric::FusionRate:
            return bin.fusionRatePlaceholder ? 0.0 : bin.fusionRate_m3_s;
    }
    return 0.0;
}

const char* RadialMetricUnits(RadialChartMetric metric) {
    switch (metric) {
        case RadialChartMetric::Density:
            return "m^-3";
        case RadialChartMetric::Temperature:
            return "keV";
        case RadialChartMetric::FusionRate:
            return "m^-3 s^-1";
    }
    return "";
}

const char* ToroidalUnwrapMetricName(ToroidalUnwrapMetric metric) {
    switch (metric) {
        case ToroidalUnwrapMetric::MagneticMagnitude:
            return "|B|";
        case ToroidalUnwrapMetric::ElectricMagnitude:
            return "|E|";
    }
    return "Unknown";
}

double ToroidalUnwrapMetricValue(const ReplayFieldProbe& probe, ToroidalUnwrapMetric metric) {
    switch (metric) {
        case ToroidalUnwrapMetric::MagneticMagnitude:
            return probe.magneticMagnitude_T;
        case ToroidalUnwrapMetric::ElectricMagnitude:
            return probe.electricMagnitude_VPerM;
    }
    return 0.0;
}

ParticleViewContext BuildParticleViewContext(
    ParticleViewMode mode,
    const ReplayFrame& frame,
    const ParticleSliceSettings& slice,
    const std::vector<ReplayRadialProfileBin>* radialProfile,
    const std::vector<ReplayRadialProfileBin>* compareRadialProfile,
    const ReplayWallInteractionPoint* wallPoint,
    const ReplayWallInteractionPoint* compareWallPoint,
    float majorRadius_m,
    float minorRadius_m,
    bool diffAgainstCompare) {
    ParticleViewContext context;
    context.mode = mode;
    context.radialProfile = radialProfile;
    context.compareRadialProfile = compareRadialProfile;
    context.wallInteraction = wallPoint;
    context.compareWallInteraction = compareWallPoint;
    context.majorRadius_m = majorRadius_m;
    context.minorRadius_m = minorRadius_m;
    context.diffAgainstCompare = diffAgainstCompare &&
        (mode == ParticleViewMode::FusionRate || mode == ParticleViewMode::WallLossRisk);
    context.slice = slice;
    context.range = ComputeParticleViewRange(frame, context);
    return context;
}

double CompareRadialMetricValue(
    const std::vector<ReplayRadialProfileBin>* compareProfile,
    int binIndex,
    RadialChartMetric metric) {
    if (compareProfile == nullptr) {
        return 0.0;
    }
    for (const ReplayRadialProfileBin& bin : *compareProfile) {
        if (bin.binIndex == binIndex) {
            return RadialMetricValue(bin, metric);
        }
    }
    return 0.0;
}

void ColorForValue(double value, double minValue, double maxValue, bool diverging, float* r, float* g, float* b) {
    if (!std::isfinite(value)) {
        *r = 0.30f;
        *g = 0.32f;
        *b = 0.35f;
        return;
    }
    if (diverging) {
        const double amplitude = std::max(std::fabs(minValue), std::fabs(maxValue));
        const float t = (amplitude > 1.0e-30)
            ? Clamp01(static_cast<float>((value / amplitude) * 0.5 + 0.5))
            : 0.5f;
        DivergingPalette(t, r, g, b);
        return;
    }
    const double span = maxValue - minValue;
    const float t = (span > 1.0e-30)
        ? Clamp01(static_cast<float>((value - minValue) / span))
        : 1.0f;
    ContinuousPalette(t, r, g, b);
}

struct HeatmapInteraction {
    int hoveredBin = -1;
    int clickedBin = -1;
};

HeatmapInteraction DrawRadialCrossSectionHeatmap(
    const std::vector<ReplayRadialProfileBin>& radialProfile,
    const std::vector<ReplayRadialProfileBin>* compareProfile,
    RadialChartMetric metric,
    bool showDelta,
    const ParticleSliceSettings& slice,
    int selectedBin) {
    HeatmapInteraction interaction;
    ImGui::Text("Poloidal cross-section heatmap");
    const ImVec2 size(248.0f, 248.0f);
    const ImVec2 topLeft = ImGui::GetCursorScreenPos();
    ImGui::InvisibleButton("##poloidal_heatmap", size);
    const ImVec2 rectMin = topLeft;
    const ImVec2 rectMax(topLeft.x + size.x, topLeft.y + size.y);
    ImDrawList* drawList = ImGui::GetWindowDrawList();
    const ImVec2 center((rectMin.x + rectMax.x) * 0.5f, (rectMin.y + rectMax.y) * 0.5f);
    const float maxRadiusPx = (std::min(size.x, size.y) * 0.5f) - 10.0f;

    std::vector<double> values;
    values.reserve(radialProfile.size());
    double minValue = std::numeric_limits<double>::infinity();
    double maxValue = -std::numeric_limits<double>::infinity();
    double maxOuterRadius_m = 0.0;
    for (const ReplayRadialProfileBin& bin : radialProfile) {
        const double currentValue = RadialMetricValue(bin, metric);
        const double compareValue = showDelta ? CompareRadialMetricValue(compareProfile, bin.binIndex, metric) : 0.0;
        const double value = showDelta ? (currentValue - compareValue) : currentValue;
        values.push_back(value);
        minValue = std::min(minValue, value);
        maxValue = std::max(maxValue, value);
        maxOuterRadius_m = std::max(maxOuterRadius_m, bin.rOuter_m);
    }
    if (!std::isfinite(minValue) || !std::isfinite(maxValue)) {
        minValue = 0.0;
        maxValue = 1.0;
    }

    drawList->AddRectFilled(rectMin, rectMax, IM_COL32(14, 20, 28, 220), 8.0f);
    for (int i = static_cast<int>(radialProfile.size()) - 1; i >= 0; --i) {
        const ReplayRadialProfileBin& bin = radialProfile[static_cast<std::size_t>(i)];
        const float outerRadius = maxRadiusPx * static_cast<float>(bin.rOuter_m / std::max(1.0e-6, maxOuterRadius_m));
        float r = 0.0f;
        float g = 0.0f;
        float b = 0.0f;
        ColorForValue(values[static_cast<std::size_t>(i)], minValue, maxValue, showDelta, &r, &g, &b);
        drawList->AddCircleFilled(center, outerRadius, ToImColor(r, g, b, 0.95f), 64);
        drawList->AddCircle(center, outerRadius, IM_COL32(18, 26, 34, 220), 64, 1.2f);
    }
    drawList->AddCircle(center, maxRadiusPx, IM_COL32(180, 210, 220, 160), 80, 1.5f);

    if (slice.enablePoloidalSlice) {
        const float centerRad = DegreesToRadians(slice.poloidalCenter_deg);
        const float lowRad = DegreesToRadians(slice.poloidalCenter_deg - slice.poloidalHalfWidth_deg);
        const float highRad = DegreesToRadians(slice.poloidalCenter_deg + slice.poloidalHalfWidth_deg);
        const ImU32 centerColor = IM_COL32(255, 220, 120, 220);
        const ImU32 edgeColor = IM_COL32(120, 220, 255, 180);
        drawList->AddLine(
            center,
            ImVec2(center.x + (std::cos(centerRad) * maxRadiusPx), center.y - (std::sin(centerRad) * maxRadiusPx)),
            centerColor,
            2.0f);
        drawList->AddLine(
            center,
            ImVec2(center.x + (std::cos(lowRad) * maxRadiusPx), center.y - (std::sin(lowRad) * maxRadiusPx)),
            edgeColor,
            1.2f);
        drawList->AddLine(
            center,
            ImVec2(center.x + (std::cos(highRad) * maxRadiusPx), center.y - (std::sin(highRad) * maxRadiusPx)),
            edgeColor,
            1.2f);
    }

    if (selectedBin >= 0 && selectedBin < static_cast<int>(radialProfile.size())) {
        const ReplayRadialProfileBin& bin = radialProfile[static_cast<std::size_t>(selectedBin)];
        const float outerRadius = maxRadiusPx * static_cast<float>(bin.rOuter_m / std::max(1.0e-6, maxOuterRadius_m));
        const float innerRadius = maxRadiusPx * static_cast<float>(bin.rInner_m / std::max(1.0e-6, maxOuterRadius_m));
        drawList->AddCircle(center, outerRadius, IM_COL32(255, 220, 120, 255), 80, 2.0f);
        if (innerRadius > 1.0f) {
            drawList->AddCircle(center, innerRadius, IM_COL32(255, 220, 120, 180), 80, 1.4f);
        }
    }

    if (ImGui::IsItemHovered()) {
        const ImVec2 mouse = ImGui::GetIO().MousePos;
        const float dx = mouse.x - center.x;
        const float dy = mouse.y - center.y;
        const float radiusPx = std::sqrt((dx * dx) + (dy * dy));
        const double radius_m = (radiusPx / std::max(maxRadiusPx, 1.0f)) * maxOuterRadius_m;
        for (std::size_t i = 0; i < radialProfile.size(); ++i) {
            const ReplayRadialProfileBin& bin = radialProfile[i];
            if (radius_m >= bin.rInner_m && radius_m <= bin.rOuter_m) {
                interaction.hoveredBin = static_cast<int>(i);
                if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                    interaction.clickedBin = interaction.hoveredBin;
                }
                break;
            }
        }
    }

    ImGui::Text(
        "%s%s",
        RadialChartMetricName(metric),
        showDelta ? " delta vs compare" : "");
    ImGui::Text("%s range: %.3e -> %.3e %s", showDelta ? "Delta" : "Value", minValue, maxValue, RadialMetricUnits(metric));
    return interaction;
}

struct ToroidalUnwrapInteraction {
    int hoveredProbeIndex = -1;
    int clickedProbeIndex = -1;
};

ToroidalUnwrapInteraction DrawToroidalUnwrap(
    const std::vector<ReplayFieldProbe>& fieldProbes,
    const std::vector<double>& rhoBuckets,
    int selectedRhoBucket,
    ToroidalUnwrapMetric metric,
    const std::optional<ReplayFieldProbe>& pinnedProbe) {
    ToroidalUnwrapInteraction interaction;
    ImGui::Text("Toroidal unwrap");
    const ImVec2 size(360.0f, 220.0f);
    const ImVec2 topLeft = ImGui::GetCursorScreenPos();
    ImGui::InvisibleButton("##toroidal_unwrap", size);
    const ImVec2 rectMin = topLeft;
    const ImVec2 rectMax(topLeft.x + size.x, topLeft.y + size.y);
    ImDrawList* drawList = ImGui::GetWindowDrawList();
    drawList->AddRectFilled(rectMin, rectMax, IM_COL32(14, 20, 28, 220), 8.0f);

    if (fieldProbes.empty() || rhoBuckets.empty()) {
        return interaction;
    }

    const int rhoIndex = std::max(0, std::min(selectedRhoBucket, static_cast<int>(rhoBuckets.size() - 1)));
    const double rhoTarget = rhoBuckets[static_cast<std::size_t>(rhoIndex)];
    std::vector<int> matchingIndices;
    matchingIndices.reserve(fieldProbes.size());
    double minValue = std::numeric_limits<double>::infinity();
    double maxValue = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < fieldProbes.size(); ++i) {
        const ReplayFieldProbe& probe = fieldProbes[i];
        if (std::fabs(probe.rho - rhoTarget) > 1.0e-4) {
            continue;
        }
        matchingIndices.push_back(static_cast<int>(i));
        const double value = ToroidalUnwrapMetricValue(probe, metric);
        minValue = std::min(minValue, value);
        maxValue = std::max(maxValue, value);
    }

    if (matchingIndices.empty()) {
        ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.25f, 1.0f), "No probes exist on the selected rho shell.");
        return interaction;
    }
    if (!std::isfinite(minValue) || !std::isfinite(maxValue)) {
        minValue = 0.0;
        maxValue = 1.0;
    }

    for (int grid = 0; grid <= 6; ++grid) {
        const float x = rectMin.x + ((rectMax.x - rectMin.x) * static_cast<float>(grid) / 6.0f);
        drawList->AddLine(ImVec2(x, rectMin.y), ImVec2(x, rectMax.y), IM_COL32(44, 60, 78, 120), 1.0f);
        const float y = rectMin.y + ((rectMax.y - rectMin.y) * static_cast<float>(grid) / 6.0f);
        drawList->AddLine(ImVec2(rectMin.x, y), ImVec2(rectMax.x, y), IM_COL32(44, 60, 78, 120), 1.0f);
    }

    if (ImGui::IsItemHovered()) {
        const ImVec2 mouse = ImGui::GetIO().MousePos;
        float nearestDistanceSq = 18.0f * 18.0f;
        for (int index : matchingIndices) {
            const ReplayFieldProbe& probe = fieldProbes[static_cast<std::size_t>(index)];
            const float x = rectMin.x + (static_cast<float>(probe.phi_deg) / 360.0f) * (rectMax.x - rectMin.x);
            const float y = rectMax.y - (static_cast<float>(probe.theta_deg) / 360.0f) * (rectMax.y - rectMin.y);
            const float dx = mouse.x - x;
            const float dy = mouse.y - y;
            const float distanceSq = (dx * dx) + (dy * dy);
            if (distanceSq < nearestDistanceSq) {
                nearestDistanceSq = distanceSq;
                interaction.hoveredProbeIndex = index;
            }
        }
        if (interaction.hoveredProbeIndex >= 0 && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
            interaction.clickedProbeIndex = interaction.hoveredProbeIndex;
        }
    }

    for (int index : matchingIndices) {
        const ReplayFieldProbe& probe = fieldProbes[static_cast<std::size_t>(index)];
        const float x = rectMin.x + (static_cast<float>(probe.phi_deg) / 360.0f) * (rectMax.x - rectMin.x);
        const float y = rectMax.y - (static_cast<float>(probe.theta_deg) / 360.0f) * (rectMax.y - rectMin.y);
        float r = 0.0f;
        float g = 0.0f;
        float b = 0.0f;
        ColorForValue(ToroidalUnwrapMetricValue(probe, metric), minValue, maxValue, false, &r, &g, &b);
        const bool pinned = pinnedProbe.has_value() && pinnedProbe->probeIndex == probe.probeIndex;
        const bool hovered = interaction.hoveredProbeIndex == index;
        const float halfSize = hovered ? 9.0f : 7.0f;
        drawList->AddRectFilled(
            ImVec2(x - halfSize, y - halfSize),
            ImVec2(x + halfSize, y + halfSize),
            ToImColor(r, g, b, 0.95f),
            2.0f);
        drawList->AddRect(
            ImVec2(x - halfSize, y - halfSize),
            ImVec2(x + halfSize, y + halfSize),
            pinned ? IM_COL32(255, 220, 120, 255) : IM_COL32(18, 26, 34, 220),
            2.0f,
            0,
            pinned ? 2.0f : 1.0f);
    }

    ImGui::Text("%s on rho = %.2f", ToroidalUnwrapMetricName(metric), rhoTarget);
    ImGui::Text("Phi 0 -> 360 deg | Theta 0 -> 360 deg");
    return interaction;
}

bool DrawAngleDial(const char* label, float* center_deg, float* halfWidth_deg, float maxHalfWidth_deg) {
    if (center_deg == nullptr || halfWidth_deg == nullptr) {
        return false;
    }

    bool changed = false;
    ImGui::PushID(label);
    ImGui::TextUnformatted(label);
    const ImVec2 size(160.0f, 160.0f);
    const ImVec2 topLeft = ImGui::GetCursorScreenPos();
    ImGui::InvisibleButton("dial", size);
    const ImVec2 rectMin = topLeft;
    const ImVec2 rectMax(topLeft.x + size.x, topLeft.y + size.y);
    const ImVec2 center((rectMin.x + rectMax.x) * 0.5f, (rectMin.y + rectMax.y) * 0.5f);
    const float radius = (std::min(size.x, size.y) * 0.5f) - 16.0f;
    const bool hovered = ImGui::IsItemHovered();
    const bool active = ImGui::IsItemActive();
    ImDrawList* drawList = ImGui::GetWindowDrawList();

    if ((hovered || active) && ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
        const ImVec2 mouse = ImGui::GetIO().MousePos;
        const float angle = WrapDegrees(RadiansToDegrees(std::atan2(center.y - mouse.y, mouse.x - center.x)));
        if (ImGui::GetIO().KeyShift) {
            *halfWidth_deg = std::max(4.0f, std::min(maxHalfWidth_deg, AngularDistanceDegrees(angle, *center_deg)));
        } else {
            *center_deg = angle;
        }
        changed = true;
    }

    drawList->AddCircleFilled(center, radius + 10.0f, IM_COL32(16, 22, 30, 220), 64);
    drawList->AddCircle(center, radius, IM_COL32(90, 120, 138, 200), 64, 1.5f);

    const float centerRad = DegreesToRadians(*center_deg);
    const float lowRad = DegreesToRadians(*center_deg - *halfWidth_deg);
    const float highRad = DegreesToRadians(*center_deg + *halfWidth_deg);
    const ImVec2 handle(center.x + (std::cos(centerRad) * radius), center.y - (std::sin(centerRad) * radius));
    drawList->AddLine(center, handle, IM_COL32(255, 220, 120, 235), 2.4f);
    drawList->AddCircleFilled(handle, 6.0f, IM_COL32(255, 220, 120, 255), 24);

    for (int i = 0; i < 48; ++i) {
        const float a0 = lowRad + ((highRad - lowRad) * static_cast<float>(i) / 48.0f);
        const float a1 = lowRad + ((highRad - lowRad) * static_cast<float>(i + 1) / 48.0f);
        drawList->AddLine(
            ImVec2(center.x + (std::cos(a0) * (radius - 8.0f)), center.y - (std::sin(a0) * (radius - 8.0f))),
            ImVec2(center.x + (std::cos(a1) * (radius - 8.0f)), center.y - (std::sin(a1) * (radius - 8.0f))),
            IM_COL32(100, 220, 255, 180),
            3.0f);
    }
    drawList->AddLine(
        center,
        ImVec2(center.x + (std::cos(lowRad) * radius), center.y - (std::sin(lowRad) * radius)),
        IM_COL32(120, 220, 255, 180),
        1.2f);
    drawList->AddLine(
        center,
        ImVec2(center.x + (std::cos(highRad) * radius), center.y - (std::sin(highRad) * radius)),
        IM_COL32(120, 220, 255, 180),
        1.2f);

    ImGui::Text("Center %.1f deg | +/- %.1f deg", *center_deg, *halfWidth_deg);
    ImGui::TextDisabled("Drag to rotate. Hold Shift while dragging to widen.");
    ImGui::PopID();
    return changed;
}

void DrawParticleViewLegend(
    const ParticleViewContext& context,
    const ReplayFrame& currentFrame) {
    const ParticleViewMode mode = context.mode;
    const ParticleViewRange& range = context.range;
    ImGui::Text("Particle view: %s", ParticleViewModeName(mode));

    if (mode == ParticleViewMode::Species) {
        struct SpeciesEntry {
            ReplaySpecies species;
            const char* label;
        };
        constexpr SpeciesEntry kEntries[] = {
            {ReplaySpecies::Deuterium, "Deuterium"},
            {ReplaySpecies::Tritium, "Tritium"},
            {ReplaySpecies::Helium, "Helium"},
            {ReplaySpecies::Unknown, "Unknown"},
        };
        for (const SpeciesEntry& entry : kEntries) {
            float r = 0.0f;
            float g = 0.0f;
            float b = 0.0f;
            ReplayParticle particle;
            particle.species = entry.species;
            ComputeParticleViewColor(particle, context, &r, &g, &b);
            ImGui::ColorButton(entry.label, ImVec4(r, g, b, 1.0f), ImGuiColorEditFlags_NoTooltip, ImVec2(18.0f, 18.0f));
            ImGui::SameLine();
            ImGui::TextUnformatted(entry.label);
        }
        return;
    }

    ImDrawList* drawList = ImGui::GetWindowDrawList();
    const ImVec2 start = ImGui::GetCursorScreenPos();
    constexpr float kLegendWidth = 220.0f;
    constexpr float kLegendHeight = 14.0f;
    for (int i = 0; i < 120; ++i) {
        const float t0 = static_cast<float>(i) / 120.0f;
        const float t1 = static_cast<float>(i + 1) / 120.0f;
        ReplayParticle particle;
        float r = 0.0f;
        float g = 0.0f;
        float b = 0.0f;
        particle.kineticEnergy_keV = range.minValue + ((range.maxValue - range.minValue) * t0);
        particle.pitchAngle_deg = range.minValue + ((range.maxValue - range.minValue) * t0);
        particle.magneticMagnitude_T = range.minValue + ((range.maxValue - range.minValue) * t0);
        particle.electricMagnitude_VPerM = range.minValue + ((range.maxValue - range.minValue) * t0);
        particle.lorentzAccelerationMagnitude_mPerS2 = range.minValue + ((range.maxValue - range.minValue) * t0);
        particle.position_m = Vec3(
            context.majorRadius_m + (context.minorRadius_m * t0),
            0.0f,
            0.0f);
        ComputeParticleViewColor(particle, context, &r, &g, &b);
        drawList->AddRectFilled(
            ImVec2(start.x + (kLegendWidth * t0), start.y),
            ImVec2(start.x + (kLegendWidth * t1), start.y + kLegendHeight),
            IM_COL32(
                static_cast<int>(r * 255.0f),
                static_cast<int>(g * 255.0f),
                static_cast<int>(b * 255.0f),
                255));
    }
    ImGui::Dummy(ImVec2(kLegendWidth, kLegendHeight + 4.0f));

    char minText[64];
    char maxText[64];
    std::snprintf(minText, sizeof(minText), "%.3e", range.minValue);
    std::snprintf(maxText, sizeof(maxText), "%.3e", range.maxValue);
    ImGui::Text("%s  ->  %s", minText, maxText);

    if (context.diffAgainstCompare &&
        (mode == ParticleViewMode::FusionRate || mode == ParticleViewMode::WallLossRisk)) {
        ImGui::TextUnformatted("Blue = lower than compare run, orange = higher than compare run.");
    } else if (mode == ParticleViewMode::PitchAngle) {
        ImGui::TextUnformatted("0 deg = field-aligned, 180 deg = counter-aligned.");
    } else if (!range.valid && currentFrame.particles.empty()) {
        ImGui::TextUnformatted("Current frame has no visible particles.");
    }
}

}  // namespace

ViewerApp::ViewerApp(ViewerCliOptions options)
    : options_(std::move(options)) {}

bool ParseViewerCliArgs(int argc, char** argv, ViewerCliOptions* outOptions, std::string* errorOut) {
    if (outOptions == nullptr) {
        if (errorOut != nullptr) {
            *errorOut = "Internal error: outOptions is null";
        }
        return false;
    }

    ViewerCliOptions parsed;

    for (int i = 1; i < argc; ++i) {
        const std::string_view arg(argv[i]);
        auto needValue = [&](const char* optionName) -> const char* {
            if (i + 1 >= argc) {
                if (errorOut != nullptr) {
                    *errorOut = std::string("Missing value for ") + optionName;
                }
                return nullptr;
            }
            return argv[++i];
        };

        if (arg == "--manifest") {
            const char* value = needValue("--manifest");
            if (value == nullptr) {
                return false;
            }
            parsed.manifestPath = value;
            continue;
        }

        if (arg == "--run-dir") {
            const char* value = needValue("--run-dir");
            if (value == nullptr) {
                return false;
            }
            parsed.runDirectory = value;
            continue;
        }

        if (arg == "--compare-manifest") {
            const char* value = needValue("--compare-manifest");
            if (value == nullptr) {
                return false;
            }
            parsed.compareManifestPath = value;
            continue;
        }

        if (arg == "--compare-run-dir") {
            const char* value = needValue("--compare-run-dir");
            if (value == nullptr) {
                return false;
            }
            parsed.compareRunDirectory = value;
            continue;
        }

        if (arg == "--point-size") {
            const char* value = needValue("--point-size");
            if (value == nullptr) {
                return false;
            }
            float parsedValue = 0.0f;
            if (!ParseFloatValue(value, &parsedValue) || parsedValue <= 0.0f) {
                if (errorOut != nullptr) {
                    *errorOut = std::string("Invalid --point-size: ") + value;
                }
                return false;
            }
            parsed.pointSizePixels = parsedValue;
            continue;
        }

        if (arg == "--playback-rate") {
            const char* value = needValue("--playback-rate");
            if (value == nullptr) {
                return false;
            }
            float parsedValue = 0.0f;
            if (!ParseFloatValue(value, &parsedValue) || parsedValue <= 0.0f) {
                if (errorOut != nullptr) {
                    *errorOut = std::string("Invalid --playback-rate: ") + value;
                }
                return false;
            }
            parsed.playbackRate = parsedValue;
            continue;
        }

        if (arg == "--start-step") {
            const char* value = needValue("--start-step");
            if (value == nullptr) {
                return false;
            }
            int parsedValue = 0;
            if (!ParseIntValue(value, &parsedValue) || parsedValue < 0) {
                if (errorOut != nullptr) {
                    *errorOut = std::string("Invalid --start-step: ") + value;
                }
                return false;
            }
            parsed.startStep = parsedValue;
            continue;
        }

        if (arg == "--help" || arg == "-h") {
            if (errorOut != nullptr) {
                *errorOut = "help";
            }
            return false;
        }

        if (errorOut != nullptr) {
            *errorOut = std::string("Unknown option: ") + std::string(arg);
        }
        return false;
    }

    const bool hasManifest = !parsed.manifestPath.empty();
    const bool hasRunDir = !parsed.runDirectory.empty();
    if (hasManifest == hasRunDir) {
        if (errorOut != nullptr) {
            *errorOut = "Specify exactly one of --manifest <path> or --run-dir <path>";
        }
        return false;
    }
    const bool hasCompareManifest = !parsed.compareManifestPath.empty();
    const bool hasCompareRunDir = !parsed.compareRunDirectory.empty();
    if (hasCompareManifest && hasCompareRunDir) {
        if (errorOut != nullptr) {
            *errorOut = "Specify at most one of --compare-manifest <path> or --compare-run-dir <path>";
        }
        return false;
    }

    *outOptions = parsed;
    return true;
}

void PrintViewerUsage(const char* argv0) {
    std::cout << "Usage: " << argv0 << " [options]\n"
              << "  --manifest <path/to/manifest_v2.json>\n"
              << "  --run-dir <path/to/run_directory>\n"
              << "  --compare-manifest <path/to/manifest_v2.json>\n"
              << "  --compare-run-dir <path/to/run_directory>\n"
              << "  --point-size <float>\n"
              << "  --playback-rate <float>\n"
              << "  --start-step <int>\n"
              << "  --help\n";
}

int ViewerApp::Run() {
    const bool opened = !options_.manifestPath.empty()
        ? loader_.OpenFromManifest(options_.manifestPath)
        : loader_.OpenFromRunDirectory(options_.runDirectory);

    if (!opened) {
        std::cerr << "Failed to open replay: " << loader_.LastError() << "\n";
        return 1;
    }

    if (!loader_.HasData()) {
        std::cerr << "No replay frames are available in manifest.\n";
        return 1;
    }

    const auto& orderedSteps = loader_.OrderedSteps();
    std::size_t currentFrameIndex = 0;
    if (options_.startStep >= 0) {
        const auto it = std::lower_bound(orderedSteps.begin(), orderedSteps.end(), options_.startStep);
        if (it == orderedSteps.end() || *it != options_.startStep) {
            std::cerr << "Requested --start-step not found: " << options_.startStep << "\n";
            return 1;
        }
        currentFrameIndex = static_cast<std::size_t>(std::distance(orderedSteps.begin(), it));
    }

    ReplayFrame currentFrame;
    if (!loader_.LoadFrameByOrderedIndex(currentFrameIndex, &currentFrame)) {
        std::cerr << "Failed to load initial frame: " << loader_.LastError() << "\n";
        return 1;
    }

    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        return 1;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE);
#endif

    GLFWwindow* window = glfwCreateWindow(1600, 900, "Tokamak Viewer (Replay)", nullptr, nullptr);
    if (window == nullptr) {
        glfwTerminate();
        std::cerr << "Failed to create GLFW window\n";
        return 1;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    if (!gladLoadGLLoader(reinterpret_cast<GLADloadproc>(glfwGetProcAddress))) {
        std::cerr << "Failed to initialize GLAD\n";
        glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");

    WindowInputState inputState;
    glfwSetWindowUserPointer(window, &inputState);
    glfwSetScrollCallback(window, ScrollCallback);

    const ReplayRunConfig& runConfig = loader_.RunConfig();
    const float majorRadius_m = runConfig.hasTokamakGeometry ? runConfig.majorRadius_m : 2.0f;
    const float minorRadius_m = runConfig.hasTokamakGeometry ? runConfig.minorRadius_m : 0.5f;

    ReplayLoader compareLoader;
    const bool hasCompareSource = !options_.compareManifestPath.empty() || !options_.compareRunDirectory.empty();
    bool hasCompareReplay = false;
    if (hasCompareSource) {
        hasCompareReplay = !options_.compareManifestPath.empty()
            ? compareLoader.OpenFromManifest(options_.compareManifestPath)
            : compareLoader.OpenFromRunDirectory(options_.compareRunDirectory);
        if (!hasCompareReplay) {
            std::cerr << "Failed to open compare replay: " << compareLoader.LastError() << "\n";
            ImGui_ImplOpenGL3_Shutdown();
            ImGui_ImplGlfw_Shutdown();
            ImGui::DestroyContext();
            glfwDestroyWindow(window);
            glfwTerminate();
            return 1;
        }
        if (!compareLoader.HasData()) {
            std::cerr << "Compare replay does not contain any frames.\n";
            ImGui_ImplOpenGL3_Shutdown();
            ImGui_ImplGlfw_Shutdown();
            ImGui::DestroyContext();
            glfwDestroyWindow(window);
            glfwTerminate();
            return 1;
        }
    }

    std::size_t initialParticleCapacity = std::max<std::size_t>(
        std::max<std::size_t>(currentFrame.particles.size(), ClampToSizeT(currentFrame.sampledParticles)),
        1);
    if (runConfig.hasMaxParticlesPerSnapshot) {
        initialParticleCapacity = std::max(initialParticleCapacity, ClampToSizeT(runConfig.maxParticlesPerSnapshot));
    }
    if (orderedSteps.size() <= 256 && initialParticleCapacity <= 12000) {
        loader_.SetCacheCapacity(orderedSteps.size());
        for (std::size_t i = 0; i < orderedSteps.size(); ++i) {
            if (!loader_.PrefetchFrameByOrderedIndex(i)) {
                std::cerr << "Frame prefetch error: " << loader_.LastError() << "\n";
                break;
            }
        }
    }
    if (hasCompareReplay && compareLoader.OrderedSteps().size() <= 256 && initialParticleCapacity <= 12000) {
        compareLoader.SetCacheCapacity(compareLoader.OrderedSteps().size());
        for (std::size_t i = 0; i < compareLoader.OrderedSteps().size(); ++i) {
            if (!compareLoader.PrefetchFrameByOrderedIndex(i)) {
                std::cerr << "Compare frame prefetch error: " << compareLoader.LastError() << "\n";
                break;
            }
        }
    }

    OrbitCamera camera;
    GlRenderer renderer;
    std::string rendererError;
    if (!renderer.Initialize(initialParticleCapacity, &rendererError)) {
        std::cerr << "Failed to initialize renderer: " << rendererError << "\n";
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();
        glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }

    if (runConfig.hasTokamakGeometry) {
        renderer.SetTorusGeometry(runConfig.majorRadius_m, runConfig.minorRadius_m);
    } else {
        renderer.SetTorusGeometry(2.0f, 0.5f);
        std::cerr << "Warning: tokamak geometry not found in run_config_v2.json, using default torus geometry.\n";
    }

    bool paused = false;
    float playbackRate = options_.playbackRate;
    float pointSize = options_.pointSizePixels;
    PhysicsOverlaySettings overlaySettings;
    GraphicsStyleSettings graphicsStyle;
    AnalyticsUiState analyticsUi;
    analyticsUi.compareEnabled = hasCompareReplay;
    analyticsUi.focusEndOrderedIndex = orderedSteps.empty() ? 0 : static_cast<int>(orderedSteps.size() - 1);
    SceneSelectionState selection;
    double displayFrameProgress_s = 0.0;
    ReplayFrame nextFrame;
    bool hasNextFrame = false;
    ReplayFrame compareFrame;
    bool hasCompareFrame = false;
    std::size_t compareFrameIndex = 0;
    int compareMatchedStep = -1;
    double mousePressX = 0.0;
    double mousePressY = 0.0;

    auto loadAdjacentFrame = [&](std::size_t baseIndex) -> bool {
        if (baseIndex + 1 >= orderedSteps.size()) {
            hasNextFrame = false;
            nextFrame = ReplayFrame();
            return true;
        }
        if (!loader_.LoadFrameByOrderedIndex(baseIndex + 1, &nextFrame)) {
            std::cerr << "Frame load error: " << loader_.LastError() << "\n";
            hasNextFrame = false;
            return false;
        }
        hasNextFrame = true;
        return true;
    };

    auto prefetchUpcomingFrames = [&](std::size_t baseIndex) {
        constexpr std::size_t kPrefetchDepth = 12;
        for (std::size_t offset = 2; offset < (2 + kPrefetchDepth); ++offset) {
            const std::size_t prefetchIndex = baseIndex + offset;
            if (prefetchIndex >= orderedSteps.size()) {
                break;
            }
            if (!loader_.PrefetchFrameByOrderedIndex(prefetchIndex)) {
                std::cerr << "Frame prefetch error: " << loader_.LastError() << "\n";
                break;
            }
        }
    };

    auto syncCompareFrame = [&](int step) -> bool {
        if (!hasCompareReplay) {
            hasCompareFrame = false;
            compareFrame = ReplayFrame();
            compareMatchedStep = -1;
            return true;
        }
        const auto& compareSteps = compareLoader.OrderedSteps();
        if (compareSteps.empty()) {
            hasCompareFrame = false;
            compareMatchedStep = -1;
            return true;
        }
        const int matchIndex = FindNearestSeriesIndexForStep(compareSteps, step);
        if (matchIndex < 0) {
            hasCompareFrame = false;
            compareMatchedStep = -1;
            return true;
        }
        const std::size_t orderedIndex = static_cast<std::size_t>(matchIndex);
        compareMatchedStep = compareSteps[orderedIndex];
        if (hasCompareFrame && orderedIndex == compareFrameIndex && compareFrame.step == compareMatchedStep) {
            return true;
        }
        ReplayFrame loadedCompareFrame;
        if (!compareLoader.LoadFrameByOrderedIndex(orderedIndex, &loadedCompareFrame)) {
            std::cerr << "Compare frame load error: " << compareLoader.LastError() << "\n";
            return false;
        }
        compareFrameIndex = orderedIndex;
        compareFrame = std::move(loadedCompareFrame);
        hasCompareFrame = true;
        return true;
    };

    auto jumpToFrame = [&](std::size_t orderedIndex, bool pausePlayback) -> bool {
        ReplayFrame loadedFrame;
        if (!loader_.LoadFrameByOrderedIndex(orderedIndex, &loadedFrame)) {
            std::cerr << "Frame load error: " << loader_.LastError() << "\n";
            return false;
        }
        currentFrameIndex = orderedIndex;
        currentFrame = std::move(loadedFrame);
        displayFrameProgress_s = 0.0;
        if (!loadAdjacentFrame(currentFrameIndex)) {
            return false;
        }
        if (!syncCompareFrame(currentFrame.step)) {
            return false;
        }
        prefetchUpcomingFrames(currentFrameIndex);
        paused = pausePlayback || !hasNextFrame;
        return true;
    };

    if (!jumpToFrame(currentFrameIndex, false)) {
        renderer.Shutdown();
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();
        glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }

    auto lastFrameTime = std::chrono::steady_clock::now();

    bool leftDownPrev = false;

    while (!glfwWindowShouldClose(window)) {
        const auto frameStart = std::chrono::steady_clock::now();
        const double deltaSeconds = std::chrono::duration<double>(frameStart - lastFrameTime).count();
        lastFrameTime = frameStart;

        glfwPollEvents();

        int frameBufferW = 0;
        int frameBufferH = 0;
        glfwGetFramebufferSize(window, &frameBufferW, &frameBufferH);
        glViewport(0, 0, frameBufferW, frameBufferH);
        camera.SetViewport(frameBufferW, frameBufferH);

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        const ImGuiIO& io = ImGui::GetIO();
        double mouseX = 0.0;
        double mouseY = 0.0;
        glfwGetCursorPos(window, &mouseX, &mouseY);
        const bool leftDown = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        bool sceneClickReleased = false;

        if (!io.WantCaptureMouse) {
            if (leftDown && !leftDownPrev) {
                mousePressX = mouseX;
                mousePressY = mouseY;
                camera.BeginRotate(mouseX, mouseY);
            }
            if (!leftDown && leftDownPrev) {
                camera.EndRotate();
                const double dx = mouseX - mousePressX;
                const double dy = mouseY - mousePressY;
                sceneClickReleased = ((dx * dx) + (dy * dy)) <= (7.0 * 7.0);
            }
            camera.OnCursorMove(mouseX, mouseY);
            if (inputState.scrollDeltaY != 0.0) {
                camera.OnScroll(inputState.scrollDeltaY);
                inputState.scrollDeltaY = 0.0;
            }
        } else if (!leftDown && leftDownPrev) {
            camera.EndRotate();
        }
        leftDownPrev = leftDown;

        if (!paused && hasNextFrame) {
            const double physicalFrameSpan_s = std::max(0.0, nextFrame.time_s - currentFrame.time_s);
            const double displayFrameSpan_s = std::max(
                kMinDisplayFrameDuration_s,
                physicalFrameSpan_s) / std::max(0.1, static_cast<double>(playbackRate));
            displayFrameProgress_s = std::min(displayFrameProgress_s + deltaSeconds, displayFrameSpan_s);

            if (displayFrameProgress_s >= displayFrameSpan_s) {
                displayFrameProgress_s = 0.0;
                ++currentFrameIndex;
                currentFrame = std::move(nextFrame);
                if (!loadAdjacentFrame(currentFrameIndex)) {
                    paused = true;
                } else {
                    prefetchUpcomingFrames(currentFrameIndex);
                    if (!syncCompareFrame(currentFrame.step)) {
                        paused = true;
                    }
                    if (!hasNextFrame) {
                        paused = true;
                    }
                }
            }
        }

        float interpolationAlpha = 0.0f;
        double displayedTime_s = currentFrame.time_s;
        double blendFraction = 0.0;
        if (!paused && hasNextFrame) {
            const double physicalFrameSpan_s = std::max(0.0, nextFrame.time_s - currentFrame.time_s);
            const double displayFrameSpan_s = std::max(
                kMinDisplayFrameDuration_s,
                physicalFrameSpan_s) / std::max(0.1, static_cast<double>(playbackRate));
            if (displayFrameSpan_s > 0.0) {
                blendFraction = std::max(0.0, std::min(displayFrameProgress_s / displayFrameSpan_s, 1.0));
            }
            interpolationAlpha = static_cast<float>(blendFraction);
            displayedTime_s = currentFrame.time_s + (blendFraction * physicalFrameSpan_s);
        } else {
            if (!hasNextFrame) {
                paused = true;
            }
            displayFrameProgress_s = 0.0;
        }

        const ReplaySummaryPoint* summary = loader_.SummaryForStep(currentFrame.step);
        const ReplayAnalytics& analytics = loader_.Analytics();
        const std::size_t fusionWindowStartIndex =
            (currentFrameIndex > 4) ? (currentFrameIndex - 4) : 0;
        const ReplaySummaryPoint* fusionWindowStartSummary =
            loader_.SummaryForStep(orderedSteps[fusionWindowStartIndex]);
        const PowerOnVisualState powerOnState = ComputePowerOnVisualState(
            displayedTime_s,
            runConfig,
            summary,
            fusionWindowStartSummary);
        const float ignitionIntensity = powerOnState.ignitionIntensity;
        const std::vector<ReplayFieldProbe>* fieldProbes = loader_.FieldProbesForStep(currentFrame.step);
        const std::vector<ReplayRadialProfileBin>* radialProfile = analytics.RadialProfileForStep(currentFrame.step);
        const ReplayElectrostaticPoint* electrostaticPoint = analytics.ElectrostaticForStep(currentFrame.step);
        const std::vector<ReplaySpeedHistogramBin>* speedHistogram = analytics.SpeedHistogramForStep(currentFrame.step);
        const std::vector<ReplayPitchHistogramBin>* pitchHistogram = analytics.PitchHistogramForStep(currentFrame.step);
        const ReplaySolverResidualPoint* solverResidual = analytics.SolverResidualForStep(currentFrame.step);
        const ReplayFusionReactivityPoint* fusionPoint = analytics.FusionReactivityForStep(currentFrame.step);
        const ReplayWallInteractionPoint* wallPoint = analytics.WallInteractionForStep(currentFrame.step);
        const ReplaySummaryPoint* compareSummary =
            (hasCompareReplay && compareMatchedStep >= 0) ? compareLoader.SummaryForStep(compareMatchedStep) : nullptr;
        const ReplayAnalytics* compareAnalytics =
            (hasCompareReplay && compareMatchedStep >= 0) ? &compareLoader.Analytics() : nullptr;
        const std::vector<ReplayRadialProfileBin>* compareRadialProfile =
            (compareAnalytics != nullptr) ? compareAnalytics->RadialProfileForStep(compareMatchedStep) : nullptr;
        const ReplayWallInteractionPoint* compareWallPoint =
            (compareAnalytics != nullptr) ? compareAnalytics->WallInteractionForStep(compareMatchedStep) : nullptr;
        const ReplayFusionReactivityPoint* compareFusionPoint =
            (compareAnalytics != nullptr) ? compareAnalytics->FusionReactivityForStep(compareMatchedStep) : nullptr;
        const ParticleViewContext particleViewContext = BuildParticleViewContext(
            analyticsUi.particleViewMode,
            currentFrame,
            analyticsUi.slice,
            radialProfile,
            compareRadialProfile,
            wallPoint,
            compareWallPoint,
            majorRadius_m,
            minorRadius_m,
            analyticsUi.compareEnabled && analyticsUi.diffAgainstCompare);
        const Mat4 viewProjection = camera.ViewProjectionMatrix();

        if (selection.pinnedParticle.has_value()) {
            const auto it = std::find_if(
                currentFrame.particles.begin(),
                currentFrame.particles.end(),
                [&](const ReplayParticle& particle) {
                    return particle.particleIndex == selection.pinnedParticle->particleIndex;
                });
            if (it != currentFrame.particles.end()) {
                selection.pinnedParticle = *it;
            } else {
                selection.pinnedParticle.reset();
            }
        }
        if (selection.pinnedProbe.has_value()) {
            bool foundPinnedProbe = false;
            if (fieldProbes != nullptr) {
                for (const ReplayFieldProbe& probe : *fieldProbes) {
                    if (probe.probeIndex == selection.pinnedProbe->probeIndex) {
                        selection.pinnedProbe = probe;
                        foundPinnedProbe = true;
                        break;
                    }
                }
            }
            if (!foundPinnedProbe) {
                selection.pinnedProbe.reset();
            }
        }

        selection.hoveredParticle.reset();
        selection.hoveredProbe.reset();
        if (!io.WantCaptureMouse) {
            const ImVec2 mouse(io.MousePos.x, io.MousePos.y);
            if (const ReplayParticle* hoveredParticle = NearestByScreenDistance(
                    &currentFrame.particles,
                    viewProjection,
                    frameBufferW,
                    frameBufferH,
                    mouse,
                    14.0f,
                    [](const ReplayParticle& particle) { return particle.position_m; })) {
                if (ParticleMatchesSlice(*hoveredParticle, majorRadius_m, analyticsUi.slice) ||
                    (!analyticsUi.slice.enableToroidalSlice && !analyticsUi.slice.enablePoloidalSlice) ||
                    !analyticsUi.slice.dimOutsideSlice) {
                    selection.hoveredParticle = *hoveredParticle;
                }
            }
            if (const ReplayFieldProbe* hoveredProbe = NearestByScreenDistance(
                    fieldProbes,
                    viewProjection,
                    frameBufferW,
                    frameBufferH,
                    mouse,
                    18.0f,
                    [](const ReplayFieldProbe& probe) { return probe.position_m; })) {
                selection.hoveredProbe = *hoveredProbe;
            }
        }

        if (sceneClickReleased && !io.WantCaptureMouse) {
            const ImVec2 mouse(static_cast<float>(mouseX), static_cast<float>(mouseY));
            float particleDistance = std::numeric_limits<float>::infinity();
            float probeDistance = std::numeric_limits<float>::infinity();
            const ReplayParticle* clickedParticle = NearestByScreenDistance(
                &currentFrame.particles,
                viewProjection,
                frameBufferW,
                frameBufferH,
                mouse,
                14.0f,
                [](const ReplayParticle& particle) { return particle.position_m; },
                &particleDistance);
            const ReplayFieldProbe* clickedProbe = NearestByScreenDistance(
                fieldProbes,
                viewProjection,
                frameBufferW,
                frameBufferH,
                mouse,
                18.0f,
                [](const ReplayFieldProbe& probe) { return probe.position_m; },
                &probeDistance);

            if (clickedParticle != nullptr &&
                (clickedProbe == nullptr || particleDistance <= probeDistance)) {
                selection.pinnedParticle = *clickedParticle;
                selection.pinnedProbe.reset();
            } else if (clickedProbe != nullptr) {
                selection.pinnedProbe = *clickedProbe;
                selection.pinnedParticle.reset();
            } else {
                selection.pinnedParticle.reset();
                selection.pinnedProbe.reset();
            }
        }

        std::vector<float> overlayVertices;
        PhysicsOverlayStats overlayStats;
        overlayVertices.reserve(8192);
        AppendFieldProbeArrows(&overlayVertices, fieldProbes, overlaySettings, &overlayStats);
        AppendParticleForceArrows(&overlayVertices, currentFrame, overlaySettings, &overlayStats);
        AppendMagneticGuideLines(&overlayVertices, runConfig, displayedTime_s, overlaySettings, &overlayStats);
        AppendSeededMagneticStreamlines(&overlayVertices, runConfig, displayedTime_s, overlaySettings, selection, &overlayStats);
        AppendSliceGuides(&overlayVertices, runConfig, overlaySettings, analyticsUi.slice);
        renderer.UploadOverlayLines(overlayVertices);
        renderer.UploadFrame(
            currentFrame,
            hasNextFrame ? &nextFrame : nullptr,
            particleViewContext);

        if (ImGui::Begin("Replay Controls")) {
            double sampledWeight = 0.0;
            double sampledWeightedEnergyKeV = 0.0;
            double sampledEnergeticWeight = 0.0;
            double sampledEnergyKeVUnweighted = 0.0;
            double maxSpeed = 0.0;
            std::size_t energeticSamples = 0;
            for (const ReplayParticle& particle : currentFrame.particles) {
                sampledWeight += particle.weight;
                maxSpeed = std::max(maxSpeed, particle.speed_mPerS);
                if (particle.kineticEnergy_keV > 0.0) {
                    sampledEnergyKeVUnweighted += particle.kineticEnergy_keV;
                    ++energeticSamples;
                    if (particle.weight > 0.0) {
                        sampledWeightedEnergyKeV += particle.kineticEnergy_keV * particle.weight;
                        sampledEnergeticWeight += particle.weight;
                    }
                }
            }

            ImGui::TextWrapped(
                "Replaying sampled macro-particles from artifact snapshots. Choose an explicit particle view below "
                "so color always maps to one physical quantity at a time.");
            ImGui::TextWrapped(
                "Drag to orbit, scroll to zoom, click timeline charts to jump in time, and use the focus window "
                "to compare startup, confinement, and burn phases.");
            ImGui::Separator();
            ImGui::Text("Run: %s", loader_.Manifest().runId.c_str());
            if (hasCompareReplay) {
                ImGui::Text("Compare run: %s", compareLoader.Manifest().runId.c_str());
                ImGui::Text("Compare step match: %d", compareMatchedStep);
            }
            ImGui::Text("Snapshots: %zu total, showing %zu / %zu", orderedSteps.size(), currentFrameIndex + 1, orderedSteps.size());
            ImGui::Text("Solver step: %d", currentFrame.step);
            ImGui::Text("Replay time: %.6f s", displayedTime_s);
            ImGui::Text("Snapshot time: %.6f s", currentFrame.time_s);
            ImGui::Text("Visible macro-particles: %zu", currentFrame.particles.size());
            ImGui::Text("Particle slots in snapshot: %llu", static_cast<unsigned long long>(currentFrame.totalParticles));
            ImGui::Text("Represented ions in sample: %.3e", sampledWeight);
            ImGui::Text(
                "Average sampled energy: %.3f keV",
                sampledEnergeticWeight > 0.0
                    ? (sampledWeightedEnergyKeV / sampledEnergeticWeight)
                    : (energeticSamples > 0 ? (sampledEnergyKeVUnweighted / static_cast<double>(energeticSamples)) : 0.0));
            ImGui::Text("Max sampled speed: %.3e m/s", maxSpeed);
            ImGui::Text("Blend to next snapshot: %.0f%%", blendFraction * 100.0);
            ImGui::Text("Core state: %s", powerOnState.label.c_str());
            ImGui::Text("Power ramp complete: %.0f%%", powerOnState.startupProgress * 100.0f);
            ImGui::Text("Fusion onset gate: %.0f%%", powerOnState.fusionGateProgress * 100.0f);
            ImGui::Text("Ignition glow: %.0f%%", ignitionIntensity * 100.0f);

            if (summary != nullptr) {
                ImGui::Separator();
                ImGui::Text("Whole-plasma ions: %llu", static_cast<unsigned long long>(summary->totalIons));
                ImGui::Text("Whole-plasma avg energy: %.6f keV", summary->avgEnergy_keV);
                ImGui::Text("Whole-plasma fusion events: %llu", static_cast<unsigned long long>(summary->fusionEventsTotal));
                ImGui::Text(
                    "Recent fusion events (%zu snapshots): %llu",
                    std::min<std::size_t>(5, currentFrameIndex + 1),
                    static_cast<unsigned long long>(powerOnState.recentFusionEvents));
            } else {
                ImGui::TextUnformatted("Summary row not found for this step.");
            }

            if (ImGui::Button(paused ? "Play" : "Pause")) {
                paused = !paused;
            }
            ImGui::SameLine();
            if (ImGui::Button("Restart")) {
                jumpToFrame(0, false);
            }

            ImGui::SliderFloat("Playback rate", &playbackRate, 0.1f, 8.0f, "%.2fx");
            ImGui::SliderFloat("Point size", &pointSize, 1.0f, 8.0f, "%.1f px");
            const char* particleViewItems[] = {
                "Species",
                "Energy",
                "Pitch Angle",
                "|B|",
                "|E|",
                "Lorentz Accel",
                "Fusion Rate",
                "Wall Loss Risk"};
            int particleViewIndex = static_cast<int>(analyticsUi.particleViewMode);
            if (ImGui::Combo("Particle view", &particleViewIndex, particleViewItems, IM_ARRAYSIZE(particleViewItems))) {
                analyticsUi.particleViewMode = static_cast<ParticleViewMode>(particleViewIndex);
            }
            if (hasCompareReplay) {
                ImGui::Checkbox("Enable compare context", &analyticsUi.compareEnabled);
                ImGui::SameLine();
                ImGui::Checkbox("Diff color mode", &analyticsUi.diffAgainstCompare);
            }
            DrawParticleViewLegend(particleViewContext, currentFrame);
            ImGui::TextWrapped(
                "Fusion rate and wall-loss risk now color the 3D particle cloud directly. When compare mode is on, "
                "diff coloring makes low-vs-high regions immediately visible.");

            ImGui::Separator();
            ImGui::TextWrapped(
                "Click a particle or field probe in the 3D view to pin details and seed local streamlines. "
                "Pinned objects follow the replay when the same id exists in later frames.");
            if (selection.hoveredParticle.has_value()) {
                const ReplayParticle& particle = *selection.hoveredParticle;
                ImGui::Text(
                    "Hover particle %llu: %s | E %.3f keV | pitch %.1f deg | |B| %.3f T | |E| %.3e V/m | |a_L| %.3e m/s^2",
                    static_cast<unsigned long long>(particle.particleIndex),
                    particle.speciesName.c_str(),
                    particle.kineticEnergy_keV,
                    particle.pitchAngle_deg,
                    particle.magneticMagnitude_T,
                    particle.electricMagnitude_VPerM,
                    particle.lorentzAccelerationMagnitude_mPerS2);
            }
            if (selection.hoveredProbe.has_value()) {
                const ReplayFieldProbe& probe = *selection.hoveredProbe;
                ImGui::Text(
                    "Hover probe %llu: phi %.1f deg | theta %.1f deg | rho %.2f | |B| %.3f T | |E| %.3e V/m",
                    static_cast<unsigned long long>(probe.probeIndex),
                    probe.phi_deg,
                    probe.theta_deg,
                    probe.rho,
                    probe.magneticMagnitude_T,
                    probe.electricMagnitude_VPerM);
            }
            if (selection.pinnedParticle.has_value()) {
                const ReplayParticle& particle = *selection.pinnedParticle;
                ImGui::TextColored(
                    ImVec4(1.0f, 0.86f, 0.42f, 1.0f),
                    "Pinned particle %llu: %s | E %.3f keV | pitch %.1f deg | |B| %.3f T | |E| %.3e V/m | |a_L| %.3e m/s^2",
                    static_cast<unsigned long long>(particle.particleIndex),
                    particle.speciesName.c_str(),
                    particle.kineticEnergy_keV,
                    particle.pitchAngle_deg,
                    particle.magneticMagnitude_T,
                    particle.electricMagnitude_VPerM,
                    particle.lorentzAccelerationMagnitude_mPerS2);
            }
            if (selection.pinnedProbe.has_value()) {
                const ReplayFieldProbe& probe = *selection.pinnedProbe;
                ImGui::TextColored(
                    ImVec4(1.0f, 0.86f, 0.42f, 1.0f),
                    "Pinned probe %llu: phi %.1f deg | theta %.1f deg | rho %.2f | |B| %.3f T | |E| %.3e V/m",
                    static_cast<unsigned long long>(probe.probeIndex),
                    probe.phi_deg,
                    probe.theta_deg,
                    probe.rho,
                    probe.magneticMagnitude_T,
                    probe.electricMagnitude_VPerM);
            }
            if ((selection.pinnedParticle.has_value() || selection.pinnedProbe.has_value()) &&
                ImGui::Button("Clear pinned selection")) {
                selection.pinnedParticle.reset();
                selection.pinnedProbe.reset();
            }

            int frameSlider = static_cast<int>(currentFrameIndex);
            if (ImGui::SliderInt("Frame", &frameSlider, 0, static_cast<int>(orderedSteps.size() - 1))) {
                jumpToFrame(static_cast<std::size_t>(frameSlider), true);
            }

            if (!runConfig.hasTokamakGeometry) {
                ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.2f, 1.0f), "Geometry fallback active (default torus).");
            }
        }
        ImGui::End();

        if (ImGui::Begin("Physics Overlays")) {
            const TokamakConfig effectiveTokamak = EffectiveTokamakConfigForViewer(runConfig, displayedTime_s);

            ImGui::TextWrapped(
                "These overlays expose the fields and Lorentz forces that drive confinement and fusion. "
                "Magnetic guide lines show the helical field geometry, field arrows show local strength/direction, "
                "and force arrows show how visible macro-particles are being accelerated.");
            ImGui::Separator();
            ImGui::Checkbox("Show magnetic arrows", &overlaySettings.showMagneticArrows);
            ImGui::Checkbox("Show electric arrows", &overlaySettings.showElectricArrows);
            ImGui::Checkbox("Show magnetic guide lines", &overlaySettings.showMagneticGuideLines);
            ImGui::Checkbox("Show particle force arrows", &overlaySettings.showParticleForceArrows);
            ImGui::Checkbox("Show seeded streamlines", &overlaySettings.showSeededStreamlines);
            ImGui::Checkbox("Show slice guides", &overlaySettings.showSliceGuides);
            ImGui::SliderFloat("Magnetic arrow scale", &overlaySettings.magneticArrowScale, 0.05f, 0.60f, "%.2f");
            ImGui::SliderFloat("Electric arrow scale", &overlaySettings.electricArrowScale, 0.05f, 0.60f, "%.2f");
            ImGui::SliderFloat("Force arrow scale", &overlaySettings.forceArrowScale, 0.05f, 0.80f, "%.2f");
            ImGui::SliderInt("Force arrow budget", &overlaySettings.maxParticleForceArrows, 8, 160);

            ImGui::Separator();
            ImGui::Text("Field probes in frame: %zu", overlayStats.fieldProbeCount);
            ImGui::Text("Magnetic guide segments: %zu", overlayStats.guideLineCount);
            ImGui::Text("Seeded streamline segments: %zu", overlayStats.seededGuideLineCount);
            ImGui::Text("Particle force arrows: %zu", overlayStats.forceArrowCount);
            ImGui::Text("Peak |B| in probes: %.3f T", overlayStats.maxMagneticMagnitude_T);
            ImGui::Text("Peak |E| in probes: %.3e V/m", overlayStats.maxElectricMagnitude_VPerM);
            ImGui::Text("Peak |a_L| on visible particles: %.3e m/s^2", overlayStats.maxLorentzAcceleration_mPerS2);

            ImGui::Separator();
            ImGui::Text("Toroidal current: %.2f MA", effectiveTokamak.toroidalCurrent_A / 1.0e6f);
            ImGui::Text("Plasma current: %.2f MA", effectiveTokamak.plasmaCurrent_A / 1.0e6f);
            if (runConfig.hasNbiConfig) {
                ImGui::Text("Beam energy: %.1f keV", runConfig.nbiConfig.beamEnergy_keV);
                ImGui::Text(
                    "Injector throughput: %.0f%% of %d pairs/step",
                    powerOnState.startupProgress * 100.0f,
                    runConfig.nbiConfig.particlesPerStep);
            }
            if (runConfig.hasPlasmaCurrentProfile) {
                ImGui::Text(
                    "Current profile: %s",
                    PlasmaCurrentProfileKindName(runConfig.plasmaCurrentProfile.kind));
            }
            if (runConfig.hasElectricFieldMode) {
                ImGui::Text(
                    "Electric field mode: %s",
                    ElectricFieldModeName(runConfig.electricFieldMode));
            }

            if (fieldProbes == nullptr) {
                ImGui::TextColored(
                    ImVec4(1.0f, 0.75f, 0.25f, 1.0f),
                    "This replay does not include exported field probes for the current frame.");
            }
            if (!CanDrawAnalyticMagneticGuideLines(runConfig)) {
                ImGui::TextColored(
                    ImVec4(1.0f, 0.75f, 0.25f, 1.0f),
                    "Analytic magnetic guide lines are unavailable for this replay configuration.");
            }
        }
        ImGui::End();

        if (ImGui::Begin("Graphics Tuning")) {
            ImGui::TextWrapped(
                "These controls keep the upgraded look without washing out the reactor physics. "
                "If color encodings or electric-field overlays feel muted, lower the shell/fog and keep overlays in x-ray mode.");
            if (ImGui::Button("Restore physics-first defaults")) {
                graphicsStyle = GraphicsStyleSettings{};
            }
            ImGui::Separator();
            ImGui::Checkbox("Show vessel shell", &graphicsStyle.showVesselShell);
            ImGui::Checkbox("Show particle trails", &graphicsStyle.showParticleTrails);
            ImGui::Checkbox("Show density splat", &graphicsStyle.showDensitySplat);
            ImGui::Checkbox("X-ray overlays", &graphicsStyle.xrayOverlays);
            ImGui::SliderFloat("Shell opacity", &graphicsStyle.shellOpacity, 0.0f, 0.16f, "%.2f");
            ImGui::SliderFloat("Scene fog", &graphicsStyle.fogStrength, 0.0f, 0.50f, "%.2f");
            ImGui::SliderFloat("Particle fog", &graphicsStyle.pointFogStrength, 0.0f, 0.35f, "%.2f");
            ImGui::SliderFloat("Scene line opacity", &graphicsStyle.sceneLineOpacity, 0.10f, 1.00f, "%.2f");
            ImGui::SliderFloat("Overlay opacity", &graphicsStyle.overlayLineOpacity, 0.25f, 1.00f, "%.2f");
            ImGui::SliderFloat("Overlay fog", &graphicsStyle.overlayFogStrength, 0.0f, 0.20f, "%.2f");
            ImGui::SliderFloat("Trail opacity", &graphicsStyle.trailOpacity, 0.0f, 1.0f, "%.2f");
            ImGui::SliderFloat("Density splat opacity", &graphicsStyle.densitySplatOpacity, 0.0f, 0.20f, "%.2f");
            ImGui::SliderFloat("Density splat scale", &graphicsStyle.densitySplatScale, 1.0f, 4.0f, "%.1f");
            ImGui::SliderFloat("Particle color boost", &graphicsStyle.particleColorBoost, 0.8f, 2.4f, "%.2f");
            ImGui::SliderFloat("Scene line width", &graphicsStyle.sceneLineWidth, 1.0f, 3.0f, "%.1f");
            ImGui::SliderFloat("Overlay line width", &graphicsStyle.overlayLineWidth, 1.0f, 4.0f, "%.1f");
            ImGui::SliderFloat("Trail line width", &graphicsStyle.trailLineWidth, 1.0f, 4.0f, "%.1f");
        }
        ImGui::End();

        if (ImGui::Begin("Slices + Diff")) {
            ImGui::TextWrapped(
                "Slice controls turn the torus into a teachable cross-section. The dials are draggable: move the "
                "center angle directly, or hold Shift while dragging to widen the slice window.");
            ImGui::Separator();
            ImGui::Checkbox("Enable toroidal slice", &analyticsUi.slice.enableToroidalSlice);
            ImGui::Checkbox("Enable poloidal slice", &analyticsUi.slice.enablePoloidalSlice);
            ImGui::Checkbox("Dim outside slice", &analyticsUi.slice.dimOutsideSlice);
            if (analyticsUi.slice.enableToroidalSlice) {
                DrawAngleDial("Toroidal plane", &analyticsUi.slice.toroidalCenter_deg, &analyticsUi.slice.toroidalHalfWidth_deg, 90.0f);
            }
            if (analyticsUi.slice.enablePoloidalSlice) {
                DrawAngleDial("Poloidal plane", &analyticsUi.slice.poloidalCenter_deg, &analyticsUi.slice.poloidalHalfWidth_deg, 90.0f);
            }
            if (hasCompareReplay) {
                ImGui::Separator();
                ImGui::TextWrapped(
                    "Run-vs-run diff mode aligns the compare replay by nearest solver step and feeds the same data into "
                    "the 3D coloring, radial heatmap, and delta readouts.");
                ImGui::Checkbox("Use compare replay", &analyticsUi.compareEnabled);
                ImGui::Checkbox("Show heatmap delta", &analyticsUi.heatmapShowsDelta);
                ImGui::Text("Compare matched step: %d", compareMatchedStep);
                if (summary != nullptr && compareSummary != nullptr) {
                    ImGui::Text(
                        "Avg energy delta: %.3f keV | fusion total delta: %+lld",
                        summary->avgEnergy_keV - compareSummary->avgEnergy_keV,
                        static_cast<long long>(summary->fusionEventsTotal) -
                            static_cast<long long>(compareSummary->fusionEventsTotal));
                }
                if (fusionPoint != nullptr && compareFusionPoint != nullptr) {
                    ImGui::Text(
                        "Fusion accepted / step delta: %+lld | probability delta %.3e",
                        static_cast<long long>(fusionPoint->fusionAcceptedStep) -
                            static_cast<long long>(compareFusionPoint->fusionAcceptedStep),
                        fusionPoint->avgProbabilityStep - compareFusionPoint->avgProbabilityStep);
                }
                if (wallPoint != nullptr && compareWallPoint != nullptr) {
                    ImGui::Text(
                        "Wall loss delta: %.3e | wall hits delta: %+lld",
                        wallPoint->wallLossWeightStep - compareWallPoint->wallLossWeightStep,
                        static_cast<long long>(wallPoint->wallHitCountStep) -
                            static_cast<long long>(compareWallPoint->wallHitCountStep));
                }
            }
        }
        ImGui::End();

        if (analyticsUi.focusEndOrderedIndex < 0) {
            analyticsUi.focusEndOrderedIndex = orderedSteps.empty() ? 0 : static_cast<int>(orderedSteps.size() - 1);
        }
        analyticsUi.focusStartOrderedIndex = std::max(0, std::min(analyticsUi.focusStartOrderedIndex, static_cast<int>(orderedSteps.size() - 1)));
        analyticsUi.focusEndOrderedIndex = std::max(analyticsUi.focusStartOrderedIndex, std::min(
            analyticsUi.focusEndOrderedIndex,
            static_cast<int>(orderedSteps.size() - 1)));

        if (ImGui::Begin("Physics Analytics")) {
            ImGui::TextWrapped(
                "Linked charts turn the replay into a reactor-learning surface. Click a timeline to jump to that "
                "phase of the run. The highlighted window lets you compare startup, confinement, and burn periods.");
            ImGui::Separator();

            ImGui::SliderInt("Focus start", &analyticsUi.focusStartOrderedIndex, 0, static_cast<int>(orderedSteps.size() - 1));
            ImGui::SliderInt(
                "Focus end",
                &analyticsUi.focusEndOrderedIndex,
                analyticsUi.focusStartOrderedIndex,
                static_cast<int>(orderedSteps.size() - 1));
            const int focusStartStep = orderedSteps[analyticsUi.focusStartOrderedIndex];
            const int focusEndStep = orderedSteps[analyticsUi.focusEndOrderedIndex];
            ImGui::Text("Focus window: steps %d -> %d", focusStartStep, focusEndStep);

            auto jumpToStep = [&](int step) {
                const auto it = std::lower_bound(orderedSteps.begin(), orderedSteps.end(), step);
                if (it == orderedSteps.end()) {
                    return;
                }
                jumpToFrame(static_cast<std::size_t>(std::distance(orderedSteps.begin(), it)), true);
            };

            auto drawTimelineTooltip =
                [&](const char* label, const std::vector<int>& steps, const std::vector<double>& values, const PlotInteraction& interaction, const char* units) {
                    if (interaction.hoveredIndex < 0 ||
                        interaction.hoveredIndex >= static_cast<int>(steps.size()) ||
                        interaction.hoveredIndex >= static_cast<int>(values.size())) {
                        return;
                    }
                    ImGui::BeginTooltip();
                    ImGui::Text("%s", label);
                    ImGui::Text("Step: %d", steps[interaction.hoveredIndex]);
                    ImGui::Text("Value: %.6e %s", values[interaction.hoveredIndex], units);
                    ImGui::EndTooltip();
                };

            std::vector<int> magneticSteps;
            std::vector<double> magneticValues;
            BuildTimelineSeries(
                analytics.MagneticFieldSeries(),
                [](const ReplayMagneticFieldStepSummary& point) { return point.stepMaxField_T; },
                &magneticSteps,
                &magneticValues);
            std::vector<int> electricSteps;
            std::vector<double> electricValues;
            BuildTimelineSeries(
                analytics.ElectrostaticSeries(),
                [](const ReplayElectrostaticPoint& point) { return point.maxElectricField_VPerM; },
                &electricSteps,
                &electricValues);
            std::vector<int> residualSteps;
            std::vector<double> residualValues;
            BuildTimelineSeries(
                analytics.SolverResidualSeries(),
                [](const ReplaySolverResidualPoint& point) { return point.residualL2; },
                &residualSteps,
                &residualValues);
            std::vector<int> wallSteps;
            std::vector<double> wallValues;
            BuildTimelineSeries(
                analytics.WallInteractionSeries(),
                [](const ReplayWallInteractionPoint& point) { return point.wallLossWeightStep; },
                &wallSteps,
                &wallValues);
            std::vector<int> fusionSteps;
            std::vector<double> fusionValues;
            BuildTimelineSeries(
                analytics.FusionReactivitySeries(),
                [](const ReplayFusionReactivityPoint& point) { return static_cast<double>(point.fusionAcceptedStep); },
                &fusionSteps,
                &fusionValues);

            const PlotInteraction magneticPlot = PlotSeriesWithOverlay(
                "Max |B| timeline##magnetic_timeline",
                magneticValues,
                FindNearestSeriesIndexForStep(magneticSteps, currentFrame.step),
                FindNearestSeriesIndexForStep(magneticSteps, focusStartStep),
                FindNearestSeriesIndexForStep(magneticSteps, focusEndStep),
                "Tesla",
                ImVec2(-1.0f, 76.0f));
            if (magneticPlot.clickedIndex >= 0 && magneticPlot.clickedIndex < static_cast<int>(magneticSteps.size())) {
                jumpToStep(magneticSteps[magneticPlot.clickedIndex]);
            }
            drawTimelineTooltip("Max |B|", magneticSteps, magneticValues, magneticPlot, "T");

            const PlotInteraction electricPlot = PlotSeriesWithOverlay(
                "Max |E| timeline##electric_timeline",
                electricValues,
                FindNearestSeriesIndexForStep(electricSteps, currentFrame.step),
                FindNearestSeriesIndexForStep(electricSteps, focusStartStep),
                FindNearestSeriesIndexForStep(electricSteps, focusEndStep),
                "V/m",
                ImVec2(-1.0f, 76.0f));
            if (electricPlot.clickedIndex >= 0 && electricPlot.clickedIndex < static_cast<int>(electricSteps.size())) {
                jumpToStep(electricSteps[electricPlot.clickedIndex]);
            }
            drawTimelineTooltip("Max |E|", electricSteps, electricValues, electricPlot, "V/m");

            const PlotInteraction residualPlot = PlotSeriesWithOverlay(
                "Solver residual timeline##residual_timeline",
                residualValues,
                FindNearestSeriesIndexForStep(residualSteps, currentFrame.step),
                FindNearestSeriesIndexForStep(residualSteps, focusStartStep),
                FindNearestSeriesIndexForStep(residualSteps, focusEndStep),
                "L2",
                ImVec2(-1.0f, 76.0f));
            if (residualPlot.clickedIndex >= 0 && residualPlot.clickedIndex < static_cast<int>(residualSteps.size())) {
                jumpToStep(residualSteps[residualPlot.clickedIndex]);
            }
            drawTimelineTooltip("Solver residual", residualSteps, residualValues, residualPlot, "L2");

            const PlotInteraction wallPlot = PlotSeriesWithOverlay(
                "Wall-loss trend##wall_timeline",
                wallValues,
                FindNearestSeriesIndexForStep(wallSteps, currentFrame.step),
                FindNearestSeriesIndexForStep(wallSteps, focusStartStep),
                FindNearestSeriesIndexForStep(wallSteps, focusEndStep),
                "weight",
                ImVec2(-1.0f, 76.0f));
            if (wallPlot.clickedIndex >= 0 && wallPlot.clickedIndex < static_cast<int>(wallSteps.size())) {
                jumpToStep(wallSteps[wallPlot.clickedIndex]);
            }
            drawTimelineTooltip("Wall loss", wallSteps, wallValues, wallPlot, "weight");

            const PlotInteraction fusionPlot = PlotSeriesWithOverlay(
                "Fusion accepted / step##fusion_timeline",
                fusionValues,
                FindNearestSeriesIndexForStep(fusionSteps, currentFrame.step),
                FindNearestSeriesIndexForStep(fusionSteps, focusStartStep),
                FindNearestSeriesIndexForStep(fusionSteps, focusEndStep),
                "events",
                ImVec2(-1.0f, 76.0f));
            if (fusionPlot.clickedIndex >= 0 && fusionPlot.clickedIndex < static_cast<int>(fusionSteps.size())) {
                jumpToStep(fusionSteps[fusionPlot.clickedIndex]);
            }
            drawTimelineTooltip("Fusion accepted / step", fusionSteps, fusionValues, fusionPlot, "events");

            ImGui::Separator();
            const char* radialMetricItems[] = {"Density", "Avg Ion Energy", "Fusion Rate"};
            int radialMetricIndex = static_cast<int>(analyticsUi.radialMetric);
            if (ImGui::Combo("Radial metric", &radialMetricIndex, radialMetricItems, IM_ARRAYSIZE(radialMetricItems))) {
                analyticsUi.radialMetric = static_cast<RadialChartMetric>(radialMetricIndex);
            }

            if (radialProfile != nullptr && !radialProfile->empty()) {
                std::vector<double> radialValues;
                radialValues.reserve(radialProfile->size());
                for (const ReplayRadialProfileBin& bin : *radialProfile) {
                    radialValues.push_back(RadialMetricValue(bin, analyticsUi.radialMetric));
                }

                if (analyticsUi.selectedRadialBin >= static_cast<int>(radialProfile->size())) {
                    analyticsUi.selectedRadialBin = static_cast<int>(radialProfile->size() - 1);
                }
                if (analyticsUi.selectedRadialBin < 0) {
                    analyticsUi.selectedRadialBin = 0;
                }
                const PlotInteraction radialPlot = PlotHistogramWithOverlay(
                    "Current-step radial profile##radial_profile",
                    radialValues,
                    analyticsUi.selectedRadialBin,
                    RadialMetricUnits(analyticsUi.radialMetric),
                    ImVec2(-1.0f, 88.0f));
                if (radialPlot.clickedIndex >= 0) {
                    analyticsUi.selectedRadialBin = radialPlot.clickedIndex;
                }
                if (radialPlot.hoveredIndex >= 0 && radialPlot.hoveredIndex < static_cast<int>(radialProfile->size())) {
                    const ReplayRadialProfileBin& hoveredBin = (*radialProfile)[radialPlot.hoveredIndex];
                    ImGui::BeginTooltip();
                    ImGui::Text("Bin %d", hoveredBin.binIndex);
                    ImGui::Text("r = %.3f -> %.3f m", hoveredBin.rInner_m, hoveredBin.rOuter_m);
                    ImGui::Text(
                        "%s: %.6e %s",
                        RadialChartMetricName(analyticsUi.radialMetric),
                        RadialMetricValue(hoveredBin, analyticsUi.radialMetric),
                        RadialMetricUnits(analyticsUi.radialMetric));
                    ImGui::EndTooltip();
                }

                if (analyticsUi.selectedRadialBin >= 0 && analyticsUi.selectedRadialBin < static_cast<int>(radialProfile->size())) {
                    const ReplayRadialProfileBin& pinnedBin = (*radialProfile)[analyticsUi.selectedRadialBin];
                    const double compareDensity = CompareRadialMetricValue(compareRadialProfile, pinnedBin.binIndex, RadialChartMetric::Density);
                    const double compareTemperature = CompareRadialMetricValue(compareRadialProfile, pinnedBin.binIndex, RadialChartMetric::Temperature);
                    const double compareFusion = CompareRadialMetricValue(compareRadialProfile, pinnedBin.binIndex, RadialChartMetric::FusionRate);
                    ImGui::Text(
                        "Pinned radial bin %d: r = %.3f -> %.3f m | density %.3e m^-3 | avg energy %.3f keV | fusion rate %.3e m^-3 s^-1",
                        pinnedBin.binIndex,
                        pinnedBin.rInner_m,
                        pinnedBin.rOuter_m,
                        pinnedBin.density_m3,
                        pinnedBin.avgIonEnergy_keV,
                        pinnedBin.fusionRatePlaceholder ? 0.0 : pinnedBin.fusionRate_m3_s);
                    if (analyticsUi.compareEnabled && compareRadialProfile != nullptr) {
                        ImGui::Text(
                            "Compare deltas: density %+.3e | avg energy %+.3f keV | fusion rate %+.3e",
                            pinnedBin.density_m3 - compareDensity,
                            pinnedBin.avgIonEnergy_keV - compareTemperature,
                            (pinnedBin.fusionRatePlaceholder ? 0.0 : pinnedBin.fusionRate_m3_s) - compareFusion);
                    }
                }

                const HeatmapInteraction heatmapInteraction = DrawRadialCrossSectionHeatmap(
                    *radialProfile,
                    compareRadialProfile,
                    analyticsUi.radialMetric,
                    analyticsUi.compareEnabled && analyticsUi.heatmapShowsDelta,
                    analyticsUi.slice,
                    analyticsUi.selectedRadialBin);
                if (heatmapInteraction.clickedBin >= 0) {
                    analyticsUi.selectedRadialBin = heatmapInteraction.clickedBin;
                }
                if (heatmapInteraction.hoveredBin >= 0 &&
                    heatmapInteraction.hoveredBin < static_cast<int>(radialProfile->size())) {
                    const ReplayRadialProfileBin& hoveredBin = (*radialProfile)[heatmapInteraction.hoveredBin];
                    ImGui::BeginTooltip();
                    ImGui::Text("Cross-section bin %d", hoveredBin.binIndex);
                    ImGui::Text("r = %.3f -> %.3f m", hoveredBin.rInner_m, hoveredBin.rOuter_m);
                    ImGui::Text(
                        "%s: %.6e %s",
                        RadialChartMetricName(analyticsUi.radialMetric),
                        RadialMetricValue(hoveredBin, analyticsUi.radialMetric),
                        RadialMetricUnits(analyticsUi.radialMetric));
                    ImGui::EndTooltip();
                }
            } else {
                ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.25f, 1.0f), "No radial profile data is available for this step.");
            }

            if (speedHistogram != nullptr && !speedHistogram->empty()) {
                std::vector<double> speedValues;
                speedValues.reserve(speedHistogram->size());
                for (const ReplaySpeedHistogramBin& bin : *speedHistogram) {
                    speedValues.push_back(static_cast<double>(bin.count));
                }
                const PlotInteraction speedPlot = PlotHistogramWithOverlay(
                    "Speed histogram##speed_histogram",
                    speedValues,
                    -1,
                    "counts",
                    ImVec2(-1.0f, 70.0f));
                if (speedPlot.hoveredIndex >= 0 && speedPlot.hoveredIndex < static_cast<int>(speedHistogram->size())) {
                    const ReplaySpeedHistogramBin& bin = (*speedHistogram)[speedPlot.hoveredIndex];
                    ImGui::BeginTooltip();
                    ImGui::Text("Speed bin %d", bin.binIndex);
                    ImGui::Text("%.3e -> %.3e m/s", bin.speedMin_mPerS, bin.speedMax_mPerS);
                    ImGui::Text("Count: %llu / %llu", static_cast<unsigned long long>(bin.count), static_cast<unsigned long long>(bin.totalSamples));
                    ImGui::EndTooltip();
                }
            }

            if (pitchHistogram != nullptr && !pitchHistogram->empty()) {
                std::vector<double> pitchValues;
                pitchValues.reserve(pitchHistogram->size());
                for (const ReplayPitchHistogramBin& bin : *pitchHistogram) {
                    pitchValues.push_back(static_cast<double>(bin.count));
                }
                const PlotInteraction pitchPlot = PlotHistogramWithOverlay(
                    "Pitch histogram##pitch_histogram",
                    pitchValues,
                    -1,
                    "counts",
                    ImVec2(-1.0f, 70.0f));
                if (pitchPlot.hoveredIndex >= 0 && pitchPlot.hoveredIndex < static_cast<int>(pitchHistogram->size())) {
                    const ReplayPitchHistogramBin& bin = (*pitchHistogram)[pitchPlot.hoveredIndex];
                    ImGui::BeginTooltip();
                    ImGui::Text("Pitch bin %d", bin.binIndex);
                    ImGui::Text("%.1f -> %.1f deg", bin.pitchMin_deg, bin.pitchMax_deg);
                    ImGui::Text("Count: %llu / %llu", static_cast<unsigned long long>(bin.count), static_cast<unsigned long long>(bin.totalSamples));
                    ImGui::EndTooltip();
                }
            }

            if (fieldProbes != nullptr && !fieldProbes->empty()) {
                std::vector<double> rhoBuckets;
                rhoBuckets.reserve(fieldProbes->size());
                for (const ReplayFieldProbe& probe : *fieldProbes) {
                    bool seen = false;
                    for (double existing : rhoBuckets) {
                        if (std::fabs(existing - probe.rho) <= 1.0e-4) {
                            seen = true;
                            break;
                        }
                    }
                    if (!seen) {
                        rhoBuckets.push_back(probe.rho);
                    }
                }
                std::sort(rhoBuckets.begin(), rhoBuckets.end());
                if (!rhoBuckets.empty()) {
                    analyticsUi.unwrapRhoBucket = std::max(0, std::min(analyticsUi.unwrapRhoBucket, static_cast<int>(rhoBuckets.size() - 1)));
                    const char* unwrapItems[] = {"|B|", "|E|"};
                    int unwrapMetricIndex = static_cast<int>(analyticsUi.unwrapMetric);
                    if (ImGui::Combo("Toroidal unwrap metric", &unwrapMetricIndex, unwrapItems, IM_ARRAYSIZE(unwrapItems))) {
                        analyticsUi.unwrapMetric = static_cast<ToroidalUnwrapMetric>(unwrapMetricIndex);
                    }
                    ImGui::SliderInt("Unwrap rho shell", &analyticsUi.unwrapRhoBucket, 0, static_cast<int>(rhoBuckets.size() - 1));
                    const ToroidalUnwrapInteraction unwrapInteraction = DrawToroidalUnwrap(
                        *fieldProbes,
                        rhoBuckets,
                        analyticsUi.unwrapRhoBucket,
                        analyticsUi.unwrapMetric,
                        selection.pinnedProbe);
                    if (unwrapInteraction.clickedProbeIndex >= 0 &&
                        unwrapInteraction.clickedProbeIndex < static_cast<int>(fieldProbes->size())) {
                        selection.pinnedProbe = (*fieldProbes)[static_cast<std::size_t>(unwrapInteraction.clickedProbeIndex)];
                        selection.pinnedParticle.reset();
                    }
                    if (unwrapInteraction.hoveredProbeIndex >= 0 &&
                        unwrapInteraction.hoveredProbeIndex < static_cast<int>(fieldProbes->size())) {
                        const ReplayFieldProbe& probe = (*fieldProbes)[static_cast<std::size_t>(unwrapInteraction.hoveredProbeIndex)];
                        ImGui::BeginTooltip();
                        ImGui::Text("Probe %llu", static_cast<unsigned long long>(probe.probeIndex));
                        ImGui::Text("phi %.1f deg | theta %.1f deg | rho %.2f", probe.phi_deg, probe.theta_deg, probe.rho);
                        ImGui::Text("|B| %.3f T | |E| %.3e V/m", probe.magneticMagnitude_T, probe.electricMagnitude_VPerM);
                        ImGui::EndTooltip();
                    }
                }
            }

            ImGui::Separator();
            if (electrostaticPoint != nullptr) {
                ImGui::Text(
                    "Current electrostatics: mode %s | max |E| %.3e V/m | mean |E| %.3e V/m | iterations %u",
                    electrostaticPoint->electricFieldMode.c_str(),
                    electrostaticPoint->maxElectricField_VPerM,
                    electrostaticPoint->meanElectricField_VPerM,
                    electrostaticPoint->solverIterations);
            }
            if (solverResidual != nullptr) {
                ImGui::Text(
                    "Current solver residual: %.3e | status %s | tolerance %.3e",
                    solverResidual->residualL2,
                    solverResidual->status.c_str(),
                    solverResidual->tolerance);
            }
            if (fusionPoint != nullptr) {
                ImGui::Text(
                    "Current fusion kinetics: accepted %llu / %llu attempts | avg probability %.3e | avg sigma %.3e m^2",
                    static_cast<unsigned long long>(fusionPoint->fusionAcceptedStep),
                    static_cast<unsigned long long>(fusionPoint->fusionAttemptsStep),
                    fusionPoint->avgProbabilityStep,
                    fusionPoint->avgSigma_m2Step);
            }
            if (wallPoint != nullptr) {
                ImGui::Text(
                    "Current wall trend: step hits %llu | step loss %.3e | total loss %.3e",
                    static_cast<unsigned long long>(wallPoint->wallHitCountStep),
                    wallPoint->wallLossWeightStep,
                    wallPoint->wallLossWeightTotal);
            }
        }
        ImGui::End();

        glClearColor(
            0.03f + (0.06f * ignitionIntensity),
            0.05f + (0.01f * ignitionIntensity),
            0.08f + (0.10f * ignitionIntensity),
            1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        renderer.Draw(viewProjection, pointSize, interpolationAlpha, ignitionIntensity, graphicsStyle);

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    renderer.Shutdown();
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}

}  // namespace tokamak::viewer
