#include "tokamak/viewer/gl_renderer.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#include <glad/glad.h>

namespace tokamak::viewer {
namespace {

unsigned int CompileShader(unsigned int shaderType, const char* source, std::string* errorOut) {
    const unsigned int shader = glCreateShader(shaderType);
    glShaderSource(shader, 1, &source, nullptr);
    glCompileShader(shader);

    int compiled = GL_FALSE;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
    if (compiled == GL_TRUE) {
        return shader;
    }

    int logLength = 0;
    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &logLength);
    std::string log(static_cast<std::size_t>(std::max(logLength, 1)), '\0');
    glGetShaderInfoLog(shader, logLength, nullptr, log.data());

    glDeleteShader(shader);
    if (errorOut != nullptr) {
        *errorOut = "Shader compile failed: " + log;
    }
    return 0;
}

unsigned int CreateProgram(const char* vertexSource, const char* fragmentSource, std::string* errorOut) {
    const unsigned int vs = CompileShader(GL_VERTEX_SHADER, vertexSource, errorOut);
    if (vs == 0) {
        return 0;
    }

    const unsigned int fs = CompileShader(GL_FRAGMENT_SHADER, fragmentSource, errorOut);
    if (fs == 0) {
        glDeleteShader(vs);
        return 0;
    }

    const unsigned int program = glCreateProgram();
    glAttachShader(program, vs);
    glAttachShader(program, fs);
    glLinkProgram(program);

    glDeleteShader(vs);
    glDeleteShader(fs);

    int linked = GL_FALSE;
    glGetProgramiv(program, GL_LINK_STATUS, &linked);
    if (linked == GL_TRUE) {
        return program;
    }

    int logLength = 0;
    glGetProgramiv(program, GL_INFO_LOG_LENGTH, &logLength);
    std::string log(static_cast<std::size_t>(std::max(logLength, 1)), '\0');
    glGetProgramInfoLog(program, logLength, nullptr, log.data());
    glDeleteProgram(program);

    if (errorOut != nullptr) {
        *errorOut = "Program link failed: " + log;
    }
    return 0;
}

const char* kPointVertexShader = R"GLSL(
#version 330 core
layout(location = 0) in vec3 aPosition;
layout(location = 1) in vec3 aPositionDelta;
layout(location = 2) in vec3 aColor;

uniform mat4 uViewProjection;
uniform float uPointSize;
uniform float uInterpolationAlpha;

out vec3 vColor;
out float vDepth01;

void main() {
    vec3 simulatedPosition = aPosition + (aPositionDelta * uInterpolationAlpha);
    gl_Position = uViewProjection * vec4(simulatedPosition, 1.0);
    gl_PointSize = uPointSize;
    vColor = aColor;
    vDepth01 = clamp((gl_Position.z / max(gl_Position.w, 1.0e-5)) * 0.5 + 0.5, 0.0, 1.0);
}
)GLSL";

const char* kPointFragmentShader = R"GLSL(
#version 330 core
in vec3 vColor;
in float vDepth01;
out vec4 fragColor;

uniform vec3 uFogColor;
uniform float uAlphaScale;
uniform float uColorBoost;
uniform float uFogMix;

void main() {
    vec2 centered = gl_PointCoord * 2.0 - 1.0;
    float radiusSq = dot(centered, centered);
    if (radiusSq > 1.0) {
        discard;
    }
    float glow = exp(-3.4 * radiusSq);
    float alpha = uAlphaScale * mix(0.08, 0.96, glow);
    vec3 color = vColor * (0.42 + (uColorBoost * glow));
    float fog = smoothstep(0.32, 0.96, vDepth01) * uFogMix;
    color = mix(color, uFogColor, fog);
    fragColor = vec4(color, alpha);
}
)GLSL";

const char* kLineVertexShader = R"GLSL(
#version 330 core
layout(location = 0) in vec3 aPosition;
layout(location = 1) in vec3 aColor;

uniform mat4 uViewProjection;
out vec3 vColor;
out float vDepth01;

void main() {
    gl_Position = uViewProjection * vec4(aPosition, 1.0);
    vColor = aColor;
    vDepth01 = clamp((gl_Position.z / max(gl_Position.w, 1.0e-5)) * 0.5 + 0.5, 0.0, 1.0);
}
)GLSL";

const char* kLineFragmentShader = R"GLSL(
#version 330 core
in vec3 vColor;
in float vDepth01;
out vec4 fragColor;

uniform vec3 uFogColor;
uniform float uAlphaScale;
uniform float uFogMix;

void main() {
    vec3 color = mix(vColor, uFogColor, smoothstep(0.40, 0.98, vDepth01) * uFogMix);
    fragColor = vec4(color, uAlphaScale);
}
)GLSL";

const char* kShellVertexShader = R"GLSL(
#version 330 core
layout(location = 0) in vec3 aPosition;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in vec3 aColor;

uniform mat4 uViewProjection;

out vec3 vColor;
out vec3 vNormal;
out float vDepth01;

void main() {
    gl_Position = uViewProjection * vec4(aPosition, 1.0);
    vColor = aColor;
    vNormal = normalize(aNormal);
    vDepth01 = clamp((gl_Position.z / max(gl_Position.w, 1.0e-5)) * 0.5 + 0.5, 0.0, 1.0);
}
)GLSL";

const char* kShellFragmentShader = R"GLSL(
#version 330 core
in vec3 vColor;
in vec3 vNormal;
in float vDepth01;
out vec4 fragColor;

uniform vec3 uFogColor;
uniform float uAlphaScale;
uniform float uFogMix;

void main() {
    vec3 lightDir = normalize(vec3(0.35, 0.25, 1.0));
    float diffuse = 0.35 + 0.65 * max(dot(normalize(vNormal), lightDir), 0.0);
    float rim = pow(1.0 - abs(vNormal.z), 1.8);
    vec3 color = vColor * diffuse + vec3(0.14, 0.18, 0.22) * rim;
    color = mix(color, uFogColor, smoothstep(0.42, 0.98, vDepth01) * uFogMix);
    fragColor = vec4(color, uAlphaScale);
}
)GLSL";

void PushVertex(std::vector<float>* out, const Vec3& position, float r, float g, float b) {
    out->push_back(position.x);
    out->push_back(position.y);
    out->push_back(position.z);
    out->push_back(r);
    out->push_back(g);
    out->push_back(b);
}

Vec3 TorusPoint(float majorRadius, float minorRadius, float phi, float theta) {
    const float cPhi = std::cos(phi);
    const float sPhi = std::sin(phi);
    const float cTheta = std::cos(theta);
    const float sTheta = std::sin(theta);
    const float ringRadius = majorRadius + minorRadius * cTheta;
    return Vec3(ringRadius * cPhi, ringRadius * sPhi, minorRadius * sTheta);
}

Vec3 ClampToTorusVolume(const Vec3& position, float majorRadius, float minorRadius, float padding_m) {
    const float safeMinorRadius = std::max(1.0e-4f, minorRadius - std::max(0.0f, padding_m));
    const float radialXY = std::sqrt((position.x * position.x) + (position.y * position.y));

    Vec3 coreCenter(majorRadius, 0.0f, 0.0f);
    if (radialXY > 1.0e-6f) {
        coreCenter = Vec3(
            majorRadius * (position.x / radialXY),
            majorRadius * (position.y / radialXY),
            0.0f);
    }

    Vec3 displacement = position - coreCenter;
    const float tubeRadius = displacement.Magnitude();
    if (tubeRadius <= safeMinorRadius) {
        return position;
    }

    if (tubeRadius <= 1.0e-6f) {
        displacement = Vec3(safeMinorRadius, 0.0f, 0.0f);
    } else {
        displacement = displacement * (safeMinorRadius / tubeRadius);
    }
    return coreCenter + displacement;
}

bool IsInsideTorusVolume(const Vec3& position, float majorRadius, float minorRadius, float padding_m) {
    const float safeMinorRadius = std::max(0.0f, minorRadius - std::max(0.0f, padding_m));
    const float radialXY = std::sqrt((position.x * position.x) + (position.y * position.y));
    const float tubeRadius =
        std::sqrt(((radialXY - majorRadius) * (radialXY - majorRadius)) + (position.z * position.z));
    return tubeRadius <= safeMinorRadius;
}

bool InterpolationSegmentStaysInsideTorus(
    const Vec3& start,
    const Vec3& end,
    float majorRadius,
    float minorRadius,
    float padding_m) {
    if (!IsInsideTorusVolume(start, majorRadius, minorRadius, padding_m) ||
        !IsInsideTorusVolume(end, majorRadius, minorRadius, padding_m)) {
        return false;
    }

    constexpr float kSamples[] = {0.25f, 0.5f, 0.75f};
    for (const float alpha : kSamples) {
        const Vec3 sample = start + ((end - start) * alpha);
        if (!IsInsideTorusVolume(sample, majorRadius, minorRadius, padding_m)) {
            return false;
        }
    }
    return true;
}

void PushLine(
    std::vector<float>* out,
    const Vec3& a,
    const Vec3& b,
    float r,
    float g,
    float bl) {
    PushVertex(out, a, r, g, bl);
    PushVertex(out, b, r, g, bl);
}

float Clamp01(float value) {
    return std::max(0.0f, std::min(value, 1.0f));
}

float Lerp(float a, float b, float t) {
    return a + ((b - a) * Clamp01(t));
}

Vec3 SceneFogColor(float ignitionIntensity) {
    const float t = Clamp01(ignitionIntensity);
    return Vec3(
        0.05f + (0.10f * t),
        0.07f + (0.06f * t),
        0.10f + (0.15f * t));
}

constexpr std::size_t kParticleVertexStrideFloats = 9;

void SpeciesColor(ReplaySpecies species, float* r, float* g, float* b) {
    switch (species) {
        case ReplaySpecies::Deuterium:
            *r = 0.20f;
            *g = 0.62f;
            *b = 1.00f;
            return;
        case ReplaySpecies::Tritium:
            *r = 1.00f;
            *g = 0.54f;
            *b = 0.18f;
            return;
        case ReplaySpecies::Helium:
            *r = 0.98f;
            *g = 0.93f;
            *b = 0.26f;
            return;
        case ReplaySpecies::Unknown:
            *r = 0.62f;
            *g = 0.62f;
            *b = 0.62f;
            return;
    }
}

float NormalizeLinear(double value, const ParticleViewRange& range) {
    if (!range.valid || !std::isfinite(value)) {
        return 0.0f;
    }
    if (range.maxValue <= range.minValue) {
        return (range.maxValue > 0.0) ? 1.0f : 0.0f;
    }
    return Clamp01(static_cast<float>((value - range.minValue) / (range.maxValue - range.minValue)));
}

float NormalizeLogarithmic(double value, const ParticleViewRange& range) {
    if (!range.valid || !std::isfinite(value) || value <= 0.0) {
        return 0.0f;
    }
    const double minValue = std::max(range.minValue, 1.0e-30);
    const double maxValue = std::max(range.maxValue, minValue);
    const double logValue = std::log10(std::max(value, minValue));
    const double logMin = std::log10(minValue);
    const double logMax = std::log10(maxValue);
    if (logMax <= logMin) {
        return 1.0f;
    }
    return Clamp01(static_cast<float>((logValue - logMin) / (logMax - logMin)));
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
            *r = Lerp(kStops[i - 1].r, kStops[i].r, localT);
            *g = Lerp(kStops[i - 1].g, kStops[i].g, localT);
            *b = Lerp(kStops[i - 1].b, kStops[i].b, localT);
            return;
        }
    }

    *r = kStops[3].r;
    *g = kStops[3].g;
    *b = kStops[3].b;
}

float WrapAngleDegrees(float angleDeg) {
    float wrapped = std::fmod(angleDeg, 360.0f);
    if (wrapped < 0.0f) {
        wrapped += 360.0f;
    }
    return wrapped;
}

float AngularDistanceDegrees(float aDeg, float bDeg) {
    float diff = std::fabs(WrapAngleDegrees(aDeg) - WrapAngleDegrees(bDeg));
    if (diff > 180.0f) {
        diff = 360.0f - diff;
    }
    return diff;
}

float MinorRadiusMeters(const ReplayParticle& particle, float majorRadius_m) {
    const float radialXY = std::sqrt((particle.position_m.x * particle.position_m.x) + (particle.position_m.y * particle.position_m.y));
    return std::sqrt(((radialXY - majorRadius_m) * (radialXY - majorRadius_m)) + (particle.position_m.z * particle.position_m.z));
}

float ToroidalAngleDegrees(const ReplayParticle& particle) {
    return WrapAngleDegrees(std::atan2(particle.position_m.y, particle.position_m.x) * (180.0f / 3.14159265359f));
}

float PoloidalAngleDegrees(const ReplayParticle& particle, float majorRadius_m) {
    const float radialXY = std::sqrt((particle.position_m.x * particle.position_m.x) + (particle.position_m.y * particle.position_m.y));
    const float tubeX = radialXY - majorRadius_m;
    return WrapAngleDegrees(std::atan2(particle.position_m.z, tubeX) * (180.0f / 3.14159265359f));
}

bool ParticlePassesSlice(const ReplayParticle& particle, const ParticleViewContext& context) {
    if (!context.slice.enableToroidalSlice && !context.slice.enablePoloidalSlice) {
        return true;
    }

    if (context.slice.enableToroidalSlice) {
        const float delta = AngularDistanceDegrees(ToroidalAngleDegrees(particle), context.slice.toroidalCenter_deg);
        if (delta > std::max(0.0f, context.slice.toroidalHalfWidth_deg)) {
            return false;
        }
    }
    if (context.slice.enablePoloidalSlice) {
        const float delta = AngularDistanceDegrees(PoloidalAngleDegrees(particle, context.majorRadius_m), context.slice.poloidalCenter_deg);
        if (delta > std::max(0.0f, context.slice.poloidalHalfWidth_deg)) {
            return false;
        }
    }
    return true;
}

int FindRadialBinIndex(
    const ReplayParticle& particle,
    const std::vector<ReplayRadialProfileBin>* radialProfile,
    float majorRadius_m,
    float minorRadius_m) {
    if (radialProfile == nullptr || radialProfile->empty() || minorRadius_m <= 0.0f) {
        return -1;
    }
    const float minorRadius = MinorRadiusMeters(particle, majorRadius_m);
    const float normalized = Clamp01(minorRadius / std::max(minorRadius_m, 1.0e-5f));
    int bin = static_cast<int>(normalized * static_cast<float>(radialProfile->size()));
    if (bin >= static_cast<int>(radialProfile->size())) {
        bin = static_cast<int>(radialProfile->size()) - 1;
    }
    return std::max(bin, 0);
}

double FusionRateForParticle(
    const ReplayParticle& particle,
    const std::vector<ReplayRadialProfileBin>* radialProfile,
    float majorRadius_m,
    float minorRadius_m) {
    const int binIndex = FindRadialBinIndex(particle, radialProfile, majorRadius_m, minorRadius_m);
    if (binIndex < 0 || radialProfile == nullptr || binIndex >= static_cast<int>(radialProfile->size())) {
        return 0.0;
    }
    const ReplayRadialProfileBin& bin = (*radialProfile)[static_cast<std::size_t>(binIndex)];
    return bin.fusionRatePlaceholder ? 0.0 : bin.fusionRate_m3_s;
}

double WallRiskForParticle(const ReplayParticle& particle, const ReplayWallInteractionPoint* wallInteraction, float majorRadius_m, float minorRadius_m) {
    if (minorRadius_m <= 0.0f) {
        return 0.0;
    }
    const float minorRadius = MinorRadiusMeters(particle, majorRadius_m);
    const double normalizedRadius = Clamp01(minorRadius / minorRadius_m);
    double activity = 0.30;
    if (wallInteraction != nullptr) {
        const double hitSignal = std::log10(1.0 + static_cast<double>(wallInteraction->wallHitCountStep));
        const double lossSignal = std::log10(1.0 + std::max(0.0, wallInteraction->wallLossWeightStep) * 1.0e6);
        activity = std::min(1.0, 0.25 + (0.45 * hitSignal) + (0.30 * lossSignal));
    }
    return normalizedRadius * activity;
}

double ParticleValueForContext(const ReplayParticle& particle, const ParticleViewContext& context) {
    switch (context.mode) {
        case ParticleViewMode::Species:
            return 0.0;
        case ParticleViewMode::Energy:
            return std::max(0.0, particle.kineticEnergy_keV);
        case ParticleViewMode::PitchAngle:
            return std::isfinite(particle.pitchAngle_deg) ? std::max(0.0, std::min(180.0, particle.pitchAngle_deg)) : 0.0;
        case ParticleViewMode::MagneticMagnitude:
            return std::max(0.0, particle.magneticMagnitude_T);
        case ParticleViewMode::ElectricMagnitude:
            return std::max(0.0, particle.electricMagnitude_VPerM);
        case ParticleViewMode::LorentzAcceleration:
            return std::max(0.0, particle.lorentzAccelerationMagnitude_mPerS2);
        case ParticleViewMode::FusionRate: {
            const double current = FusionRateForParticle(particle, context.radialProfile, context.majorRadius_m, context.minorRadius_m);
            if (!context.diffAgainstCompare) {
                return current;
            }
            const double compare = FusionRateForParticle(particle, context.compareRadialProfile, context.majorRadius_m, context.minorRadius_m);
            return current - compare;
        }
        case ParticleViewMode::WallLossRisk: {
            const double current = WallRiskForParticle(particle, context.wallInteraction, context.majorRadius_m, context.minorRadius_m);
            if (!context.diffAgainstCompare) {
                return current;
            }
            const double compare = WallRiskForParticle(particle, context.compareWallInteraction, context.majorRadius_m, context.minorRadius_m);
            return current - compare;
        }
    }
    return 0.0;
}

bool UsesLogarithmicScale(ParticleViewMode mode) {
    switch (mode) {
        case ParticleViewMode::Energy:
        case ParticleViewMode::MagneticMagnitude:
        case ParticleViewMode::ElectricMagnitude:
        case ParticleViewMode::LorentzAcceleration:
        case ParticleViewMode::FusionRate:
            return true;
        case ParticleViewMode::Species:
        case ParticleViewMode::PitchAngle:
        case ParticleViewMode::WallLossRisk:
            return false;
    }
    return false;
}

void DiffPalette(float value, const ParticleViewRange& range, float* r, float* g, float* b) {
    if (!range.valid) {
        *r = 0.5f;
        *g = 0.5f;
        *b = 0.5f;
        return;
    }
    const double amplitude = std::max(std::fabs(range.minValue), std::fabs(range.maxValue));
    if (amplitude <= 1.0e-30) {
        *r = 0.5f;
        *g = 0.5f;
        *b = 0.5f;
        return;
    }
    const float t = Clamp01(static_cast<float>((value / amplitude) * 0.5 + 0.5));
    struct Stop {
        float t;
        float r;
        float g;
        float b;
    };
    constexpr Stop kStops[] = {
        {0.0f, 0.12f, 0.32f, 0.84f},
        {0.5f, 0.70f, 0.72f, 0.74f},
        {1.0f, 0.96f, 0.28f, 0.14f},
    };
    if (t <= 0.5f) {
        const float local = t / 0.5f;
        *r = Lerp(kStops[0].r, kStops[1].r, local);
        *g = Lerp(kStops[0].g, kStops[1].g, local);
        *b = Lerp(kStops[0].b, kStops[1].b, local);
        return;
    }
    const float local = (t - 0.5f) / 0.5f;
    *r = Lerp(kStops[1].r, kStops[2].r, local);
    *g = Lerp(kStops[1].g, kStops[2].g, local);
    *b = Lerp(kStops[1].b, kStops[2].b, local);
}

Vec3 TorusNormal(float phi, float theta) {
    const float cPhi = std::cos(phi);
    const float sPhi = std::sin(phi);
    const float cTheta = std::cos(theta);
    const float sTheta = std::sin(theta);
    return Vec3(cTheta * cPhi, cTheta * sPhi, sTheta).Normalized();
}

void PushShellVertex(std::vector<float>* out, const Vec3& position, const Vec3& normal, float r, float g, float b) {
    out->push_back(position.x);
    out->push_back(position.y);
    out->push_back(position.z);
    out->push_back(normal.x);
    out->push_back(normal.y);
    out->push_back(normal.z);
    out->push_back(r);
    out->push_back(g);
    out->push_back(b);
}

void PushShellTriangle(
    std::vector<float>* out,
    const Vec3& a,
    const Vec3& b,
    const Vec3& c,
    const Vec3& na,
    const Vec3& nb,
    const Vec3& nc,
    float r,
    float g,
    float bl) {
    PushShellVertex(out, a, na, r, g, bl);
    PushShellVertex(out, b, nb, r, g, bl);
    PushShellVertex(out, c, nc, r, g, bl);
}

}  // namespace

const char* ParticleViewModeName(ParticleViewMode mode) {
    switch (mode) {
        case ParticleViewMode::Species:
            return "Species";
        case ParticleViewMode::Energy:
            return "Energy";
        case ParticleViewMode::PitchAngle:
            return "Pitch Angle";
        case ParticleViewMode::MagneticMagnitude:
            return "|B|";
        case ParticleViewMode::ElectricMagnitude:
            return "|E|";
        case ParticleViewMode::LorentzAcceleration:
            return "Lorentz Accel";
        case ParticleViewMode::FusionRate:
            return "Fusion Rate";
        case ParticleViewMode::WallLossRisk:
            return "Wall Loss Risk";
    }
    return "Unknown";
}

ParticleViewRange ComputeParticleViewRange(const ReplayFrame& frame, const ParticleViewContext& context) {
    ParticleViewRange range;
    if (context.mode == ParticleViewMode::Species) {
        return range;
    }
    if (context.mode == ParticleViewMode::PitchAngle) {
        range.minValue = 0.0;
        range.maxValue = 180.0;
        range.valid = true;
        return range;
    }

    double minValue = std::numeric_limits<double>::infinity();
    double maxValue = -std::numeric_limits<double>::infinity();
    for (const ReplayParticle& particle : frame.particles) {
        if (!ParticlePassesSlice(particle, context) && context.slice.dimOutsideSlice) {
            continue;
        }
        const double value = ParticleValueForContext(particle, context);
        if (!std::isfinite(value)) {
            continue;
        }
        if (!context.diffAgainstCompare && value < 0.0) {
            continue;
        }
        if (UsesLogarithmicScale(context.mode) && value <= 0.0) {
            continue;
        }
        minValue = std::min(minValue, value);
        maxValue = std::max(maxValue, value);
    }

    if (!std::isfinite(minValue)) {
        return range;
    }
    if (maxValue <= minValue) {
        if (context.diffAgainstCompare) {
            maxValue = std::fabs(maxValue);
            minValue = -maxValue;
        } else {
            minValue = 0.0;
        }
    }
    range.minValue = minValue;
    range.maxValue = maxValue;
    range.valid = true;
    return range;
}

void ComputeParticleViewColor(
    const ReplayParticle& particle,
    const ParticleViewContext& context,
    float* r,
    float* g,
    float* b) {
    if (context.mode == ParticleViewMode::Species) {
        SpeciesColor(particle.species, r, g, b);
    } else {
        const double value = ParticleValueForContext(particle, context);
        if (!std::isfinite(value)) {
            *r = 0.45f;
            *g = 0.45f;
            *b = 0.45f;
        } else if (context.diffAgainstCompare &&
                   (context.mode == ParticleViewMode::FusionRate || context.mode == ParticleViewMode::WallLossRisk)) {
            DiffPalette(value, context.range, r, g, b);
        } else {
            const float normalized = UsesLogarithmicScale(context.mode)
                ? NormalizeLogarithmic(value, context.range)
                : NormalizeLinear(value, context.range);
            ContinuousPalette(normalized, r, g, b);
        }
    }

    if (!ParticlePassesSlice(particle, context) && context.slice.dimOutsideSlice) {
        *r = Lerp(*r, 0.10f, 0.80f);
        *g = Lerp(*g, 0.10f, 0.80f);
        *b = Lerp(*b, 0.12f, 0.80f);
    }
}

bool GlRenderer::Initialize(std::size_t maxParticles, std::string* errorOut) {
    maxParticles_ = std::max<std::size_t>(1, maxParticles);
    particleVertexBuffer_.reserve(maxParticles_ * kParticleVertexStrideFloats);

    if (!InitializeParticlePipeline(errorOut)) {
        Shutdown();
        return false;
    }

    if (!InitializeLinePipeline(errorOut)) {
        Shutdown();
        return false;
    }

    if (!InitializeShellPipeline(errorOut)) {
        Shutdown();
        return false;
    }

    RebuildSceneGeometry(2.0f, 0.5f);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    return true;
}

bool GlRenderer::InitializeParticlePipeline(std::string* errorOut) {
    pointProgram_ = CreateProgram(kPointVertexShader, kPointFragmentShader, errorOut);
    if (pointProgram_ == 0) {
        return false;
    }

    pointViewProjectionLocation_ = glGetUniformLocation(pointProgram_, "uViewProjection");
    pointSizeLocation_ = glGetUniformLocation(pointProgram_, "uPointSize");
    pointInterpolationAlphaLocation_ = glGetUniformLocation(pointProgram_, "uInterpolationAlpha");
    pointFogColorLocation_ = glGetUniformLocation(pointProgram_, "uFogColor");
    pointAlphaScaleLocation_ = glGetUniformLocation(pointProgram_, "uAlphaScale");
    pointColorBoostLocation_ = glGetUniformLocation(pointProgram_, "uColorBoost");
    pointFogMixLocation_ = glGetUniformLocation(pointProgram_, "uFogMix");

    glGenVertexArrays(1, &pointVao_);
    glGenBuffers(1, &pointVbo_);

    glBindVertexArray(pointVao_);
    glBindBuffer(GL_ARRAY_BUFFER, pointVbo_);
    glBufferData(
        GL_ARRAY_BUFFER,
        static_cast<GLsizeiptr>(maxParticles_ * kParticleVertexStrideFloats * sizeof(float)),
        nullptr,
        GL_STREAM_DRAW);

    constexpr GLsizei stride = static_cast<GLsizei>(kParticleVertexStrideFloats * sizeof(float));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(0));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(3 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(6 * sizeof(float)));

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    return true;
}

bool GlRenderer::InitializeLinePipeline(std::string* errorOut) {
    lineProgram_ = CreateProgram(kLineVertexShader, kLineFragmentShader, errorOut);
    if (lineProgram_ == 0) {
        return false;
    }

    lineViewProjectionLocation_ = glGetUniformLocation(lineProgram_, "uViewProjection");
    lineFogColorLocation_ = glGetUniformLocation(lineProgram_, "uFogColor");
    lineAlphaScaleLocation_ = glGetUniformLocation(lineProgram_, "uAlphaScale");
    lineFogMixLocation_ = glGetUniformLocation(lineProgram_, "uFogMix");

    glGenVertexArrays(1, &lineVao_);
    glGenBuffers(1, &lineVbo_);

    glBindVertexArray(lineVao_);
    glBindBuffer(GL_ARRAY_BUFFER, lineVbo_);
    glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STATIC_DRAW);

    constexpr GLsizei stride = static_cast<GLsizei>(6 * sizeof(float));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(0));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(3 * sizeof(float)));

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    glGenVertexArrays(1, &overlayLineVao_);
    glGenBuffers(1, &overlayLineVbo_);
    glBindVertexArray(overlayLineVao_);
    glBindBuffer(GL_ARRAY_BUFFER, overlayLineVbo_);
    glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STREAM_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(0));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(3 * sizeof(float)));
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    return true;
}

bool GlRenderer::InitializeShellPipeline(std::string* errorOut) {
    shellProgram_ = CreateProgram(kShellVertexShader, kShellFragmentShader, errorOut);
    if (shellProgram_ == 0) {
        return false;
    }

    shellViewProjectionLocation_ = glGetUniformLocation(shellProgram_, "uViewProjection");
    shellFogColorLocation_ = glGetUniformLocation(shellProgram_, "uFogColor");
    shellAlphaScaleLocation_ = glGetUniformLocation(shellProgram_, "uAlphaScale");
    shellFogMixLocation_ = glGetUniformLocation(shellProgram_, "uFogMix");

    glGenVertexArrays(1, &shellVao_);
    glGenBuffers(1, &shellVbo_);
    glBindVertexArray(shellVao_);
    glBindBuffer(GL_ARRAY_BUFFER, shellVbo_);
    glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STATIC_DRAW);

    constexpr GLsizei stride = static_cast<GLsizei>(9 * sizeof(float));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(0));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(3 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(6 * sizeof(float)));
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    glGenVertexArrays(1, &trailVao_);
    glGenBuffers(1, &trailVbo_);
    glBindVertexArray(trailVao_);
    glBindBuffer(GL_ARRAY_BUFFER, trailVbo_);
    glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STREAM_DRAW);
    constexpr GLsizei lineStride = static_cast<GLsizei>(6 * sizeof(float));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, lineStride, reinterpret_cast<void*>(0));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, lineStride, reinterpret_cast<void*>(3 * sizeof(float)));
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    return true;
}

void GlRenderer::Shutdown() {
    if (pointVbo_ != 0) {
        glDeleteBuffers(1, &pointVbo_);
        pointVbo_ = 0;
    }
    if (pointVao_ != 0) {
        glDeleteVertexArrays(1, &pointVao_);
        pointVao_ = 0;
    }
    if (pointProgram_ != 0) {
        glDeleteProgram(pointProgram_);
        pointProgram_ = 0;
    }

    if (lineVbo_ != 0) {
        glDeleteBuffers(1, &lineVbo_);
        lineVbo_ = 0;
    }
    if (lineVao_ != 0) {
        glDeleteVertexArrays(1, &lineVao_);
        lineVao_ = 0;
    }
    if (lineProgram_ != 0) {
        glDeleteProgram(lineProgram_);
        lineProgram_ = 0;
    }
    if (overlayLineVbo_ != 0) {
        glDeleteBuffers(1, &overlayLineVbo_);
        overlayLineVbo_ = 0;
    }
    if (overlayLineVao_ != 0) {
        glDeleteVertexArrays(1, &overlayLineVao_);
        overlayLineVao_ = 0;
    }
    if (shellVbo_ != 0) {
        glDeleteBuffers(1, &shellVbo_);
        shellVbo_ = 0;
    }
    if (shellVao_ != 0) {
        glDeleteVertexArrays(1, &shellVao_);
        shellVao_ = 0;
    }
    if (shellProgram_ != 0) {
        glDeleteProgram(shellProgram_);
        shellProgram_ = 0;
    }
    if (trailVbo_ != 0) {
        glDeleteBuffers(1, &trailVbo_);
        trailVbo_ = 0;
    }
    if (trailVao_ != 0) {
        glDeleteVertexArrays(1, &trailVao_);
        trailVao_ = 0;
    }

    uploadedParticles_ = 0;
    lineVertexCount_ = 0;
    overlayLineVertexCount_ = 0;
    shellVertexCount_ = 0;
    trailVertexCount_ = 0;
    pointInterpolationAlphaLocation_ = -1;
    pointFogColorLocation_ = -1;
    pointAlphaScaleLocation_ = -1;
    pointColorBoostLocation_ = -1;
    pointFogMixLocation_ = -1;
    lineFogColorLocation_ = -1;
    lineAlphaScaleLocation_ = -1;
    lineFogMixLocation_ = -1;
    shellAlphaScaleLocation_ = -1;
    shellFogMixLocation_ = -1;
    particleVertexBuffer_.clear();
    trailVertexBuffer_.clear();
    overlayLineVertexBuffer_.clear();
}

void GlRenderer::SetTorusGeometry(float majorRadius_m, float minorRadius_m) {
    torusMajorRadius_m_ = majorRadius_m;
    torusMinorRadius_m_ = minorRadius_m;
    RebuildSceneGeometry(majorRadius_m, minorRadius_m);
}

void GlRenderer::EnsureParticleCapacity(std::size_t requiredParticles) {
    if (requiredParticles <= maxParticles_) {
        return;
    }

    maxParticles_ = std::max(requiredParticles, std::max<std::size_t>(maxParticles_ * 2, 1));
    particleVertexBuffer_.reserve(maxParticles_ * kParticleVertexStrideFloats);

    glBindBuffer(GL_ARRAY_BUFFER, pointVbo_);
    glBufferData(
        GL_ARRAY_BUFFER,
        static_cast<GLsizeiptr>(maxParticles_ * kParticleVertexStrideFloats * sizeof(float)),
        nullptr,
        GL_STREAM_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void GlRenderer::RebuildSceneGeometry(float majorRadius_m, float minorRadius_m) {
    constexpr int kMajorSegments = 96;
    constexpr int kMinorSegments = 32;
    constexpr float kTwoPi = 2.0f * 3.14159265359f;

    std::vector<float> lineVertices;
    lineVertices.reserve(static_cast<std::size_t>(kMajorSegments * kMinorSegments * 12));
    std::vector<float> shellVertices;
    shellVertices.reserve(static_cast<std::size_t>(kMajorSegments * kMinorSegments * 54));

    for (int i = 0; i < kMajorSegments; ++i) {
        const float phi0 = kTwoPi * static_cast<float>(i) / static_cast<float>(kMajorSegments);
        const float phi1 = kTwoPi * static_cast<float>(i + 1) / static_cast<float>(kMajorSegments);

        for (int j = 0; j < kMinorSegments; ++j) {
            const float theta0 = kTwoPi * static_cast<float>(j) / static_cast<float>(kMinorSegments);
            const float theta1 = kTwoPi * static_cast<float>(j + 1) / static_cast<float>(kMinorSegments);

            const Vec3 p00 = TorusPoint(majorRadius_m, minorRadius_m, phi0, theta0);
            const Vec3 p10 = TorusPoint(majorRadius_m, minorRadius_m, phi1, theta0);
            const Vec3 p01 = TorusPoint(majorRadius_m, minorRadius_m, phi0, theta1);
            const Vec3 p11 = TorusPoint(majorRadius_m, minorRadius_m, phi1, theta1);
            const Vec3 n00 = TorusNormal(phi0, theta0);
            const Vec3 n10 = TorusNormal(phi1, theta0);
            const Vec3 n01 = TorusNormal(phi0, theta1);
            const Vec3 n11 = TorusNormal(phi1, theta1);

            PushLine(&lineVertices, p00, p10, 0.22f, 0.42f, 0.44f);
            PushLine(&lineVertices, p00, p01, 0.16f, 0.28f, 0.30f);

            PushShellTriangle(&shellVertices, p00, p10, p11, n00, n10, n11, 0.22f, 0.52f, 0.58f);
            PushShellTriangle(&shellVertices, p00, p11, p01, n00, n11, n01, 0.20f, 0.46f, 0.54f);
        }
    }

    constexpr int kCoilCount = 18;
    constexpr int kCoilSegments = 42;
    const float coilRadius = minorRadius_m * 1.28f;
    for (int coil = 0; coil < kCoilCount; ++coil) {
        const float phi = kTwoPi * static_cast<float>(coil) / static_cast<float>(kCoilCount);
        const Vec3 center(majorRadius_m * std::cos(phi), majorRadius_m * std::sin(phi), 0.0f);
        const Vec3 radial(std::cos(phi), std::sin(phi), 0.0f);
        const Vec3 up(0.0f, 0.0f, 1.0f);
        for (int seg = 0; seg < kCoilSegments; ++seg) {
            const float a0 = kTwoPi * static_cast<float>(seg) / static_cast<float>(kCoilSegments);
            const float a1 = kTwoPi * static_cast<float>(seg + 1) / static_cast<float>(kCoilSegments);
            const Vec3 p0 = center + (radial * (coilRadius * std::cos(a0))) + (up * (coilRadius * std::sin(a0)));
            const Vec3 p1 = center + (radial * (coilRadius * std::cos(a1))) + (up * (coilRadius * std::sin(a1)));
            PushLine(&lineVertices, p0, p1, 0.42f, 0.38f, 0.28f);
        }
    }

    const Vec3 injectorStart(majorRadius_m + (minorRadius_m * 1.65f), 0.0f, minorRadius_m * 0.10f);
    const Vec3 injectorEnd(majorRadius_m + (minorRadius_m * 0.40f), 0.0f, 0.0f);
    PushLine(&lineVertices, injectorStart, injectorEnd, 0.86f, 0.70f, 0.30f);
    PushLine(&lineVertices, injectorStart + Vec3(0.0f, minorRadius_m * 0.10f, 0.0f), injectorStart - Vec3(0.0f, minorRadius_m * 0.10f, 0.0f), 0.86f, 0.70f, 0.30f);
    PushLine(&lineVertices, injectorStart + Vec3(0.0f, 0.0f, minorRadius_m * 0.10f), injectorStart - Vec3(0.0f, 0.0f, minorRadius_m * 0.10f), 0.86f, 0.70f, 0.30f);

    lineVertexCount_ = static_cast<int>(lineVertices.size() / 6);
    shellVertexCount_ = static_cast<int>(shellVertices.size() / 9);

    glBindBuffer(GL_ARRAY_BUFFER, lineVbo_);
    glBufferData(
        GL_ARRAY_BUFFER,
        static_cast<GLsizeiptr>(lineVertices.size() * sizeof(float)),
        lineVertices.data(),
        GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindBuffer(GL_ARRAY_BUFFER, shellVbo_);
    glBufferData(
        GL_ARRAY_BUFFER,
        static_cast<GLsizeiptr>(shellVertices.size() * sizeof(float)),
        shellVertices.data(),
        GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void GlRenderer::UploadFrame(
    const ReplayFrame& frame,
    const ReplayFrame* nextFrame,
    const ParticleViewContext& viewContext) {
    constexpr float kRenderInteriorPadding_m = 0.003f;
    const std::size_t count = frame.particles.size();
    EnsureParticleCapacity(std::max<std::size_t>(count, 1));
    uploadedParticles_ = count;
    particleVertexBuffer_.resize(count * kParticleVertexStrideFloats);
    trailVertexBuffer_.clear();
    trailVertexBuffer_.reserve(count * 12);

    std::unordered_map<uint64_t, Vec3> nextPositionsById;
    if (nextFrame != nullptr) {
        nextPositionsById.reserve(nextFrame->particles.size());
        for (const ReplayParticle& nextParticle : nextFrame->particles) {
            nextPositionsById.emplace(nextParticle.particleIndex, nextParticle.position_m);
        }
    }

    for (std::size_t i = 0; i < count; ++i) {
        const ReplayParticle& particle = frame.particles[i];
        float r = 0.0f;
        float g = 0.0f;
        float b = 0.0f;
        ComputeParticleViewColor(particle, viewContext, &r, &g, &b);

        // A torus is non-convex, so a straight Cartesian blend between two valid
        // in-vessel positions can still cut outside the chamber. When that would
        // happen, keep the particle on its current snapshot until the next frame.
        const Vec3 currentPosition = ClampToTorusVolume(
            particle.position_m,
            torusMajorRadius_m_,
            torusMinorRadius_m_,
            kRenderInteriorPadding_m);
        Vec3 positionDelta(0.0f, 0.0f, 0.0f);
        if (!nextPositionsById.empty()) {
            const auto it = nextPositionsById.find(particle.particleIndex);
            if (it != nextPositionsById.end()) {
                const Vec3 nextPosition = ClampToTorusVolume(
                    it->second,
                    torusMajorRadius_m_,
                    torusMinorRadius_m_,
                    kRenderInteriorPadding_m);
                if (InterpolationSegmentStaysInsideTorus(
                        currentPosition,
                        nextPosition,
                        torusMajorRadius_m_,
                        torusMinorRadius_m_,
                        kRenderInteriorPadding_m)) {
                    positionDelta = nextPosition - currentPosition;
                }
            }
        }

        Vec3 trailStart = currentPosition;
        float trailFade = 0.0f;
        const float deltaMagnitude = positionDelta.Magnitude();
        if (deltaMagnitude > 1.0e-5f) {
            trailStart = ClampToTorusVolume(
                currentPosition - (positionDelta * 0.58f),
                torusMajorRadius_m_,
                torusMinorRadius_m_,
                kRenderInteriorPadding_m);
            trailFade = Clamp01(deltaMagnitude / std::max(0.08f * torusMinorRadius_m_, 1.0e-4f));
        } else if (particle.velocity_mPerS.Magnitude() > 1.0e-5f) {
            const Vec3 fallbackTrail = particle.velocity_mPerS.Normalized() * std::max(0.04f, torusMinorRadius_m_ * 0.06f);
            trailStart = ClampToTorusVolume(
                currentPosition - fallbackTrail,
                torusMajorRadius_m_,
                torusMinorRadius_m_,
                kRenderInteriorPadding_m);
            trailFade = 0.35f;
        }
        if ((trailStart - currentPosition).Magnitude() > 1.0e-4f) {
            PushVertex(&trailVertexBuffer_, trailStart, r * (0.14f + (0.12f * trailFade)), g * (0.14f + (0.12f * trailFade)), b * (0.16f + (0.12f * trailFade)));
            PushVertex(&trailVertexBuffer_, currentPosition, r * (0.46f + (0.20f * trailFade)), g * (0.46f + (0.20f * trailFade)), b * (0.50f + (0.20f * trailFade)));
        }

        const std::size_t base = i * kParticleVertexStrideFloats;
        particleVertexBuffer_[base + 0] = currentPosition.x;
        particleVertexBuffer_[base + 1] = currentPosition.y;
        particleVertexBuffer_[base + 2] = currentPosition.z;
        particleVertexBuffer_[base + 3] = positionDelta.x;
        particleVertexBuffer_[base + 4] = positionDelta.y;
        particleVertexBuffer_[base + 5] = positionDelta.z;
        particleVertexBuffer_[base + 6] = r;
        particleVertexBuffer_[base + 7] = g;
        particleVertexBuffer_[base + 8] = b;
    }

    if (uploadedParticles_ == 0) {
        return;
    }

    glBindBuffer(GL_ARRAY_BUFFER, pointVbo_);
    glBufferData(
        GL_ARRAY_BUFFER,
        static_cast<GLsizeiptr>(particleVertexBuffer_.size() * sizeof(float)),
        particleVertexBuffer_.data(),
        GL_STREAM_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    trailVertexCount_ = static_cast<int>(trailVertexBuffer_.size() / 6);
    glBindBuffer(GL_ARRAY_BUFFER, trailVbo_);
    glBufferData(
        GL_ARRAY_BUFFER,
        static_cast<GLsizeiptr>(trailVertexBuffer_.size() * sizeof(float)),
        trailVertexBuffer_.empty() ? nullptr : trailVertexBuffer_.data(),
        GL_STREAM_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void GlRenderer::UploadOverlayLines(const std::vector<float>& overlayVertices) {
    overlayLineVertexBuffer_ = overlayVertices;
    overlayLineVertexCount_ = static_cast<int>(overlayLineVertexBuffer_.size() / 6);

    glBindBuffer(GL_ARRAY_BUFFER, overlayLineVbo_);
    glBufferData(
        GL_ARRAY_BUFFER,
        static_cast<GLsizeiptr>(overlayLineVertexBuffer_.size() * sizeof(float)),
        overlayLineVertexBuffer_.empty() ? nullptr : overlayLineVertexBuffer_.data(),
        GL_STREAM_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void GlRenderer::Draw(
    const Mat4& viewProjection,
    float pointSizePixels,
    float interpolationAlpha,
    float ignitionIntensity,
    const GraphicsStyleSettings& graphicsStyle) const {
    const Vec3 fogColor = SceneFogColor(ignitionIntensity);
    const float fogStrength = Clamp01(graphicsStyle.fogStrength);

    glDepthMask(GL_FALSE);
    if (graphicsStyle.showVesselShell) {
        glUseProgram(shellProgram_);
        glUniformMatrix4fv(shellViewProjectionLocation_, 1, GL_FALSE, viewProjection.Data());
        glUniform3f(shellFogColorLocation_, fogColor.x, fogColor.y, fogColor.z);
        glUniform1f(shellAlphaScaleLocation_, Clamp01(graphicsStyle.shellOpacity));
        glUniform1f(shellFogMixLocation_, fogStrength * 0.85f);
        glBindVertexArray(shellVao_);
        glDrawArrays(GL_TRIANGLES, 0, shellVertexCount_);
        glBindVertexArray(0);
    }

    glUseProgram(lineProgram_);
    glUniformMatrix4fv(lineViewProjectionLocation_, 1, GL_FALSE, viewProjection.Data());
    glUniform3f(lineFogColorLocation_, fogColor.x, fogColor.y, fogColor.z);
    glUniform1f(lineAlphaScaleLocation_, Clamp01(graphicsStyle.sceneLineOpacity));
    glUniform1f(lineFogMixLocation_, fogStrength);
    glLineWidth(std::max(1.0f, graphicsStyle.sceneLineWidth));
    glBindVertexArray(lineVao_);
    glDrawArrays(GL_LINES, 0, lineVertexCount_);
    glBindVertexArray(0);

    if (graphicsStyle.showParticleTrails && trailVertexCount_ > 0) {
        glUniform1f(lineAlphaScaleLocation_, Clamp01(graphicsStyle.trailOpacity));
        glUniform1f(lineFogMixLocation_, fogStrength * 0.35f);
        glLineWidth(std::max(1.0f, graphicsStyle.trailLineWidth));
        glBindVertexArray(trailVao_);
        glDrawArrays(GL_LINES, 0, trailVertexCount_);
        glBindVertexArray(0);
    }

    glUseProgram(pointProgram_);
    glUniformMatrix4fv(pointViewProjectionLocation_, 1, GL_FALSE, viewProjection.Data());
    glUniform1f(pointInterpolationAlphaLocation_, std::max(0.0f, std::min(1.0f, interpolationAlpha)));
    glUniform3f(pointFogColorLocation_, fogColor.x, fogColor.y, fogColor.z);
    glUniform1f(pointFogMixLocation_, Clamp01(graphicsStyle.pointFogStrength));
    glBindVertexArray(pointVao_);

    if (graphicsStyle.showDensitySplat) {
        glBlendFunc(GL_SRC_ALPHA, GL_ONE);
        glUniform1f(pointSizeLocation_, std::max(1.5f, pointSizePixels * std::max(1.1f, graphicsStyle.densitySplatScale)));
        glUniform1f(pointAlphaScaleLocation_, Clamp01(graphicsStyle.densitySplatOpacity));
        glUniform1f(pointColorBoostLocation_, std::max(0.2f, graphicsStyle.particleColorBoost * 0.52f));
        glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(uploadedParticles_));
    }

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(GL_TRUE);
    glUniform1f(pointSizeLocation_, pointSizePixels);
    glUniform1f(pointAlphaScaleLocation_, 1.00f);
    glUniform1f(pointColorBoostLocation_, std::max(0.5f, graphicsStyle.particleColorBoost));
    glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(uploadedParticles_));
    glBindVertexArray(0);

    if (overlayLineVertexCount_ > 0) {
        glUseProgram(lineProgram_);
        glUniformMatrix4fv(lineViewProjectionLocation_, 1, GL_FALSE, viewProjection.Data());
        glUniform3f(lineFogColorLocation_, fogColor.x, fogColor.y, fogColor.z);
        glUniform1f(lineAlphaScaleLocation_, Clamp01(graphicsStyle.overlayLineOpacity));
        glUniform1f(lineFogMixLocation_, Clamp01(graphicsStyle.overlayFogStrength));
        glLineWidth(std::max(1.0f, graphicsStyle.overlayLineWidth));
        const GLboolean depthWasEnabled = glIsEnabled(GL_DEPTH_TEST);
        if (graphicsStyle.xrayOverlays && depthWasEnabled == GL_TRUE) {
            glDisable(GL_DEPTH_TEST);
        }
        glBindVertexArray(overlayLineVao_);
        glDrawArrays(GL_LINES, 0, overlayLineVertexCount_);
        glBindVertexArray(0);
        if (graphicsStyle.xrayOverlays && depthWasEnabled == GL_TRUE) {
            glEnable(GL_DEPTH_TEST);
        }
    }

    glLineWidth(1.0f);
    glUseProgram(0);
}

}  // namespace tokamak::viewer
