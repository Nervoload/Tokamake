#include "tokamak/viewer/gl_renderer.hpp"

#include <algorithm>
#include <cmath>
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
uniform float uIgnitionIntensity;

out vec3 vColor;

void main() {
    vec3 simulatedPosition = aPosition + (aPositionDelta * uInterpolationAlpha);
    gl_Position = uViewProjection * vec4(simulatedPosition, 1.0);
    gl_PointSize = uPointSize;
    vColor = aColor;
}
)GLSL";

const char* kPointFragmentShader = R"GLSL(
#version 330 core
in vec3 vColor;
out vec4 fragColor;

uniform float uIgnitionIntensity;

void main() {
    vec2 centered = gl_PointCoord * 2.0 - 1.0;
    float radiusSq = dot(centered, centered);
    if (radiusSq > 1.0) {
        discard;
    }
    float glow = exp(-3.4 * radiusSq);
    float alpha = mix(0.10, 0.96, glow);
    vec3 color = vColor * (0.55 + 0.75 * glow);
    vec3 plasmaGlow = vec3(0.72, 0.20, 1.00) * uIgnitionIntensity * (0.30 + 0.70 * glow);
    color += plasmaGlow;
    alpha = min(1.0, alpha + (0.18 * uIgnitionIntensity));
    fragColor = vec4(color, alpha);
}
)GLSL";

const char* kLineVertexShader = R"GLSL(
#version 330 core
layout(location = 0) in vec3 aPosition;
layout(location = 1) in vec3 aColor;

uniform mat4 uViewProjection;
uniform float uIgnitionIntensity;
out vec3 vColor;
out float vIgnitionIntensity;

void main() {
    gl_Position = uViewProjection * vec4(aPosition, 1.0);
    vColor = aColor;
    vIgnitionIntensity = uIgnitionIntensity;
}
)GLSL";

const char* kLineFragmentShader = R"GLSL(
#version 330 core
in vec3 vColor;
in float vIgnitionIntensity;
out vec4 fragColor;

void main() {
    vec3 ignitionColor = mix(vColor, vec3(0.76, 0.24, 1.00), 0.70 * vIgnitionIntensity);
    float alpha = 0.55 + (0.35 * vIgnitionIntensity);
    fragColor = vec4(ignitionColor, alpha);
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

void EnrichParticleColor(const ReplayParticle& particle, float* r, float* g, float* b) {
    SpeciesColor(particle.species, r, g, b);

    const double kineticEnergyKeV = std::max(0.0, particle.kineticEnergy_keV);
    const double speed = std::max(0.0, particle.speed_mPerS);
    const float heat = Clamp01(static_cast<float>(std::log10(1.0 + kineticEnergyKeV) / 3.0));
    const float speedGlow = Clamp01(static_cast<float>(std::log10(1.0 + speed) / 7.0));
    const float pitchAccent = std::isfinite(particle.pitchAngle_deg)
        ? Clamp01(static_cast<float>(std::fabs(90.0 - particle.pitchAngle_deg) / 90.0))
        : 0.0f;

    const float hotR = 1.00f;
    const float hotG = Lerp(0.92f, 0.48f, heat);
    const float hotB = Lerp(0.78f, 0.18f, heat);

    *r = Clamp01(Lerp(*r, hotR, 0.48f * heat) * (0.72f + 0.28f * speedGlow));
    *g = Clamp01(Lerp(*g, hotG, 0.42f * heat) * (0.72f + 0.28f * speedGlow));
    *b = Clamp01(Lerp(*b, hotB, 0.42f * heat) + (0.14f * pitchAccent));
}

}  // namespace

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

    RebuildTorusLines(2.0f, 0.5f);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_PROGRAM_POINT_SIZE);
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
    pointIgnitionIntensityLocation_ = glGetUniformLocation(pointProgram_, "uIgnitionIntensity");

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
    lineIgnitionIntensityLocation_ = glGetUniformLocation(lineProgram_, "uIgnitionIntensity");

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

    uploadedParticles_ = 0;
    lineVertexCount_ = 0;
    pointInterpolationAlphaLocation_ = -1;
    pointIgnitionIntensityLocation_ = -1;
    lineIgnitionIntensityLocation_ = -1;
    particleVertexBuffer_.clear();
}

void GlRenderer::SetTorusGeometry(float majorRadius_m, float minorRadius_m) {
    RebuildTorusLines(majorRadius_m, minorRadius_m);
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

void GlRenderer::RebuildTorusLines(float majorRadius_m, float minorRadius_m) {
    constexpr int kMajorSegments = 96;
    constexpr int kMinorSegments = 32;
    constexpr float kTwoPi = 2.0f * 3.14159265359f;

    std::vector<float> vertices;
    vertices.reserve(static_cast<std::size_t>(kMajorSegments * kMinorSegments * 12));

    for (int i = 0; i < kMajorSegments; ++i) {
        const float phi0 = kTwoPi * static_cast<float>(i) / static_cast<float>(kMajorSegments);
        const float phi1 = kTwoPi * static_cast<float>(i + 1) / static_cast<float>(kMajorSegments);

        for (int j = 0; j < kMinorSegments; ++j) {
            const float theta0 = kTwoPi * static_cast<float>(j) / static_cast<float>(kMinorSegments);
            const float theta1 = kTwoPi * static_cast<float>(j + 1) / static_cast<float>(kMinorSegments);

            const Vec3 p00 = TorusPoint(majorRadius_m, minorRadius_m, phi0, theta0);
            const Vec3 p10 = TorusPoint(majorRadius_m, minorRadius_m, phi1, theta0);
            const Vec3 p01 = TorusPoint(majorRadius_m, minorRadius_m, phi0, theta1);

            PushLine(&vertices, p00, p10, 0.25f, 0.94f, 0.70f);
            PushLine(&vertices, p00, p01, 0.12f, 0.62f, 0.56f);
        }
    }

    lineVertexCount_ = static_cast<int>(vertices.size() / 6);

    glBindBuffer(GL_ARRAY_BUFFER, lineVbo_);
    glBufferData(
        GL_ARRAY_BUFFER,
        static_cast<GLsizeiptr>(vertices.size() * sizeof(float)),
        vertices.data(),
        GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void GlRenderer::UploadFrame(const ReplayFrame& frame, const ReplayFrame* nextFrame) {
    const std::size_t count = frame.particles.size();
    EnsureParticleCapacity(std::max<std::size_t>(count, 1));
    uploadedParticles_ = count;
    particleVertexBuffer_.resize(count * kParticleVertexStrideFloats);

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
        EnrichParticleColor(particle, &r, &g, &b);

        Vec3 positionDelta(0.0f, 0.0f, 0.0f);
        if (!nextPositionsById.empty()) {
            const auto it = nextPositionsById.find(particle.particleIndex);
            if (it != nextPositionsById.end()) {
                positionDelta = it->second - particle.position_m;
            }
        }

        const std::size_t base = i * kParticleVertexStrideFloats;
        particleVertexBuffer_[base + 0] = particle.position_m.x;
        particleVertexBuffer_[base + 1] = particle.position_m.y;
        particleVertexBuffer_[base + 2] = particle.position_m.z;
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
}

void GlRenderer::Draw(
    const Mat4& viewProjection,
    float pointSizePixels,
    float interpolationAlpha,
    float ignitionIntensity) const {
    glUseProgram(lineProgram_);
    glUniformMatrix4fv(lineViewProjectionLocation_, 1, GL_FALSE, viewProjection.Data());
    glUniform1f(lineIgnitionIntensityLocation_, std::max(0.0f, std::min(1.0f, ignitionIntensity)));
    glBindVertexArray(lineVao_);
    glDrawArrays(GL_LINES, 0, lineVertexCount_);
    glBindVertexArray(0);

    glUseProgram(pointProgram_);
    glUniformMatrix4fv(pointViewProjectionLocation_, 1, GL_FALSE, viewProjection.Data());
    glUniform1f(pointSizeLocation_, pointSizePixels);
    glUniform1f(pointInterpolationAlphaLocation_, std::max(0.0f, std::min(1.0f, interpolationAlpha)));
    glUniform1f(pointIgnitionIntensityLocation_, std::max(0.0f, std::min(1.0f, ignitionIntensity)));
    glBindVertexArray(pointVao_);
    glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(uploadedParticles_));
    glBindVertexArray(0);

    glUseProgram(0);
}

}  // namespace tokamak::viewer
