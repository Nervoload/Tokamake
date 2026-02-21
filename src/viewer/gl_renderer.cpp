#include "tokamak/viewer/gl_renderer.hpp"

#include <algorithm>
#include <cmath>
#include <string>
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
layout(location = 1) in vec3 aColor;

uniform mat4 uViewProjection;
uniform float uPointSize;

out vec3 vColor;

void main() {
    gl_Position = uViewProjection * vec4(aPosition, 1.0);
    gl_PointSize = uPointSize;
    vColor = aColor;
}
)GLSL";

const char* kPointFragmentShader = R"GLSL(
#version 330 core
in vec3 vColor;
out vec4 fragColor;

void main() {
    vec2 centered = gl_PointCoord * 2.0 - 1.0;
    if (dot(centered, centered) > 1.0) {
        discard;
    }
    fragColor = vec4(vColor, 0.92);
}
)GLSL";

const char* kLineVertexShader = R"GLSL(
#version 330 core
layout(location = 0) in vec3 aPosition;
layout(location = 1) in vec3 aColor;

uniform mat4 uViewProjection;
out vec3 vColor;

void main() {
    gl_Position = uViewProjection * vec4(aPosition, 1.0);
    vColor = aColor;
}
)GLSL";

const char* kLineFragmentShader = R"GLSL(
#version 330 core
in vec3 vColor;
out vec4 fragColor;

void main() {
    fragColor = vec4(vColor, 1.0);
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

}  // namespace

bool GlRenderer::Initialize(std::size_t maxParticles, std::string* errorOut) {
    maxParticles_ = std::max<std::size_t>(1, maxParticles);

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

    glGenVertexArrays(1, &pointVao_);
    glGenBuffers(1, &pointVbo_);

    glBindVertexArray(pointVao_);
    glBindBuffer(GL_ARRAY_BUFFER, pointVbo_);
    glBufferData(
        GL_ARRAY_BUFFER,
        static_cast<GLsizeiptr>(maxParticles_ * 6 * sizeof(float)),
        nullptr,
        GL_DYNAMIC_DRAW);

    constexpr GLsizei stride = static_cast<GLsizei>(6 * sizeof(float));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(0));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(3 * sizeof(float)));

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

    glGenVertexArrays(1, &lineVao_);
    glGenBuffers(1, &lineVbo_);

    glBindVertexArray(lineVao_);
    glBindBuffer(GL_ARRAY_BUFFER, lineVbo_);
    glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);

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
}

void GlRenderer::SetTorusGeometry(float majorRadius_m, float minorRadius_m) {
    RebuildTorusLines(majorRadius_m, minorRadius_m);
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
        GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void GlRenderer::UploadFrame(const ReplayFrame& frame) {
    const std::size_t count = std::min<std::size_t>(frame.particles.size(), maxParticles_);
    uploadedParticles_ = count;

    std::vector<float> vertices;
    vertices.reserve(count * 6);

    for (std::size_t i = 0; i < count; ++i) {
        const ReplayParticle& particle = frame.particles[i];
        float r = 0.0f;
        float g = 0.0f;
        float b = 0.0f;
        SpeciesColor(particle.species, &r, &g, &b);

        vertices.push_back(particle.position_m.x);
        vertices.push_back(particle.position_m.y);
        vertices.push_back(particle.position_m.z);
        vertices.push_back(r);
        vertices.push_back(g);
        vertices.push_back(b);
    }

    glBindBuffer(GL_ARRAY_BUFFER, pointVbo_);
    glBufferSubData(
        GL_ARRAY_BUFFER,
        0,
        static_cast<GLsizeiptr>(vertices.size() * sizeof(float)),
        vertices.data());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void GlRenderer::Draw(const Mat4& viewProjection, float pointSizePixels) const {
    glUseProgram(lineProgram_);
    glUniformMatrix4fv(lineViewProjectionLocation_, 1, GL_FALSE, viewProjection.Data());
    glBindVertexArray(lineVao_);
    glDrawArrays(GL_LINES, 0, lineVertexCount_);
    glBindVertexArray(0);

    glUseProgram(pointProgram_);
    glUniformMatrix4fv(pointViewProjectionLocation_, 1, GL_FALSE, viewProjection.Data());
    glUniform1f(pointSizeLocation_, pointSizePixels);
    glBindVertexArray(pointVao_);
    glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(uploadedParticles_));
    glBindVertexArray(0);

    glUseProgram(0);
}

}  // namespace tokamak::viewer
