#include "render/scene_renderer.h"

#include <cmath>
#include <iostream>
#include <string>

#include <glad/glad.h>

namespace {

unsigned int CompileShader(unsigned int type, const char* source) {
    const unsigned int shader = glCreateShader(type);
    glShaderSource(shader, 1, &source, nullptr);
    glCompileShader(shader);

    int ok = 0;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &ok);
    if (ok == GL_FALSE) {
        int logLength = 0;
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &logLength);
        std::string log(static_cast<size_t>(logLength > 0 ? logLength : 1), '\0');
        glGetShaderInfoLog(shader, logLength, nullptr, log.data());
        std::cerr << "Shader compile failure: " << log << std::endl;
        glDeleteShader(shader);
        return 0;
    }

    return shader;
}

unsigned int CreateProgram(const char* vertexSource, const char* fragmentSource) {
    const unsigned int vs = CompileShader(GL_VERTEX_SHADER, vertexSource);
    if (vs == 0) {
        return 0;
    }

    const unsigned int fs = CompileShader(GL_FRAGMENT_SHADER, fragmentSource);
    if (fs == 0) {
        glDeleteShader(vs);
        return 0;
    }

    const unsigned int program = glCreateProgram();
    glAttachShader(program, vs);
    glAttachShader(program, fs);
    glLinkProgram(program);

    int ok = 0;
    glGetProgramiv(program, GL_LINK_STATUS, &ok);
    glDeleteShader(vs);
    glDeleteShader(fs);

    if (ok == GL_FALSE) {
        int logLength = 0;
        glGetProgramiv(program, GL_INFO_LOG_LENGTH, &logLength);
        std::string log(static_cast<size_t>(logLength > 0 ? logLength : 1), '\0');
        glGetProgramInfoLog(program, logLength, nullptr, log.data());
        std::cerr << "Program link failure: " << log << std::endl;
        glDeleteProgram(program);
        return 0;
    }

    return program;
}

const char* kSceneVertexShader = R"GLSL(
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

const char* kSceneFragmentShader = R"GLSL(
#version 330 core
in vec3 vColor;
out vec4 fragColor;

void main() {
    fragColor = vec4(vColor, 1.0);
}
)GLSL";

Vec3 TorusPoint(float majorRadius, float minorRadius, float phi, float theta) {
    const float cPhi = std::cos(phi);
    const float sPhi = std::sin(phi);
    const float cTheta = std::cos(theta);
    const float sTheta = std::sin(theta);

    const float ring = majorRadius + (minorRadius * cTheta);
    return Vec3(ring * cPhi, ring * sPhi, minorRadius * sTheta);
}

void PushVertex(std::vector<float>& out, const Vec3& position, float r, float g, float b) {
    out.push_back(position.x);
    out.push_back(position.y);
    out.push_back(position.z);
    out.push_back(r);
    out.push_back(g);
    out.push_back(b);
}

void PushLine(
    std::vector<float>& out,
    const Vec3& a,
    const Vec3& b,
    float r,
    float g,
    float bl) {
    PushVertex(out, a, r, g, bl);
    PushVertex(out, b, r, g, bl);
}

}  // namespace

bool SceneRenderer::Initialize() {
    program_ = CreateProgram(kSceneVertexShader, kSceneFragmentShader);
    if (program_ == 0) {
        return false;
    }

    viewProjectionLoc_ = glGetUniformLocation(program_, "uViewProjection");

    glGenVertexArrays(1, &vao_);
    glGenBuffers(1, &vbo_);

    glBindVertexArray(vao_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
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

void SceneRenderer::Shutdown() {
    if (vbo_ != 0) {
        glDeleteBuffers(1, &vbo_);
        vbo_ = 0;
    }
    if (vao_ != 0) {
        glDeleteVertexArrays(1, &vao_);
        vao_ = 0;
    }
    if (program_ != 0) {
        glDeleteProgram(program_);
        program_ = 0;
    }

    vertices_.clear();
    vertexCount_ = 0;
}

void SceneRenderer::BuildTorusWire(const TokamakConfig& config) {
    constexpr int kMajorSegments = 96;
    constexpr int kMinorSegments = 32;

    const float r = config.minorRadius;
    const float R = config.majorRadius;
    const float twoPi = 2.0f * 3.14159265359f;

    for (int i = 0; i < kMajorSegments; ++i) {
        const float phi0 = twoPi * static_cast<float>(i) / static_cast<float>(kMajorSegments);
        const float phi1 = twoPi * static_cast<float>(i + 1) / static_cast<float>(kMajorSegments);

        for (int j = 0; j < kMinorSegments; ++j) {
            const float theta0 = twoPi * static_cast<float>(j) / static_cast<float>(kMinorSegments);
            const float theta1 = twoPi * static_cast<float>(j + 1) / static_cast<float>(kMinorSegments);

            const Vec3 p00 = TorusPoint(R, r, phi0, theta0);
            const Vec3 p10 = TorusPoint(R, r, phi1, theta0);
            const Vec3 p01 = TorusPoint(R, r, phi0, theta1);

            PushLine(vertices_, p00, p10, 0.28f, 0.95f, 0.72f);
            PushLine(vertices_, p00, p01, 0.14f, 0.65f, 0.58f);
        }
    }
}

void SceneRenderer::BuildInjectorMarker(const NBIConfig& nbi) {
    const Vec3 c = nbi.injectorPos;
    const float s = 0.09f;

    PushLine(vertices_, Vec3(c.x - s, c.y, c.z), Vec3(c.x + s, c.y, c.z), 1.0f, 0.82f, 0.25f);
    PushLine(vertices_, Vec3(c.x, c.y - s, c.z), Vec3(c.x, c.y + s, c.z), 1.0f, 0.82f, 0.25f);
    PushLine(vertices_, Vec3(c.x, c.y, c.z - s), Vec3(c.x, c.y, c.z + s), 1.0f, 0.82f, 0.25f);
}

void SceneRenderer::UpdateGeometry(const TokamakConfig& config, const NBIConfig& nbi) {
    vertices_.clear();
    vertices_.reserve(96 * 32 * 24);

    BuildTorusWire(config);
    BuildInjectorMarker(nbi);

    vertexCount_ = static_cast<int>(vertices_.size() / 6);

    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferData(
        GL_ARRAY_BUFFER,
        static_cast<GLsizeiptr>(vertices_.size() * sizeof(float)),
        vertices_.data(),
        GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void SceneRenderer::Draw(const Mat4& viewProjection) const {
    if (program_ == 0 || vao_ == 0 || vertexCount_ == 0) {
        return;
    }

    glUseProgram(program_);
    glUniformMatrix4fv(viewProjectionLoc_, 1, GL_FALSE, viewProjection.Data());

    glBindVertexArray(vao_);
    glDrawArrays(GL_LINES, 0, vertexCount_);
    glBindVertexArray(0);

    glUseProgram(0);
}
