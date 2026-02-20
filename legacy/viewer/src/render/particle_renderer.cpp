#include "render/particle_renderer.h"

#include <algorithm>
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
        std::string log(static_cast<size_t>(std::max(logLength, 1)), '\0');
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
        std::string log(static_cast<size_t>(std::max(logLength, 1)), '\0');
        glGetProgramInfoLog(program, logLength, nullptr, log.data());
        std::cerr << "Program link failure: " << log << std::endl;
        glDeleteProgram(program);
        return 0;
    }

    return program;
}

const char* kParticleVertexShader = R"GLSL(
#version 330 core
layout(location = 0) in vec3 aPosition;
layout(location = 1) in float aSpecies;

uniform mat4 uViewProjection;
uniform float uPointSize;

out vec3 vColor;

vec3 SpeciesColor(float species) {
    if (species < 0.5) {
        return vec3(0.22, 0.62, 1.0); // Deuterium
    }
    if (species < 1.5) {
        return vec3(1.0, 0.53, 0.22); // Tritium
    }
    if (species < 2.5) {
        return vec3(0.95, 0.95, 0.35); // Helium ash
    }
    return vec3(0.4, 0.4, 0.4);
}

void main() {
    gl_Position = uViewProjection * vec4(aPosition, 1.0);
    gl_PointSize = uPointSize;
    vColor = SpeciesColor(aSpecies);
}
)GLSL";

const char* kParticleFragmentShader = R"GLSL(
#version 330 core
in vec3 vColor;
out vec4 fragColor;

void main() {
    vec2 centered = gl_PointCoord * 2.0 - 1.0;
    float r2 = dot(centered, centered);
    if (r2 > 1.0) {
        discard;
    }
    fragColor = vec4(vColor, 0.92);
}
)GLSL";

}  // namespace

bool ParticleRenderer::Initialize(size_t maxParticles) {
    maxParticles_ = std::max<size_t>(1, maxParticles);

    program_ = CreateProgram(kParticleVertexShader, kParticleFragmentShader);
    if (program_ == 0) {
        return false;
    }

    viewProjectionLoc_ = glGetUniformLocation(program_, "uViewProjection");
    pointSizeLoc_ = glGetUniformLocation(program_, "uPointSize");

    glGenVertexArrays(1, &vao_);
    glGenBuffers(1, &vbo_);

    glBindVertexArray(vao_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);

    glBufferData(
        GL_ARRAY_BUFFER,
        static_cast<GLsizeiptr>(maxParticles_ * 4 * sizeof(float)),
        nullptr,
        GL_DYNAMIC_DRAW);

    constexpr GLsizei stride = static_cast<GLsizei>(4 * sizeof(float));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(0));

    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(3 * sizeof(float)));

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    interleaved_.reserve(maxParticles_ * 4);
    uploadedCount_ = 0;
    return true;
}

void ParticleRenderer::Shutdown() {
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
    interleaved_.clear();
    uploadedCount_ = 0;
}

void ParticleRenderer::Upload(const std::vector<ParticleSample>& particles) {
    const size_t count = std::min(particles.size(), maxParticles_);
    uploadedCount_ = count;

    interleaved_.resize(count * 4);
    for (size_t i = 0; i < count; ++i) {
        const ParticleSample& p = particles[i];
        interleaved_[i * 4 + 0] = p.x;
        interleaved_[i * 4 + 1] = p.y;
        interleaved_[i * 4 + 2] = p.z;
        interleaved_[i * 4 + 3] = static_cast<float>(p.species);
    }

    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferSubData(
        GL_ARRAY_BUFFER,
        0,
        static_cast<GLsizeiptr>(interleaved_.size() * sizeof(float)),
        interleaved_.data());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void ParticleRenderer::Draw(const Mat4& viewProjection, float pointSizePixels) const {
    if (program_ == 0 || vao_ == 0 || uploadedCount_ == 0) {
        return;
    }

    glUseProgram(program_);
    glUniformMatrix4fv(viewProjectionLoc_, 1, GL_FALSE, viewProjection.Data());
    glUniform1f(pointSizeLoc_, pointSizePixels);

    glBindVertexArray(vao_);
    glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(uploadedCount_));
    glBindVertexArray(0);

    glUseProgram(0);
}
