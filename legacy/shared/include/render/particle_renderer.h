#pragma once

#include <cstddef>
#include <vector>

#include "render/render_math.h"
#include "sim_snapshot.h"

class ParticleRenderer {
public:
    bool Initialize(size_t maxParticles);
    void Shutdown();

    void Upload(const std::vector<ParticleSample>& particles);
    void Draw(const Mat4& viewProjection, float pointSizePixels) const;

    size_t UploadedCount() const { return uploadedCount_; }

private:
    unsigned int program_ = 0;
    unsigned int vao_ = 0;
    unsigned int vbo_ = 0;

    size_t maxParticles_ = 0;
    size_t uploadedCount_ = 0;
    int viewProjectionLoc_ = -1;
    int pointSizeLoc_ = -1;
    std::vector<float> interleaved_;
};
