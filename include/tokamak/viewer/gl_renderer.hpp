#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

#include "tokamak/viewer/camera.hpp"
#include "tokamak/viewer/replay_snapshot.hpp"

namespace tokamak::viewer {

class GlRenderer {
public:
    bool Initialize(std::size_t maxParticles, std::string* errorOut);
    void Shutdown();

    void SetTorusGeometry(float majorRadius_m, float minorRadius_m);
    void UploadFrame(const ReplayFrame& frame);

    void Draw(const Mat4& viewProjection, float pointSizePixels) const;

private:
    bool InitializeParticlePipeline(std::string* errorOut);
    bool InitializeLinePipeline(std::string* errorOut);

    void RebuildTorusLines(float majorRadius_m, float minorRadius_m);

    unsigned int pointProgram_ = 0;
    unsigned int pointVao_ = 0;
    unsigned int pointVbo_ = 0;
    int pointViewProjectionLocation_ = -1;
    int pointSizeLocation_ = -1;

    unsigned int lineProgram_ = 0;
    unsigned int lineVao_ = 0;
    unsigned int lineVbo_ = 0;
    int lineViewProjectionLocation_ = -1;

    std::size_t maxParticles_ = 0;
    std::size_t uploadedParticles_ = 0;
    int lineVertexCount_ = 0;
};

}  // namespace tokamak::viewer
