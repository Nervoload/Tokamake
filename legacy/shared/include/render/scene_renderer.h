#pragma once

#include <vector>

#include "render/render_math.h"
#include "tokamak_types.h"

class SceneRenderer {
public:
    bool Initialize();
    void Shutdown();

    void UpdateGeometry(const TokamakConfig& config, const NBIConfig& nbi);
    void Draw(const Mat4& viewProjection) const;

private:
    void BuildTorusWire(const TokamakConfig& config);
    void BuildInjectorMarker(const NBIConfig& nbi);

    unsigned int program_ = 0;
    unsigned int vao_ = 0;
    unsigned int vbo_ = 0;

    std::vector<float> vertices_;
    int vertexCount_ = 0;
    int viewProjectionLoc_ = -1;
};
