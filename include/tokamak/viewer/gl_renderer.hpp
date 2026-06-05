#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "tokamak/viewer/camera.hpp"
#include "tokamak/viewer/replay_analytics.hpp"
#include "tokamak/viewer/replay_snapshot.hpp"

namespace tokamak::viewer {

enum class ParticleViewMode : uint8_t {
    Species = 0,
    Energy = 1,
    PitchAngle = 2,
    MagneticMagnitude = 3,
    ElectricMagnitude = 4,
    LorentzAcceleration = 5,
    FusionRate = 6,
    WallLossRisk = 7,
};

struct ParticleViewRange {
    double minValue = 0.0;
    double maxValue = 0.0;
    bool valid = false;
};

struct ParticleSliceSettings {
    bool enableToroidalSlice = false;
    bool enablePoloidalSlice = false;
    bool dimOutsideSlice = true;
    float toroidalCenter_deg = 0.0f;
    float toroidalHalfWidth_deg = 18.0f;
    float poloidalCenter_deg = 0.0f;
    float poloidalHalfWidth_deg = 28.0f;
};

struct ParticleViewContext {
    ParticleViewMode mode = ParticleViewMode::Species;
    ParticleViewRange range;
    const std::vector<ReplayRadialProfileBin>* radialProfile = nullptr;
    const std::vector<ReplayRadialProfileBin>* compareRadialProfile = nullptr;
    const ReplayWallInteractionPoint* wallInteraction = nullptr;
    const ReplayWallInteractionPoint* compareWallInteraction = nullptr;
    float majorRadius_m = 2.0f;
    float minorRadius_m = 0.5f;
    bool diffAgainstCompare = false;
    ParticleSliceSettings slice;
};

struct GraphicsStyleSettings {
    bool showVesselShell = true;
    bool showParticleTrails = true;
    bool showDensitySplat = true;
    bool xrayOverlays = true;
    float shellOpacity = 0.06f;
    float fogStrength = 0.18f;
    float overlayFogStrength = 0.02f;
    float pointFogStrength = 0.10f;
    float sceneLineOpacity = 0.44f;
    float overlayLineOpacity = 0.96f;
    float trailOpacity = 0.74f;
    float densitySplatOpacity = 0.08f;
    float densitySplatScale = 2.4f;
    float particleColorBoost = 1.52f;
    float sceneLineWidth = 1.2f;
    float overlayLineWidth = 2.0f;
    float trailLineWidth = 1.9f;
};

const char* ParticleViewModeName(ParticleViewMode mode);
ParticleViewRange ComputeParticleViewRange(const ReplayFrame& frame, const ParticleViewContext& context);
void ComputeParticleViewColor(
    const ReplayParticle& particle,
    const ParticleViewContext& context,
    float* r,
    float* g,
    float* b);

class GlRenderer {
public:
    bool Initialize(std::size_t maxParticles, std::string* errorOut);
    void Shutdown();

    void SetTorusGeometry(float majorRadius_m, float minorRadius_m);
    void UploadFrame(
        const ReplayFrame& frame,
        const ReplayFrame* nextFrame,
        const ParticleViewContext& viewContext);
    void UploadOverlayLines(const std::vector<float>& overlayVertices);

    void Draw(
        const Mat4& viewProjection,
        float pointSizePixels,
        float interpolationAlpha,
        float ignitionIntensity,
        const GraphicsStyleSettings& graphicsStyle) const;

private:
    bool InitializeParticlePipeline(std::string* errorOut);
    bool InitializeLinePipeline(std::string* errorOut);
    bool InitializeShellPipeline(std::string* errorOut);
    void EnsureParticleCapacity(std::size_t requiredParticles);

    void RebuildSceneGeometry(float majorRadius_m, float minorRadius_m);

    unsigned int pointProgram_ = 0;
    unsigned int pointVao_ = 0;
    unsigned int pointVbo_ = 0;
    int pointViewProjectionLocation_ = -1;
    int pointSizeLocation_ = -1;
    int pointInterpolationAlphaLocation_ = -1;
    int pointFogColorLocation_ = -1;
    int pointAlphaScaleLocation_ = -1;
    int pointColorBoostLocation_ = -1;
    int pointFogMixLocation_ = -1;

    unsigned int lineProgram_ = 0;
    unsigned int lineVao_ = 0;
    unsigned int lineVbo_ = 0;
    int lineViewProjectionLocation_ = -1;
    int lineFogColorLocation_ = -1;
    int lineAlphaScaleLocation_ = -1;
    int lineFogMixLocation_ = -1;

    unsigned int overlayLineVao_ = 0;
    unsigned int overlayLineVbo_ = 0;
    int overlayLineVertexCount_ = 0;

    unsigned int shellProgram_ = 0;
    unsigned int shellVao_ = 0;
    unsigned int shellVbo_ = 0;
    int shellVertexCount_ = 0;
    int shellViewProjectionLocation_ = -1;
    int shellFogColorLocation_ = -1;
    int shellAlphaScaleLocation_ = -1;
    int shellFogMixLocation_ = -1;

    unsigned int trailVao_ = 0;
    unsigned int trailVbo_ = 0;
    int trailVertexCount_ = 0;

    std::size_t maxParticles_ = 0;
    std::size_t uploadedParticles_ = 0;
    int lineVertexCount_ = 0;
    float torusMajorRadius_m_ = 2.0f;
    float torusMinorRadius_m_ = 0.5f;
    std::vector<float> particleVertexBuffer_;
    std::vector<float> trailVertexBuffer_;
    std::vector<float> overlayLineVertexBuffer_;
};

}  // namespace tokamak::viewer
