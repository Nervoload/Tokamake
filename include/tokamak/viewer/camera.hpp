#pragma once

#include <array>

#include "tokamak/types.hpp"

namespace tokamak::viewer {

struct Mat4 {
    std::array<float, 16> elements{};

    const float* Data() const { return elements.data(); }
};

class OrbitCamera {
public:
    OrbitCamera();

    void SetViewport(int width, int height);
    void BeginRotate(double x, double y);
    void EndRotate();
    void OnCursorMove(double x, double y);
    void OnScroll(double deltaY);

    Mat4 ViewProjectionMatrix() const;

private:
    int viewportWidth_ = 1600;
    int viewportHeight_ = 900;

    bool rotating_ = false;
    double lastCursorX_ = 0.0;
    double lastCursorY_ = 0.0;

    float yawRadians_ = 0.35f;
    float pitchRadians_ = 0.45f;
    float distance_ = 4.2f;
    tokamak::Vec3 target_ = tokamak::Vec3(0.0f, 0.0f, 0.0f);
};

Mat4 Multiply(const Mat4& a, const Mat4& b);
Mat4 Perspective(float fovYRadians, float aspect, float nearPlane, float farPlane);
Mat4 LookAt(const tokamak::Vec3& eye, const tokamak::Vec3& center, const tokamak::Vec3& up);

}  // namespace tokamak::viewer
