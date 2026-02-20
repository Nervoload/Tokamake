#pragma once

#include "render/render_math.h"

class OrbitCamera {
public:
    OrbitCamera();

    void SetViewport(int width, int height);
    void SetTarget(const Vec3& target);

    void BeginDrag(int button, double x, double y);
    void EndDrag(int button);
    void OnCursorMove(double x, double y);
    void OnScroll(double yOffset);

    Mat4 ViewMatrix() const;
    Mat4 ProjectionMatrix() const;
    Mat4 ViewProjectionMatrix() const;
    Vec3 Position() const;

private:
    Vec3 target_{0.0f, 0.0f, 0.0f};
    float yawDegrees_ = 35.0f;
    float pitchDegrees_ = 25.0f;
    float distance_ = 6.0f;

    float fovYRadians_ = 45.0f * 3.14159265359f / 180.0f;
    int viewportWidth_ = 1280;
    int viewportHeight_ = 720;

    bool leftDragging_ = false;
    bool rightDragging_ = false;
    double lastX_ = 0.0;
    double lastY_ = 0.0;
};
