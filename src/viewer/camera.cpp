#include "tokamak/viewer/camera.hpp"

#include <algorithm>
#include <cmath>

namespace tokamak::viewer {
namespace {

Mat4 Identity() {
    Mat4 out;
    out.elements = {
        1.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f,
    };
    return out;
}

}  // namespace

OrbitCamera::OrbitCamera() = default;

void OrbitCamera::SetViewport(int width, int height) {
    viewportWidth_ = std::max(1, width);
    viewportHeight_ = std::max(1, height);
}

void OrbitCamera::BeginRotate(double x, double y) {
    rotating_ = true;
    lastCursorX_ = x;
    lastCursorY_ = y;
}

void OrbitCamera::EndRotate() {
    rotating_ = false;
}

void OrbitCamera::OnCursorMove(double x, double y) {
    if (!rotating_) {
        return;
    }

    const double deltaX = x - lastCursorX_;
    const double deltaY = y - lastCursorY_;
    lastCursorX_ = x;
    lastCursorY_ = y;

    yawRadians_ += static_cast<float>(deltaX * 0.006);
    pitchRadians_ += static_cast<float>(-deltaY * 0.006);

    const float minPitch = -1.4f;
    const float maxPitch = 1.4f;
    pitchRadians_ = std::max(minPitch, std::min(maxPitch, pitchRadians_));
}

void OrbitCamera::OnScroll(double deltaY) {
    const float zoomFactor = std::pow(0.92f, static_cast<float>(deltaY));
    distance_ *= zoomFactor;
    distance_ = std::max(0.8f, std::min(20.0f, distance_));
}

Mat4 OrbitCamera::ViewProjectionMatrix() const {
    const float cYaw = std::cos(yawRadians_);
    const float sYaw = std::sin(yawRadians_);
    const float cPitch = std::cos(pitchRadians_);
    const float sPitch = std::sin(pitchRadians_);

    const Vec3 eye(
        target_.x + distance_ * cPitch * cYaw,
        target_.y + distance_ * cPitch * sYaw,
        target_.z + distance_ * sPitch);

    const Vec3 up(0.0f, 0.0f, 1.0f);
    const float aspect = static_cast<float>(viewportWidth_) / static_cast<float>(viewportHeight_);

    const Mat4 projection = Perspective(45.0f * (3.14159265359f / 180.0f), aspect, 0.01f, 100.0f);
    const Mat4 view = LookAt(eye, target_, up);
    return Multiply(projection, view);
}

Mat4 Multiply(const Mat4& a, const Mat4& b) {
    Mat4 out;
    for (int col = 0; col < 4; ++col) {
        for (int row = 0; row < 4; ++row) {
            float sum = 0.0f;
            for (int k = 0; k < 4; ++k) {
                sum += a.elements[k * 4 + row] * b.elements[col * 4 + k];
            }
            out.elements[col * 4 + row] = sum;
        }
    }
    return out;
}

Mat4 Perspective(float fovYRadians, float aspect, float nearPlane, float farPlane) {
    Mat4 out{};
    const float tanHalf = std::tan(fovYRadians * 0.5f);

    out.elements[0] = 1.0f / (aspect * tanHalf);
    out.elements[5] = 1.0f / tanHalf;
    out.elements[10] = -(farPlane + nearPlane) / (farPlane - nearPlane);
    out.elements[11] = -1.0f;
    out.elements[14] = -(2.0f * farPlane * nearPlane) / (farPlane - nearPlane);
    return out;
}

Mat4 LookAt(const Vec3& eye, const Vec3& center, const Vec3& up) {
    const Vec3 f = (center - eye).Normalized();
    const Vec3 s = Vec3::Cross(f, up).Normalized();
    const Vec3 u = Vec3::Cross(s, f);

    Mat4 out = Identity();
    out.elements[0] = s.x;
    out.elements[4] = s.y;
    out.elements[8] = s.z;

    out.elements[1] = u.x;
    out.elements[5] = u.y;
    out.elements[9] = u.z;

    out.elements[2] = -f.x;
    out.elements[6] = -f.y;
    out.elements[10] = -f.z;

    out.elements[12] = -Vec3::Dot(s, eye);
    out.elements[13] = -Vec3::Dot(u, eye);
    out.elements[14] = Vec3::Dot(f, eye);
    return out;
}

}  // namespace tokamak::viewer
