#include "render/camera.h"

#include <algorithm>
#include <cmath>

namespace {

constexpr float kDegToRad = 3.14159265359f / 180.0f;

}  // namespace

OrbitCamera::OrbitCamera() = default;

void OrbitCamera::SetViewport(int width, int height) {
    viewportWidth_ = std::max(width, 1);
    viewportHeight_ = std::max(height, 1);
}

void OrbitCamera::SetTarget(const Vec3& target) {
    target_ = target;
}

void OrbitCamera::BeginDrag(int button, double x, double y) {
    if (button == 0) {
        leftDragging_ = true;
    } else if (button == 1) {
        rightDragging_ = true;
    }
    lastX_ = x;
    lastY_ = y;
}

void OrbitCamera::EndDrag(int button) {
    if (button == 0) {
        leftDragging_ = false;
    } else if (button == 1) {
        rightDragging_ = false;
    }
}

void OrbitCamera::OnCursorMove(double x, double y) {
    const float dx = static_cast<float>(x - lastX_);
    const float dy = static_cast<float>(y - lastY_);
    lastX_ = x;
    lastY_ = y;

    if (leftDragging_) {
        yawDegrees_ += dx * 0.2f;
        pitchDegrees_ += dy * 0.2f;
        pitchDegrees_ = std::clamp(pitchDegrees_, -89.0f, 89.0f);
    }

    if (rightDragging_) {
        const float yawRad = yawDegrees_ * kDegToRad;
        const Vec3 right(std::cos(yawRad), std::sin(yawRad), 0.0f);
        const Vec3 up(0.0f, 0.0f, 1.0f);
        const float panScale = distance_ * 0.0015f;

        target_ = target_ - (right * (dx * panScale));
        target_ = target_ + (up * (dy * panScale));
    }
}

void OrbitCamera::OnScroll(double yOffset) {
    const float scale = 1.0f - static_cast<float>(yOffset) * 0.08f;
    distance_ *= std::max(0.1f, scale);
    distance_ = std::clamp(distance_, 0.4f, 100.0f);
}

Vec3 OrbitCamera::Position() const {
    const float yaw = yawDegrees_ * kDegToRad;
    const float pitch = pitchDegrees_ * kDegToRad;

    const float cosPitch = std::cos(pitch);
    const float sinPitch = std::sin(pitch);
    const float cosYaw = std::cos(yaw);
    const float sinYaw = std::sin(yaw);

    const Vec3 offset(
        distance_ * cosPitch * cosYaw,
        distance_ * cosPitch * sinYaw,
        distance_ * sinPitch);

    return target_ + offset;
}

Mat4 OrbitCamera::ViewMatrix() const {
    return LookAt(Position(), target_, Vec3(0.0f, 0.0f, 1.0f));
}

Mat4 OrbitCamera::ProjectionMatrix() const {
    const float aspect = static_cast<float>(viewportWidth_) / static_cast<float>(viewportHeight_);
    return Perspective(fovYRadians_, aspect, 0.01f, 200.0f);
}

Mat4 OrbitCamera::ViewProjectionMatrix() const {
    return Multiply(ProjectionMatrix(), ViewMatrix());
}
