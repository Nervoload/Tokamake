#pragma once

#include <array>
#include <cmath>

#include "math_types.h"

struct Mat4 {
    std::array<float, 16> m{};

    static Mat4 Identity() {
        Mat4 out;
        out.m = {
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1,
        };
        return out;
    }

    const float* Data() const { return m.data(); }
};

inline Mat4 Multiply(const Mat4& a, const Mat4& b) {
    Mat4 out{};
    for (int col = 0; col < 4; ++col) {
        for (int row = 0; row < 4; ++row) {
            float sum = 0.0f;
            for (int k = 0; k < 4; ++k) {
                sum += a.m[k * 4 + row] * b.m[col * 4 + k];
            }
            out.m[col * 4 + row] = sum;
        }
    }
    return out;
}

inline Mat4 Perspective(float fovYRadians, float aspect, float zNear, float zFar) {
    Mat4 out{};
    const float f = 1.0f / std::tan(fovYRadians * 0.5f);
    out.m = {
        f / aspect, 0, 0, 0,
        0, f, 0, 0,
        0, 0, (zFar + zNear) / (zNear - zFar), -1,
        0, 0, (2.0f * zFar * zNear) / (zNear - zFar), 0,
    };
    return out;
}

inline Mat4 LookAt(const Vec3& eye, const Vec3& center, const Vec3& up) {
    const Vec3 f = (center - eye).normalize();
    const Vec3 s = Vec3::cross(f, up).normalize();
    const Vec3 u = Vec3::cross(s, f);

    Mat4 out = Mat4::Identity();
    out.m[0] = s.x;
    out.m[1] = u.x;
    out.m[2] = -f.x;

    out.m[4] = s.y;
    out.m[5] = u.y;
    out.m[6] = -f.y;

    out.m[8] = s.z;
    out.m[9] = u.z;
    out.m[10] = -f.z;

    out.m[12] = -Vec3::dot(s, eye);
    out.m[13] = -Vec3::dot(u, eye);
    out.m[14] = Vec3::dot(f, eye);
    return out;
}
