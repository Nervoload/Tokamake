#pragma once

#include <cmath>

namespace tokamak {

struct Vec3 {
    float x;
    float y;
    float z;

    Vec3() : x(0.0f), y(0.0f), z(0.0f) {}
    Vec3(float xIn, float yIn, float zIn) : x(xIn), y(yIn), z(zIn) {}

    Vec3 operator+(const Vec3& rhs) const { return Vec3(x + rhs.x, y + rhs.y, z + rhs.z); }
    Vec3 operator-(const Vec3& rhs) const { return Vec3(x - rhs.x, y - rhs.y, z - rhs.z); }
    Vec3 operator*(float scale) const { return Vec3(x * scale, y * scale, z * scale); }
    Vec3 operator/(float scale) const { return Vec3(x / scale, y / scale, z / scale); }

    Vec3& operator+=(const Vec3& rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    float Magnitude() const { return std::sqrt((x * x) + (y * y) + (z * z)); }

    Vec3 Normalized() const {
        const float m = Magnitude();
        return (m > 0.0f) ? (*this / m) : Vec3();
    }

    bool IsFinite() const {
        return std::isfinite(x) && std::isfinite(y) && std::isfinite(z);
    }

    static float Dot(const Vec3& a, const Vec3& b) {
        return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
    }

    static Vec3 Cross(const Vec3& a, const Vec3& b) {
        return Vec3(
            (a.y * b.z) - (a.z * b.y),
            (a.z * b.x) - (a.x * b.z),
            (a.x * b.y) - (a.y * b.x));
    }
};

}  // namespace tokamak
