#include "tokamak/particle_push.hpp"

#include <cmath>

namespace tokamak {

Vec3 BorisVelocityStep(
    const Vec3& velocity,
    const Vec3& electricField_VPerM,
    const Vec3& magneticField_T,
    float chargeToMass_CPerKg,
    float dt_s) {
    const float halfDt = dt_s * 0.5f;
    const Vec3 vMinus = velocity + (electricField_VPerM * (chargeToMass_CPerKg * halfDt));
    const Vec3 t = magneticField_T * (chargeToMass_CPerKg * halfDt);
    const Vec3 vPrime = vMinus + Vec3::Cross(vMinus, t);
    const Vec3 s = t * (2.0f / (1.0f + Vec3::Dot(t, t)));
    const Vec3 vPlus = vMinus + Vec3::Cross(vPrime, s);
    return vPlus + (electricField_VPerM * (chargeToMass_CPerKg * halfDt));
}

bool ReflectAtTokamakWall(
    Vec3* position,
    Vec3* velocity,
    float majorRadius_m,
    float minorRadius_m,
    float wallPlacementScale) {
    if (position == nullptr || velocity == nullptr) {
        return false;
    }

    const float radialPos = std::sqrt((position->x * position->x) + (position->y * position->y));
    const float radialTube = std::sqrt(std::pow(radialPos - majorRadius_m, 2.0f) + (position->z * position->z));
    if (radialTube < minorRadius_m || radialPos <= 0.001f) {
        return false;
    }

    const Vec3 coreCenter(majorRadius_m * (position->x / radialPos), majorRadius_m * (position->y / radialPos), 0.0f);
    const Vec3 normal = (*position - coreCenter).Normalized();

    *velocity = *velocity - (normal * (2.0f * Vec3::Dot(*velocity, normal)));
    *position = coreCenter + (normal * (minorRadius_m * wallPlacementScale));
    return true;
}

}  // namespace tokamak
