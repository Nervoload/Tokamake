#pragma once

#include "tokamak/types.hpp"

namespace tokamak {

Vec3 BorisVelocityStep(
    const Vec3& velocity,
    const Vec3& electricField_VPerM,
    const Vec3& magneticField_T,
    float chargeToMass_CPerKg,
    float dt_s);

bool ReflectAtTokamakWall(
    Vec3* position,
    Vec3* velocity,
    float majorRadius_m,
    float minorRadius_m,
    float wallPlacementScale = 0.99f);

}  // namespace tokamak
