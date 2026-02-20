#include <cmath>
#include <limits>

#include "gtest/gtest.h"
#include "tokamak/magnetic_field.hpp"

namespace {

tokamak::TokamakConfig BaseConfig() {
    tokamak::TokamakConfig config;
    config.majorRadius_m = 2.0f;
    config.minorRadius_m = 0.5f;
    config.toroidalCurrent_A = 15.0e6f;
    config.toroidalCoilTurns = 18;
    config.plasmaCurrent_A = 2.0e6f;
    return config;
}

}  // namespace

TEST(Milestone2MagneticFieldTest, AxisBoundedNearMinorAxis) {
    const tokamak::TokamakConfig config = BaseConfig();
    tokamak::PlasmaCurrentProfileConfig profile;
    profile.kind = tokamak::PlasmaCurrentProfileKind::Uniform;
    profile.axisEpsilon_m = 1.0e-6f;

    const tokamak::Vec3 pos(config.majorRadius_m + 1.0e-8f, 0.0f, 1.0e-8f);
    const tokamak::MagneticFieldSample sample = tokamak::EvaluateMagneticFieldSample(config, profile, pos);

    EXPECT_TRUE(std::isfinite(sample.totalField_T.x));
    EXPECT_TRUE(std::isfinite(sample.totalField_T.y));
    EXPECT_TRUE(std::isfinite(sample.totalField_T.z));
    EXPECT_TRUE(std::isfinite(sample.totalMagnitude_T));
    EXPECT_TRUE(std::isfinite(sample.poloidalField_T.Magnitude()));
    EXPECT_LT(sample.poloidalField_T.Magnitude(), 1.0e3f);
}

TEST(Milestone2MagneticFieldTest, UniformProfileMatchesQuadraticEnclosedCurrent) {
    tokamak::PlasmaCurrentProfileConfig profile;
    profile.kind = tokamak::PlasmaCurrentProfileKind::Uniform;

    EXPECT_NEAR(tokamak::ComputeEnclosedCurrentFraction(0.0f, profile), 0.0f, 1.0e-7f);
    EXPECT_NEAR(tokamak::ComputeEnclosedCurrentFraction(0.5f, profile), 0.25f, 1.0e-6f);
    EXPECT_NEAR(tokamak::ComputeEnclosedCurrentFraction(1.0f, profile), 1.0f, 1.0e-6f);
    EXPECT_NEAR(tokamak::ComputeEnclosedCurrentFraction(1.5f, profile), 1.0f, 1.0e-6f);
}

TEST(Milestone2MagneticFieldTest, ParabolicProfileMatchesAnalyticEnclosedCurrent) {
    tokamak::PlasmaCurrentProfileConfig profile;
    profile.kind = tokamak::PlasmaCurrentProfileKind::Parabolic;

    EXPECT_NEAR(tokamak::ComputeEnclosedCurrentFraction(0.25f, profile), 0.12109375f, 1.0e-6f);
    EXPECT_NEAR(tokamak::ComputeEnclosedCurrentFraction(0.50f, profile), 0.4375f, 1.0e-6f);
    EXPECT_NEAR(tokamak::ComputeEnclosedCurrentFraction(0.75f, profile), 0.80859375f, 1.0e-6f);
    EXPECT_NEAR(tokamak::ComputeEnclosedCurrentFraction(1.00f, profile), 1.0f, 1.0e-6f);
}

TEST(Milestone2MagneticFieldTest, CustomProfileSanitizedMonotonicAndInterpolated) {
    tokamak::PlasmaCurrentProfileConfig profile;
    profile.kind = tokamak::PlasmaCurrentProfileKind::CustomTable;
    profile.customTable = {
        {1.0f, 0.6f},
        {0.5f, 0.7f},
        {0.2f, 0.3f},
        {0.8f, 0.2f},
    };

    const float at05 = tokamak::ComputeEnclosedCurrentFraction(0.5f, profile);
    const float at08 = tokamak::ComputeEnclosedCurrentFraction(0.8f, profile);
    const float at09 = tokamak::ComputeEnclosedCurrentFraction(0.9f, profile);

    EXPECT_NEAR(at05, 0.7f, 1.0e-6f);
    EXPECT_NEAR(at08, 0.7f, 1.0e-6f);
    EXPECT_NEAR(at09, 0.85f, 1.0e-6f);
    EXPECT_LE(at05, at08);
    EXPECT_LE(at08, at09);
}

TEST(Milestone2MagneticFieldTest, ContinuityNearRhoOneBoundary) {
    const tokamak::TokamakConfig config = BaseConfig();
    tokamak::PlasmaCurrentProfileConfig profile;
    profile.kind = tokamak::PlasmaCurrentProfileKind::Uniform;

    const float rInside = config.minorRadius_m * (1.0f - 1.0e-5f);
    const float rOutside = config.minorRadius_m * (1.0f + 1.0e-5f);

    const tokamak::Vec3 insidePos(config.majorRadius_m + rInside, 0.0f, 0.0f);
    const tokamak::Vec3 outsidePos(config.majorRadius_m + rOutside, 0.0f, 0.0f);

    const tokamak::MagneticFieldSample inside = tokamak::EvaluateMagneticFieldSample(config, profile, insidePos);
    const tokamak::MagneticFieldSample outside = tokamak::EvaluateMagneticFieldSample(config, profile, outsidePos);

    const float insidePol = inside.poloidalField_T.Magnitude();
    const float outsidePol = outside.poloidalField_T.Magnitude();
    EXPECT_GT(insidePol, 0.0f);
    EXPECT_GT(outsidePol, 0.0f);

    const float tolerance = outsidePol * 5.0e-3f;
    EXPECT_NEAR(insidePol, outsidePol, tolerance);
}

TEST(Milestone2MagneticFieldTest, CustomAxisBlendKeepsNearAxisCurrentSmall) {
    const tokamak::TokamakConfig config = BaseConfig();
    tokamak::PlasmaCurrentProfileConfig profile;
    profile.kind = tokamak::PlasmaCurrentProfileKind::CustomTable;
    profile.customAxisBlendRadius_m = 0.1f;
    profile.customTable = {
        {0.05f, 1.0f},
        {1.0f, 1.0f},
    };

    const tokamak::Vec3 nearAxisPos(config.majorRadius_m + 5.0e-4f, 0.0f, 0.0f);
    const tokamak::MagneticFieldSample sample =
        tokamak::EvaluateMagneticFieldSample(config, profile, nearAxisPos);

    EXPECT_LT(sample.normalizedMinorRadius, 0.01f);
    EXPECT_LT(sample.enclosedCurrentFraction, 1.0e-3f);
}

TEST(Milestone2MagneticFieldTest, RecommendDtMonotonicWithFieldStrength) {
    const double dtWeak = tokamak::RecommendDtFromMaxField(1.0, 0.05);
    const double dtStrong = tokamak::RecommendDtFromMaxField(10.0, 0.05);

    EXPECT_TRUE(std::isfinite(dtWeak));
    EXPECT_TRUE(std::isfinite(dtStrong));
    EXPECT_GT(dtWeak, dtStrong);
    EXPECT_EQ(tokamak::RecommendDtFromMaxField(0.0, 0.05), std::numeric_limits<double>::infinity());
    EXPECT_NEAR(tokamak::RecommendDtFromMaxField(1.0, -1.0), 0.0, 1.0e-12);
}
