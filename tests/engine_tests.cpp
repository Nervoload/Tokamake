#include <cstdint>
#include <cmath>
#include <vector>

#include "gtest/gtest.h"
#include "tokamak/engine.hpp"

TEST(EngineTest, NbiInsertionAtomicAtCapacity) {
    tokamak::RunConfig config;
    config.scenario = tokamak::Scenario::NbiIgnition;
    config.seed = 7;
    config.particleCap = 1001;

    tokamak::TokamakEngine engine(config);
    EXPECT_EQ(engine.Particles().Size(), static_cast<std::size_t>(1000));

    engine.Step(config.timeStep_s);

    EXPECT_EQ(engine.Particles().Size(), static_cast<std::size_t>(1000));
    EXPECT_GT(engine.Counters().rejectedInjectionPairs, static_cast<uint64_t>(0));
    EXPECT_GT(engine.Counters().particleCapHitEvents, static_cast<uint64_t>(0));
}

TEST(EngineTest, CollisionShuffleDeterministicWithSeed) {
    tokamak::RunConfig config;
    config.scenario = tokamak::Scenario::NbiIgnition;
    config.seed = 42;

    tokamak::TokamakEngine engineA(config);
    tokamak::TokamakEngine engineB(config);

    std::vector<uint64_t> fusionSeriesA;
    std::vector<uint64_t> fusionSeriesB;

    for (int step = 1; step <= 300; ++step) {
        engineA.Step(config.timeStep_s);
        engineB.Step(config.timeStep_s);

        if (step % 25 == 0) {
            fusionSeriesA.push_back(engineA.Snapshot(step).fusionEvents);
            fusionSeriesB.push_back(engineB.Snapshot(step).fusionEvents);
            EXPECT_EQ(engineA.Snapshot(step).species.AliveCount(), engineB.Snapshot(step).species.AliveCount());
        }
    }

    EXPECT_EQ(fusionSeriesA.size(), fusionSeriesB.size());
    for (std::size_t i = 0; i < fusionSeriesA.size(); ++i) {
        EXPECT_EQ(fusionSeriesA[i], fusionSeriesB[i]);
    }

    EXPECT_EQ(engineA.Counters().fusionAttempts, engineB.Counters().fusionAttempts);
    EXPECT_EQ(engineA.Counters().fusionAccepted, engineB.Counters().fusionAccepted);
    EXPECT_EQ(engineA.Counters().maxReactionsInCell, engineB.Counters().maxReactionsInCell);
}

TEST(EngineTest, EngineNoNaNInfLongRun) {
    tokamak::RunConfig config;
    config.scenario = tokamak::Scenario::ColdVacuum;
    config.seed = 99;
    config.totalSteps = 10000;

    tokamak::TokamakEngine engine(config);
    for (int step = 0; step < config.totalSteps; ++step) {
        engine.Step(config.timeStep_s);
        if (step % 1000 == 0) {
            EXPECT_TRUE(engine.HasFiniteState());
        }
    }

    EXPECT_TRUE(engine.HasFiniteState());
}

TEST(EngineTest, SortGridAccumulatesOutOfDomainClampCounters) {
    tokamak::RunConfig config;
    config.scenario = tokamak::Scenario::ColdVacuum;
    config.seed = 1234;

    tokamak::TokamakEngine engine(config);
    EXPECT_GT(engine.Particles().Size(), static_cast<std::size_t>(0));

    engine.MutableParticles().MutablePositions()[0] = tokamak::Vec3(0.0f, 0.0f, 1000.0f);
    engine.MutableParticles().MutableVelocities()[0] = tokamak::Vec3(0.0f, 0.0f, 0.0f);

    engine.Step(config.timeStep_s);
    const tokamak::TelemetrySnapshot snapshot = engine.Snapshot(1);

    EXPECT_GT(snapshot.stepCounters.outOfDomainCellClampEvents, static_cast<uint64_t>(0));
    EXPECT_GT(snapshot.counters.outOfDomainCellClampEvents, static_cast<uint64_t>(0));
}

TEST(EngineTest, WallBridgeCountersAccumulateWhileReflectBehaviorRemainsActive) {
    tokamak::RunConfig config;
    config.scenario = tokamak::Scenario::ColdVacuum;
    config.seed = 1337;
    config.wallBoundaryMode = tokamak::WallBoundaryMode::Recycle;
    config.recycleFraction = 0.5;

    tokamak::TokamakEngine engine(config);
    //ASSERT_GT(engine.Particles().Size(), static_cast<std::size_t>(0));

    // Force a near-wall outward trajectory to guarantee at least one reflective hit.
    engine.MutableParticles().MutablePositions()[0] = tokamak::Vec3(2.49f, 0.0f, 0.49f);
    engine.MutableParticles().MutableVelocities()[0] = tokamak::Vec3(3.0e5f, 0.0f, 3.0e5f);

    engine.Step(config.timeStep_s);
    const tokamak::TelemetrySnapshot snapshot = engine.Snapshot(1);

    EXPECT_GT(snapshot.stepCounters.wallHitCount, static_cast<uint64_t>(0));
    EXPECT_GT(snapshot.stepCounters.wallImpactEnergy_J, 0.0);
    EXPECT_EQ(snapshot.stepCounters.wallLossWeight, 0.0);
    EXPECT_EQ(snapshot.counters.wallLossWeight, 0.0);
}

TEST(EngineTest, MagneticDiagnosticsFiniteAndPositiveWhenFieldPresent) {
    tokamak::RunConfig config;
    config.scenario = tokamak::Scenario::ColdVacuum;
    config.seed = 1001;
    config.totalSteps = 200;
    config.magneticFieldRadialBinCount = 16;
    config.magneticFieldDtSafetyFraction = 0.05;

    tokamak::TokamakEngine engine(config);
    for (int step = 1; step <= config.totalSteps; ++step) {
        engine.Step(config.timeStep_s);
        if (step % 25 == 0) {
            const tokamak::TelemetrySnapshot snapshot = engine.Snapshot(step);
            EXPECT_TRUE(std::isfinite(snapshot.magneticField.maxField_T));
            EXPECT_GT(snapshot.magneticField.maxField_T, 0.0);
            EXPECT_TRUE(std::isfinite(snapshot.magneticField.recommendedDt_s));
            EXPECT_GT(snapshot.magneticField.recommendedDt_s, 0.0);
            EXPECT_EQ(snapshot.magneticField.radialMeanField_T.size(), config.magneticFieldRadialBinCount);
            EXPECT_EQ(snapshot.magneticField.radialSampleCounts.size(), config.magneticFieldRadialBinCount);
        }
    }
}

TEST(EngineTest, ProfileSwitchChangesMagneticDiagnosticsDeterministically) {
    tokamak::RunConfig uniformConfig;
    uniformConfig.scenario = tokamak::Scenario::NbiIgnition;
    uniformConfig.seed = 222;
    uniformConfig.totalSteps = 80;
    uniformConfig.magneticFieldRadialBinCount = 32;
    uniformConfig.plasmaCurrentProfileKind = tokamak::PlasmaCurrentProfileKind::Uniform;

    tokamak::RunConfig parabolicConfig = uniformConfig;
    parabolicConfig.plasmaCurrentProfileKind = tokamak::PlasmaCurrentProfileKind::Parabolic;

    tokamak::TokamakEngine uniformEngine(uniformConfig);
    tokamak::TokamakEngine parabolicEngine(parabolicConfig);

    for (int step = 1; step <= uniformConfig.totalSteps; ++step) {
        uniformEngine.Step(uniformConfig.timeStep_s);
        parabolicEngine.Step(parabolicConfig.timeStep_s);
    }

    const tokamak::TelemetrySnapshot uniformSnapshot = uniformEngine.Snapshot(uniformConfig.totalSteps);
    const tokamak::TelemetrySnapshot parabolicSnapshot = parabolicEngine.Snapshot(parabolicConfig.totalSteps);

    EXPECT_TRUE(std::isfinite(uniformSnapshot.magneticField.maxField_T));
    EXPECT_TRUE(std::isfinite(parabolicSnapshot.magneticField.maxField_T));
    EXPECT_NE(uniformSnapshot.magneticField.maxField_T, parabolicSnapshot.magneticField.maxField_T);
    EXPECT_NE(uniformSnapshot.magneticField.recommendedDt_s, parabolicSnapshot.magneticField.recommendedDt_s);
}

TEST(EngineTest, ElectrostaticModeProducesMeasuredResidualAndFiniteDiagnostics) {
    tokamak::RunConfig config;
    config.scenario = tokamak::Scenario::ColdVacuum;
    config.seed = 2048;
    config.electricFieldMode = tokamak::ElectricFieldMode::Electrostatic;
    config.electrostaticGridBinCount = 8;
    config.electrostaticSolverMaxIterations = 120;
    config.electrostaticSolverTolerance = 1.0e-5;
    config.electrostaticSorOmega = 1.5;
    config.chargeAssignmentScheme = tokamak::ChargeAssignmentScheme::CIC;

    tokamak::TokamakEngine engine(config);
    for (int step = 1; step <= 6; ++step) {
        engine.Step(config.timeStep_s);
    }

    const tokamak::TelemetrySnapshot snapshot = engine.Snapshot(6);
    EXPECT_EQ(engine.CurrentElectricFieldMode(), tokamak::ElectricFieldMode::Electrostatic);
    EXPECT_EQ(snapshot.solverResidual.status, tokamak::ResidualStatus::Measured);
    EXPECT_TRUE(snapshot.solverResidual.residualAvailable);
    EXPECT_TRUE(std::isfinite(snapshot.solverResidual.residualL2));
    EXPECT_GT(snapshot.electrostaticField.solveIterations, static_cast<uint32_t>(0));
    EXPECT_TRUE(std::isfinite(snapshot.electrostaticField.maxElectricField_VPerM));
    EXPECT_TRUE(std::isfinite(snapshot.electrostaticField.meanElectricField_VPerM));
    EXPECT_GE(snapshot.electrostaticField.maxElectricField_VPerM, 0.0);
    EXPECT_GE(snapshot.electrostaticField.meanElectricField_VPerM, 0.0);
}

TEST(EngineTest, PlaceholderModeRemainsDeterministicWhenExplicitlySelected) {
    tokamak::RunConfig defaultConfig;
    defaultConfig.scenario = tokamak::Scenario::NbiIgnition;
    defaultConfig.seed = 777;

    tokamak::RunConfig explicitPlaceholder = defaultConfig;
    explicitPlaceholder.electricFieldMode = tokamak::ElectricFieldMode::Placeholder;

    tokamak::TokamakEngine engineDefault(defaultConfig);
    tokamak::TokamakEngine engineExplicit(explicitPlaceholder);

    for (int step = 1; step <= 40; ++step) {
        engineDefault.Step(defaultConfig.timeStep_s);
        engineExplicit.Step(explicitPlaceholder.timeStep_s);
    }

    const tokamak::TelemetrySnapshot defaultSnapshot = engineDefault.Snapshot(40);
    const tokamak::TelemetrySnapshot explicitSnapshot = engineExplicit.Snapshot(40);

    EXPECT_EQ(defaultSnapshot.species.AliveCount(), explicitSnapshot.species.AliveCount());
    EXPECT_EQ(defaultSnapshot.fusionEvents, explicitSnapshot.fusionEvents);
    EXPECT_EQ(defaultSnapshot.counters.fusionAccepted, explicitSnapshot.counters.fusionAccepted);
    EXPECT_EQ(defaultSnapshot.solverResidual.status, tokamak::ResidualStatus::Placeholder);
    EXPECT_EQ(explicitSnapshot.solverResidual.status, tokamak::ResidualStatus::Placeholder);
}
