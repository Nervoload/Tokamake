#include "tokamak/engine.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>

#include "tokamak/particle_push.hpp"

namespace tokamak {
namespace {

double ClampUnit(double value) {
    return std::max(0.0, std::min(value, 1.0));
}

double SmoothStepUnit(double value) {
    const double clamped = ClampUnit(value);
    return clamped * clamped * (3.0 - (2.0 * clamped));
}

struct TokamakWallHit {
    Vec3 projectedPosition;
    Vec3 normal;
    Vec3 coreCenter;
};

bool ComputeTokamakWallHit(
    const Vec3& position,
    float majorRadius_m,
    float minorRadius_m,
    TokamakWallHit* outHit,
    float wallPlacementScale = 0.99f) {
    const float radialPos = std::sqrt((position.x * position.x) + (position.y * position.y));
    if (radialPos <= 1.0e-6f) {
        return false;
    }

    const Vec3 coreCenter(
        majorRadius_m * (position.x / radialPos),
        majorRadius_m * (position.y / radialPos),
        0.0f);
    const Vec3 displacement = position - coreCenter;
    const float radialTube = displacement.Magnitude();
    if (radialTube < minorRadius_m) {
        return false;
    }

    TokamakWallHit hit;
    hit.coreCenter = coreCenter;
    hit.normal = displacement.Normalized();
    hit.projectedPosition = coreCenter + (hit.normal * (minorRadius_m * wallPlacementScale));

    if (outHit != nullptr) {
        *outHit = hit;
    }
    return true;
}

Vec3 ToroidalTangent(const Vec3& coreCenter) {
    const Vec3 tangent(-coreCenter.y, coreCenter.x, 0.0f);
    const Vec3 normalized = tangent.Normalized();
    if (normalized.Magnitude() <= 1.0e-6f) {
        return Vec3(0.0f, 1.0f, 0.0f);
    }
    return normalized;
}

}  // namespace

TokamakEngine::TokamakEngine(const RunConfig& runConfig)
    : TokamakEngine(
          runConfig,
          std::make_unique<DefaultElectrostaticSourceModel>(),
          std::make_unique<DefaultWallInteractionModel>()) {}

TokamakEngine::TokamakEngine(
    const RunConfig& runConfig,
    std::unique_ptr<IElectrostaticSourceModel> sourceModel,
    std::unique_ptr<IWallInteractionModel> wallInteractionModel)
    : runConfig_(runConfig),
      electrostaticSourceModel_((sourceModel != nullptr) ? std::move(sourceModel) : std::make_unique<DefaultElectrostaticSourceModel>()),
      wallInteractionModel_((wallInteractionModel != nullptr) ? std::move(wallInteractionModel) : std::make_unique<DefaultWallInteractionModel>()),
      particles_(runConfig.particleCap),
      grid_(std::make_unique<SpatialGrid>(config_.majorRadius_m + config_.minorRadius_m, 0.2f)),
      fusionAcceptedByRadiusBins_(std::max<std::size_t>(1, runConfig.fusionDiagnosticsRadialBins), 0),
      fusionAcceptedWeightByRadiusBins_(std::max<std::size_t>(1, runConfig.fusionDiagnosticsRadialBins), 0.0),
      magneticFieldRadialMean_T_(std::max<std::size_t>(1, runConfig.magneticFieldRadialBinCount), 0.0),
      magneticFieldRadialSum_T_(std::max<std::size_t>(1, runConfig.magneticFieldRadialBinCount), 0.0),
      magneticFieldRadialCounts_(std::max<std::size_t>(1, runConfig.magneticFieldRadialBinCount), 0),
      magneticFieldRecommendedDt_s_(std::numeric_limits<double>::infinity()) {
    if (runConfig_.seed.has_value()) {
        activeSeed_ = runConfig_.seed.value();
    } else {
        activeSeed_ = static_cast<uint64_t>(std::random_device{}());
    }
    rng_.seed(static_cast<uint32_t>(activeSeed_));

    plasmaCurrentProfile_.kind = runConfig_.plasmaCurrentProfileKind;
    plasmaCurrentProfile_.axisEpsilon_m = runConfig_.plasmaCurrentAxisEpsilon_m;
    plasmaCurrentProfile_.customAxisBlendRadius_m = runConfig_.plasmaCurrentCustomAxisBlendRadius_m;
    plasmaCurrentProfile_.customTable = runConfig_.plasmaCurrentCustomTable;
    electrostaticGeometry_.gridWidth = static_cast<int>(std::max<std::size_t>(2, runConfig_.electrostaticGridBinCount));
    electrostaticGeometry_.offset_m = static_cast<double>(config_.majorRadius_m + config_.minorRadius_m);
    electrostaticGeometry_.cellSize_m =
        (2.0 * electrostaticGeometry_.offset_m) / static_cast<double>(electrostaticGeometry_.gridWidth);
    electrostaticSolveConfig_.boundaryCondition = runConfig_.electrostaticBoundaryCondition;
    electrostaticSolveConfig_.chargeAssignmentScheme = runConfig_.chargeAssignmentScheme;
    electrostaticSolveConfig_.tolerance = runConfig_.electrostaticSolverTolerance;
    electrostaticSolveConfig_.maxIterations = runConfig_.electrostaticSolverMaxIterations;
    electrostaticSolveConfig_.sorOmega = runConfig_.electrostaticSorOmega;
    electrostaticSolveConfig_.neutralizingBackgroundFraction = runConfig_.electrostaticNeutralizingBackgroundFraction;
    electrostaticChargeDensity_CPerM3_.assign(ElectrostaticCellCount(electrostaticGeometry_), 0.0);
    electrostaticPotential_V_.assign(ElectrostaticCellCount(electrostaticGeometry_), 0.0);
    electrostaticField_VPerM_.assign(ElectrostaticCellCount(electrostaticGeometry_), Vec3(0.0f, 0.0f, 0.0f));
    electrostaticSolverResidualStep_ = SolverResidualSnapshot{};
    electrostaticSolverResidualStep_.solverKind = ResidualSolverKind::None;
    electrostaticSolverResidualStep_.status = ResidualStatus::Placeholder;

    switch (runConfig_.scenario) {
        case Scenario::ColdVacuum:
            nbi_.isActive = false;
            break;
        case Scenario::NbiIgnition:
            nbi_.isActive = true;
            nbi_.beamEnergy_keV = 120.0f;
            break;
        case Scenario::MagneticFailure:
            nbi_.isActive = true;
            config_.plasmaCurrent_A = 0.0f;
            break;
    }

    const double referenceDt_s = std::max(1.0e-12, static_cast<double>(runConfig_.timeStep_s));
    nbiPairsPerSecond_ = static_cast<double>(nbi_.particlesPerStep) / referenceDt_s;

    std::uniform_real_distribution<float> angleDist(0.0f, 2.0f * constants::kPi);
    std::uniform_real_distribution<float> radiusDist(0.0f, config_.minorRadius_m * 0.8f);

    for (int i = 0; i < 1000; ++i) {
        const float a1 = angleDist(rng_);
        const float a2 = angleDist(rng_);
        const float rad = radiusDist(rng_);

        const Vec3 pos(
            (config_.majorRadius_m + (rad * std::cos(a2))) * std::cos(a1),
            (config_.majorRadius_m + (rad * std::cos(a2))) * std::sin(a1),
            rad * std::sin(a2));

        const ParticleType type = (i % 2 == 0) ? ParticleType::Deuterium : ParticleType::Tritium;
        const double mass = (type == ParticleType::Deuterium) ? constants::kMassDeuterium_kg : constants::kMassTritium_kg;

        if (!particles_.AddParticle(pos, Vec3(0.0f, 0.0f, 0.0f), mass, constants::kElementaryCharge_C, type)) {
            ++counters_.particleCapHitEvents;
            break;
        }
    }

#ifndef NDEBUG
    ValidateStateDebug();
#endif
}

double TokamakEngine::StartupRampFraction(double sampleTime_s) const {
    if (runConfig_.startupRampDuration_s <= 0.0) {
        return 1.0;
    }
    return SmoothStepUnit(sampleTime_s / runConfig_.startupRampDuration_s);
}

double TokamakEngine::FusionRampFraction(double sampleTime_s) const {
    const double startupFraction = StartupRampFraction(sampleTime_s);
    const double fusionStartTime_s =
        std::max(0.0, runConfig_.startupRampDuration_s) + std::max(0.0, runConfig_.fusionStartDelay_s);

    if (sampleTime_s < fusionStartTime_s) {
        return 0.0;
    }

    if (runConfig_.fusionRampDuration_s <= 0.0) {
        return startupFraction;
    }

    const double fusionOnsetFraction =
        SmoothStepUnit((sampleTime_s - fusionStartTime_s) / runConfig_.fusionRampDuration_s);
    return startupFraction * fusionOnsetFraction;
}

TokamakConfig TokamakEngine::EffectiveTokamakConfig(double sampleTime_s) const {
    TokamakConfig effectiveConfig = config_;
    const float rampFraction = static_cast<float>(StartupRampFraction(sampleTime_s));
    effectiveConfig.toroidalCurrent_A *= rampFraction;
    effectiveConfig.plasmaCurrent_A *= rampFraction;
    return effectiveConfig;
}

NBIConfig TokamakEngine::EffectiveNbiConfig(double sampleTime_s) const {
    NBIConfig effectiveNbi = nbi_;
    const float rampFraction = static_cast<float>(StartupRampFraction(sampleTime_s));
    effectiveNbi.particlesPerStep = static_cast<int>(std::lround(
        static_cast<double>(nbi_.particlesPerStep) * static_cast<double>(rampFraction)));
    effectiveNbi.isActive = effectiveNbi.isActive && (rampFraction > 0.0f);
    return effectiveNbi;
}

RunConfig TokamakEngine::EffectiveCollisionRunConfig(double sampleTime_s) const {
    RunConfig effectiveRunConfig = runConfig_;
    effectiveRunConfig.fusionCrossSectionScale *= FusionRampFraction(sampleTime_s);
    return effectiveRunConfig;
}

Vec3 TokamakEngine::CalculateBField(const Vec3& position) const {
    const MagneticFieldSample sample =
        EvaluateMagneticFieldSample(EffectiveTokamakConfig(time_s_), plasmaCurrentProfile_, position);
    return sample.totalField_T;
}

Vec3 TokamakEngine::CalculateEField(const Vec3& position) const {
    if (runConfig_.electricFieldMode == ElectricFieldMode::Electrostatic) {
        return CalculateElectrostaticEField(position);
    }
    return CalculatePlaceholderEField(position);
}

Vec3 TokamakEngine::CalculatePlaceholderEField(const Vec3& position) const {
    const float R = std::sqrt((position.x * position.x) + (position.y * position.y));
    if (R < 0.001f) {
        return Vec3(0.0f, 0.0f, 0.0f);
    }

    const Vec3 coreCenter(config_.majorRadius_m * (position.x / R), config_.majorRadius_m * (position.y / R), 0.0f);
    const float startupFraction = static_cast<float>(StartupRampFraction(time_s_));
    return (coreCenter - position) * (50.0f * startupFraction);
}

Vec3 TokamakEngine::CalculateElectrostaticEField(const Vec3& position) const {
    return SampleElectricField(
        electrostaticGeometry_,
        electrostaticField_VPerM_,
        position,
        electrostaticSolveConfig_.chargeAssignmentScheme);
}

void TokamakEngine::PushParticles(float dt_s) {
    const TokamakConfig effectiveConfig = EffectiveTokamakConfig(time_s_);
    auto& positions = particles_.MutablePositions();
    auto& velocities = particles_.MutableVelocities();
    auto& weights = particles_.MutableWeights();
    auto& species = particles_.MutableSpecies();
    const auto& masses = particles_.Masses();
    const auto& chargeToMass = particles_.ChargeToMass();
    std::uniform_real_distribution<float> recycleVerticalSpread(-0.25f, 0.25f);

    for (std::size_t i = 0; i < positions.size(); ++i) {
        if (species[i] == ParticleType::Dead) {
            continue;
        }

        const Vec3 pos = positions[i];
        const Vec3 vel = velocities[i];
        const double qOverM = chargeToMass[i];

        const Vec3 E = CalculateEField(pos);
        const MagneticFieldSample magneticSample =
            EvaluateMagneticFieldSample(effectiveConfig, plasmaCurrentProfile_, pos);
        AccumulateMagneticFieldDiagnosticsSample(magneticSample);
        const Vec3 B = magneticSample.totalField_T;

        Vec3 vNew = BorisVelocityStep(vel, E, B, qOverM, dt_s);
        Vec3 posNew = pos + (vNew * dt_s);
        const Vec3 preWallVel = vNew;
        TokamakWallHit wallHitGeometry;
        if (ComputeTokamakWallHit(posNew, config_.majorRadius_m, config_.minorRadius_m, &wallHitGeometry)) {
            ++counters_.wallHitCount;
            ++stepCounters_.wallHitCount;

            const double speedSq = static_cast<double>(Vec3::Dot(preWallVel, preWallVel));
            const double impactEnergy_J =
                0.5 * static_cast<double>(masses[i]) * speedSq * static_cast<double>(weights[i]);
            counters_.wallImpactEnergy_J += impactEnergy_J;
            stepCounters_.wallImpactEnergy_J += impactEnergy_J;

            double wallLossWeight = 0.0;
            const bool recycleEligible =
                (runConfig_.wallBoundaryMode == WallBoundaryMode::Recycle) &&
                (species[i] == ParticleType::Deuterium || species[i] == ParticleType::Tritium);

            if (runConfig_.wallBoundaryMode == WallBoundaryMode::Reflect) {
                ReflectAtTokamakWall(&posNew, &vNew, config_.majorRadius_m, config_.minorRadius_m);
            } else if (runConfig_.wallBoundaryMode == WallBoundaryMode::Absorb || !recycleEligible) {
                wallLossWeight = weights[i];
                weights[i] = 0.0;
                species[i] = ParticleType::Dead;
                posNew = wallHitGeometry.projectedPosition;
                vNew = Vec3(0.0f, 0.0f, 0.0f);
            } else {
                const double recycleFraction = std::max(0.0, std::min(1.0, runConfig_.recycleFraction));
                const double recycledWeight = weights[i] * recycleFraction;
                wallLossWeight = weights[i] - recycledWeight;
                weights[i] = recycledWeight;

                if (recycledWeight <= 1.0e-9) {
                    weights[i] = 0.0;
                    species[i] = ParticleType::Dead;
                    posNew = wallHitGeometry.projectedPosition;
                    vNew = Vec3(0.0f, 0.0f, 0.0f);
                } else {
                    const Vec3 toroidalTangent = ToroidalTangent(wallHitGeometry.coreCenter);
                    const float preSpeed = preWallVel.Magnitude();
                    const float recycleSpeed = std::max(2.0e4f, preSpeed * 0.10f);
                    const float toroidalSpeed = recycleSpeed * 0.35f;
                    const float verticalSpeed = recycleSpeed * recycleVerticalSpread(rng_);
                    Vec3 recycledVelocity =
                        (wallHitGeometry.normal * -recycleSpeed) +
                        (toroidalTangent * toroidalSpeed) +
                        Vec3(0.0f, 0.0f, verticalSpeed);
                    if (recycledVelocity.Magnitude() <= 1.0e-6f) {
                        recycledVelocity = wallHitGeometry.normal * -recycleSpeed;
                    }

                    posNew = wallHitGeometry.projectedPosition + (wallHitGeometry.normal * -0.002f);
                    vNew = recycledVelocity;
                }
            }

            counters_.wallLossWeight += wallLossWeight;
            stepCounters_.wallLossWeight += wallLossWeight;
        }

        velocities[i] = vNew;
        positions[i] = posNew;
    }
}

void TokamakEngine::PrepareElectrostaticFieldStep() {
    electrostaticMaxField_VPerM_ = 0.0;
    electrostaticMeanField_VPerM_ = 0.0;

    if (runConfig_.electricFieldMode != ElectricFieldMode::Electrostatic) {
        electrostaticSolverResidualStep_ = SolverResidualSnapshot{};
        electrostaticSolverResidualStep_.solverKind = ResidualSolverKind::None;
        electrostaticSolverResidualStep_.status = ResidualStatus::Placeholder;
        return;
    }

    electrostaticSolveConfig_.boundaryCondition = wallInteractionModel_->ResolveBoundaryCondition(runConfig_, config_);
    electrostaticSolveConfig_.chargeAssignmentScheme = runConfig_.chargeAssignmentScheme;
    const auto& charges = particles_.Charges();
    const auto& weights = particles_.Weights();
    const std::size_t weightedChargeCount = std::min(charges.size(), weights.size());
    electrostaticEffectiveCharges_C_.resize(weightedChargeCount);
    for (std::size_t i = 0; i < weightedChargeCount; ++i) {
        electrostaticEffectiveCharges_C_[i] = charges[i] * weights[i];
    }

    DepositChargeDensity(
        particles_.Positions(),
        electrostaticEffectiveCharges_C_,
        1.0,
        electrostaticGeometry_,
        electrostaticSolveConfig_.chargeAssignmentScheme,
        &electrostaticChargeDensity_CPerM3_);

    const double cellVolume_m3 = electrostaticGeometry_.cellSize_m * electrostaticGeometry_.cellSize_m *
                                 electrostaticGeometry_.cellSize_m;
    electrostaticSourceModel_->ApplySourceTerms(
        runConfig_,
        config_,
        particles_,
        cellVolume_m3,
        &electrostaticChargeDensity_CPerM3_);
    ApplyNeutralizingBackground(
        electrostaticSolveConfig_.neutralizingBackgroundFraction,
        &electrostaticChargeDensity_CPerM3_);

    const double dirichletBoundaryPotential_V =
        wallInteractionModel_->DirichletBoundaryPotential_V(runConfig_, config_, 0, 0, 0);
    const ElectrostaticSolveResult solveResult = SolvePoissonSor(
        electrostaticGeometry_,
        electrostaticSolveConfig_,
        dirichletBoundaryPotential_V,
        electrostaticChargeDensity_CPerM3_,
        &electrostaticPotential_V_);
    electrostaticSolverResidualStep_ = solveResult.residual;

    ReconstructElectricField(
        electrostaticGeometry_,
        electrostaticPotential_V_,
        &electrostaticField_VPerM_);

    double magnitudeSum = 0.0;
    for (const Vec3& field : electrostaticField_VPerM_) {
        const double magnitude = static_cast<double>(field.Magnitude());
        if (!std::isfinite(magnitude)) {
            continue;
        }
        electrostaticMaxField_VPerM_ = std::max(electrostaticMaxField_VPerM_, magnitude);
        magnitudeSum += magnitude;
    }
    if (!electrostaticField_VPerM_.empty()) {
        electrostaticMeanField_VPerM_ = magnitudeSum / static_cast<double>(electrostaticField_VPerM_.size());
    }
}

void TokamakEngine::ResetMagneticFieldDiagnosticsStep() {
    magneticFieldStepMax_T_ = 0.0;
    std::fill(magneticFieldRadialSum_T_.begin(), magneticFieldRadialSum_T_.end(), 0.0);
    std::fill(magneticFieldRadialCounts_.begin(), magneticFieldRadialCounts_.end(), 0);
}

void TokamakEngine::AccumulateMagneticFieldDiagnosticsSample(const MagneticFieldSample& sample) {
    magneticFieldStepMax_T_ = std::max(magneticFieldStepMax_T_, static_cast<double>(sample.totalMagnitude_T));
    if (magneticFieldRadialSum_T_.empty()) {
        return;
    }

    const float clampedRadius = std::max(0.0f, std::min(sample.normalizedMinorRadius, 1.0f));
    const float scaled = clampedRadius * static_cast<float>(magneticFieldRadialSum_T_.size());
    std::size_t binIndex = static_cast<std::size_t>(scaled);
    if (binIndex >= magneticFieldRadialSum_T_.size()) {
        binIndex = magneticFieldRadialSum_T_.size() - 1;
    }

    magneticFieldRadialSum_T_[binIndex] += static_cast<double>(sample.totalMagnitude_T);
    ++magneticFieldRadialCounts_[binIndex];
}

void TokamakEngine::FinalizeMagneticFieldDiagnosticsStep() {
    for (std::size_t i = 0; i < magneticFieldRadialMean_T_.size(); ++i) {
        const uint64_t count = magneticFieldRadialCounts_[i];
        magneticFieldRadialMean_T_[i] =
            (count == 0U) ? 0.0 : (magneticFieldRadialSum_T_[i] / static_cast<double>(count));
    }
    magneticFieldRecommendedDt_s_ =
        RecommendDtFromMaxField(magneticFieldStepMax_T_, runConfig_.magneticFieldDtSafetyFraction);
}

void TokamakEngine::RunNBI(float dt_s) {
    const double startupFraction = StartupRampFraction(time_s_);
    const NBIConfig effectiveNbi = EffectiveNbiConfig(time_s_);

    if (!effectiveNbi.isActive) {
        return;
    }
    if (dt_s <= 0.0f) {
        return;
    }

    const double eJoules = static_cast<double>(effectiveNbi.beamEnergy_keV) * 1000.0 * constants::kElementaryCharge_C;
    const float speed = static_cast<float>(std::sqrt((2.0 * eJoules) / constants::kMassDeuterium_kg));
    const Vec3 beamVel = effectiveNbi.injectionNormal * speed;
    std::uniform_real_distribution<float> spreadDist(-0.1f, 0.1f);

    nbiPairAccumulator_ += nbiPairsPerSecond_ * startupFraction * static_cast<double>(dt_s);
    const int pairCount = static_cast<int>(std::floor(nbiPairAccumulator_));
    nbiPairAccumulator_ -= static_cast<double>(pairCount);

    for (int i = 0; i < pairCount; ++i) {
        if (!particles_.CanInsert(2)) {
            ++counters_.rejectedInjectionPairs;
            ++counters_.particleCapHitEvents;
            ++stepCounters_.rejectedInjectionPairs;
            ++stepCounters_.particleCapHitEvents;
            break;
        }

        Vec3 spawnPos = effectiveNbi.injectorPos + (effectiveNbi.injectionNormal * config_.minorRadius_m) +
                        Vec3(spreadDist(rng_), spreadDist(rng_), spreadDist(rng_));

        const bool dAdded = particles_.AddParticle(
            spawnPos,
            beamVel,
            constants::kMassDeuterium_kg,
            constants::kElementaryCharge_C,
            ParticleType::Deuterium);

        spawnPos = spawnPos + Vec3(spreadDist(rng_), spreadDist(rng_), spreadDist(rng_));

        const bool tAdded = particles_.AddParticle(
            spawnPos,
            beamVel,
            constants::kMassTritium_kg,
            constants::kElementaryCharge_C,
            ParticleType::Tritium);

        if (!dAdded || !tAdded) {
            ++counters_.rejectedInjectionPairs;
            ++counters_.particleCapHitEvents;
            ++stepCounters_.rejectedInjectionPairs;
            ++stepCounters_.particleCapHitEvents;
            break;
        }

        const double beamSpeedSq = static_cast<double>(Vec3::Dot(beamVel, beamVel));
        const auto& insertedWeights = particles_.Weights();
        const std::size_t insertedCount = insertedWeights.size();
        budget_.beamInjected_J +=
            0.5 * constants::kMassDeuterium_kg * beamSpeedSq * insertedWeights[insertedCount - 2];
        budget_.beamInjected_J +=
            0.5 * constants::kMassTritium_kg * beamSpeedSq * insertedWeights[insertedCount - 1];
    }
}

void TokamakEngine::SortGrid() {
    const float gridOffset_m = config_.majorRadius_m + config_.minorRadius_m;
    const uint64_t clampCount = SortParticlesIntoGrid(particles_, *grid_, gridOffset_m);
    counters_.outOfDomainCellClampEvents += clampCount;
    stepCounters_.outOfDomainCellClampEvents += clampCount;
}

void TokamakEngine::SelectCollisionEvents(float dt_s) {
    const RunConfig effectiveRunConfig = EffectiveCollisionRunConfig(time_s_ + (0.5 * static_cast<double>(dt_s)));
    const auto summary =
        tokamak::SelectCollisionEvents(effectiveRunConfig, particles_, *grid_, rng_, dt_s, pendingFusionEvents_);
    counters_.fusionAttempts += summary.fusionAttempts;
    stepCounters_.fusionAttempts += summary.fusionAttempts;
    counters_.fusionKineticsSamples += summary.fusionKineticsSamples;
    stepCounters_.fusionKineticsSamples += summary.fusionKineticsSamples;
    counters_.fusionWeightAttempted += summary.fusionWeightAttempted;
    stepCounters_.fusionWeightAttempted += summary.fusionWeightAttempted;
    counters_.fusionSigmaSum_m2 += summary.fusionSigmaSum_m2;
    stepCounters_.fusionSigmaSum_m2 += summary.fusionSigmaSum_m2;
    counters_.fusionProbabilitySum += summary.fusionProbabilitySum;
    stepCounters_.fusionProbabilitySum += summary.fusionProbabilitySum;
    counters_.fusionRelativeSpeedSum_mPerS += summary.fusionRelativeSpeedSum_mPerS;
    stepCounters_.fusionRelativeSpeedSum_mPerS += summary.fusionRelativeSpeedSum_mPerS;
    if (summary.maxReactionsInCell > counters_.maxReactionsInCell) {
        counters_.maxReactionsInCell = summary.maxReactionsInCell;
    }
    if (summary.maxReactionsInCell > stepCounters_.maxReactionsInCell) {
        stepCounters_.maxReactionsInCell = summary.maxReactionsInCell;
    }
}

void TokamakEngine::ApplyCollisionEvents() {
    const RuntimeCounters before = counters_;
    tokamak::ApplyCollisionEvents(
        particles_,
        pendingFusionEvents_,
        counters_,
        budget_,
        &acceptedFusionPositions_,
        &acceptedFusionWeights_);
    stepCounters_.particleCapHitEvents += (counters_.particleCapHitEvents - before.particleCapHitEvents);
    stepCounters_.rejectedFusionAsh += (counters_.rejectedFusionAsh - before.rejectedFusionAsh);
    stepCounters_.fusionAccepted += (counters_.fusionAccepted - before.fusionAccepted);
    stepCounters_.fusionWeightAccepted += (counters_.fusionWeightAccepted - before.fusionWeightAccepted);
    stepCounters_.fuelWeightConsumedD += (counters_.fuelWeightConsumedD - before.fuelWeightConsumedD);
    stepCounters_.fuelWeightConsumedT += (counters_.fuelWeightConsumedT - before.fuelWeightConsumedT);
    stepCounters_.ashWeightProducedHe += (counters_.ashWeightProducedHe - before.ashWeightProducedHe);

    for (std::size_t i = 0; i < acceptedFusionPositions_.size(); ++i) {
        const Vec3& position = acceptedFusionPositions_[i];
        const std::size_t binIndex = FusionRadiusBinIndex(position);
        if (binIndex < fusionAcceptedByRadiusBins_.size()) {
            ++fusionAcceptedByRadiusBins_[binIndex];
        }
        if (binIndex < fusionAcceptedWeightByRadiusBins_.size() && i < acceptedFusionWeights_.size()) {
            fusionAcceptedWeightByRadiusBins_[binIndex] += acceptedFusionWeights_[i];
        }
    }
}

void TokamakEngine::Step(float dt_s) {
#ifndef NDEBUG
    ValidateStateDebug();
#endif

    stepCounters_ = RuntimeCounters{};
    ResetMagneticFieldDiagnosticsStep();
    RunNBI(dt_s);
    PrepareElectrostaticFieldStep();
    PushParticles(dt_s);
    FinalizeMagneticFieldDiagnosticsStep();
    SortGrid();
    SelectCollisionEvents(dt_s);
    ApplyCollisionEvents();
    particles_.Compact();

    time_s_ += dt_s;

#ifndef NDEBUG
    ValidateStateDebug();
#endif
}

SpeciesCounts TokamakEngine::CountSpecies(
    double* totalKineticEnergy_J,
    double* totalCharge_C,
    double* totalMacroWeight) const {
    SpeciesCounts counts;
    double totalKineticEnergy = 0.0;
    double totalCharge = 0.0;
    double totalWeight = 0.0;

    const auto& velocities = particles_.Velocities();
    const auto& masses = particles_.Masses();
    const auto& charges = particles_.Charges();
    const auto& weights = particles_.Weights();
    const auto& species = particles_.Species();

    for (std::size_t i = 0; i < particles_.Size(); ++i) {
        if (species[i] == ParticleType::Dead) {
            continue;
        }

        if (species[i] == ParticleType::Deuterium) {
            ++counts.deuterium;
        } else if (species[i] == ParticleType::Tritium) {
            ++counts.tritium;
        } else if (species[i] == ParticleType::Helium) {
            ++counts.helium;
        }

        const double speedSq = static_cast<double>(Vec3::Dot(velocities[i], velocities[i]));
        totalKineticEnergy += 0.5 * masses[i] * speedSq * weights[i];
        totalCharge += charges[i] * weights[i];
        totalWeight += weights[i];
    }

    if (totalKineticEnergy_J != nullptr) {
        *totalKineticEnergy_J = totalKineticEnergy;
    }
    if (totalCharge_C != nullptr) {
        *totalCharge_C = totalCharge;
    }
    if (totalMacroWeight != nullptr) {
        *totalMacroWeight = totalWeight;
    }

    return counts;
}

TelemetrySnapshot TokamakEngine::Snapshot(int step) const {
    TelemetrySnapshot snapshot;
    snapshot.step = step;
    snapshot.time_s = time_s_;
    snapshot.activeSeed = activeSeed_;

    double totalKineticEnergy_J = 0.0;
    double totalCharge_C = 0.0;
    double totalMacroWeight = 0.0;
    snapshot.species = CountSpecies(&totalKineticEnergy_J, &totalCharge_C, &totalMacroWeight);

    if (totalMacroWeight > 0.0) {
        snapshot.avgEnergy_keV =
            (totalKineticEnergy_J / totalMacroWeight) /
            (1000.0 * constants::kElementaryCharge_C);
    }

    snapshot.fusionEvents = particles_.fusionCountTotal;
    snapshot.stepCounters = stepCounters_;
    snapshot.counters = counters_;
    snapshot.budget = budget_;
    snapshot.budget.kinetic_J = totalKineticEnergy_J;
    snapshot.budget.totalCharge_C = totalCharge_C;
    snapshot.magneticField.maxField_T = magneticFieldStepMax_T_;
    snapshot.magneticField.recommendedDt_s = magneticFieldRecommendedDt_s_;
    snapshot.magneticField.radialMeanField_T = magneticFieldRadialMean_T_;
    snapshot.magneticField.radialSampleCounts = magneticFieldRadialCounts_;
    snapshot.electrostaticField.maxElectricField_VPerM = electrostaticMaxField_VPerM_;
    snapshot.electrostaticField.meanElectricField_VPerM = electrostaticMeanField_VPerM_;
    snapshot.electrostaticField.solveIterations = electrostaticSolverResidualStep_.iterations;
    snapshot.electrostaticField.solveConverged = electrostaticSolverResidualStep_.converged;
    snapshot.solverResidual = electrostaticSolverResidualStep_;

    return snapshot;
}

bool TokamakEngine::HasFiniteState() const {
    if (!particles_.IsFiniteState()) {
        return false;
    }

    const TelemetrySnapshot snapshot = Snapshot(0);
    return std::isfinite(snapshot.time_s) && std::isfinite(snapshot.avgEnergy_keV) &&
           std::isfinite(snapshot.budget.kinetic_J) && std::isfinite(snapshot.budget.beamInjected_J) &&
           std::isfinite(snapshot.budget.fusionAlphaInjected_J) && std::isfinite(snapshot.budget.totalCharge_C) &&
           std::isfinite(snapshot.electrostaticField.maxElectricField_VPerM) &&
           std::isfinite(snapshot.electrostaticField.meanElectricField_VPerM);
}

void TokamakEngine::ValidateStateDebug() const {
    assert(particles_.IsArrayLengthConsistent());
    assert(particles_.IsFiniteState());
    assert(grid_ != nullptr);

    const float gridOffset_m = config_.majorRadius_m + config_.minorRadius_m;
    const auto& positions = particles_.Positions();
    for (std::size_t i = 0; i < positions.size(); ++i) {
        const int cellIndex = grid_->GetCellIndex(positions[i], gridOffset_m);
        assert(grid_->IsValidCellIndex(cellIndex));
    }
}

std::size_t TokamakEngine::FusionRadiusBinIndex(const Vec3& position) const {
    if (fusionAcceptedByRadiusBins_.empty() || config_.minorRadius_m <= 0.0f) {
        return 0;
    }

    const float majorRadius = std::sqrt((position.x * position.x) + (position.y * position.y));
    const float minorRadius =
        std::sqrt(std::pow(majorRadius - config_.majorRadius_m, 2.0f) + (position.z * position.z));
    const float clampedMinorRadius = std::max(0.0f, std::min(minorRadius, config_.minorRadius_m));
    const float normalized = clampedMinorRadius / config_.minorRadius_m;
    const float scaled = normalized * static_cast<float>(fusionAcceptedByRadiusBins_.size());
    const std::size_t index = static_cast<std::size_t>(scaled);
    return (index >= fusionAcceptedByRadiusBins_.size()) ? (fusionAcceptedByRadiusBins_.size() - 1) : index;
}

}  // namespace tokamak
