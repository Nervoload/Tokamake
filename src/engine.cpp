#include "tokamak/engine.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>

#include "tokamak/particle_push.hpp"

namespace tokamak {

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
        const float mass = (type == ParticleType::Deuterium) ? constants::kMassDeuterium_kg : constants::kMassTritium_kg;

        if (!particles_.AddParticle(pos, Vec3(0.0f, 0.0f, 0.0f), mass, constants::kElementaryCharge_C, type)) {
            ++counters_.particleCapHitEvents;
            break;
        }
    }

#ifndef NDEBUG
    ValidateStateDebug();
#endif
}

Vec3 TokamakEngine::CalculateBField(const Vec3& position) const {
    const MagneticFieldSample sample = EvaluateMagneticFieldSample(config_, plasmaCurrentProfile_, position);
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
    return (coreCenter - position) * 50.0f;
}

Vec3 TokamakEngine::CalculateElectrostaticEField(const Vec3& position) const {
    return SampleElectricField(
        electrostaticGeometry_,
        electrostaticField_VPerM_,
        position,
        electrostaticSolveConfig_.chargeAssignmentScheme);
}

void TokamakEngine::PushParticles(float dt_s) {
    auto& positions = particles_.MutablePositions();
    auto& velocities = particles_.MutableVelocities();
    const auto& masses = particles_.Masses();
    const auto& weights = particles_.Weights();
    const auto& chargeToMass = particles_.ChargeToMass();

    for (std::size_t i = 0; i < positions.size(); ++i) {
        const Vec3 pos = positions[i];
        const Vec3 vel = velocities[i];
        const float qOverM = chargeToMass[i];

        const Vec3 E = CalculateEField(pos);
        const MagneticFieldSample magneticSample =
            EvaluateMagneticFieldSample(config_, plasmaCurrentProfile_, pos);
        AccumulateMagneticFieldDiagnosticsSample(magneticSample);
        const Vec3 B = magneticSample.totalField_T;

        Vec3 vNew = BorisVelocityStep(vel, E, B, qOverM, dt_s);
        Vec3 posNew = pos + (vNew * dt_s);
        const Vec3 preWallPos = posNew;
        const Vec3 preWallVel = vNew;

        // M5 seam: keep current behavior reflective while collecting wall-impact counters.
        ReflectAtTokamakWall(&posNew, &vNew, config_.majorRadius_m, config_.minorRadius_m);
        const bool wallHit =
            ((posNew - preWallPos).Magnitude() > 1.0e-8f) || ((vNew - preWallVel).Magnitude() > 1.0e-8f);
        if (wallHit) {
            ++counters_.wallHitCount;
            ++stepCounters_.wallHitCount;

            const double speedSq = static_cast<double>(Vec3::Dot(preWallVel, preWallVel));
            const double impactEnergy_J =
                0.5 * static_cast<double>(masses[i]) * speedSq * static_cast<double>(weights[i]);
            counters_.wallImpactEnergy_J += impactEnergy_J;
            stepCounters_.wallImpactEnergy_J += impactEnergy_J;
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
        1.0f,
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
    if (!nbi_.isActive) {
        return;
    }
    if (dt_s <= 0.0f) {
        return;
    }

    const float eJoules = nbi_.beamEnergy_keV * 1000.0f * constants::kElementaryCharge_C;
    const float speed = std::sqrt((2.0f * eJoules) / constants::kMassDeuterium_kg);
    const Vec3 beamVel = nbi_.injectionNormal * speed;
    std::uniform_real_distribution<float> spreadDist(-0.1f, 0.1f);

    nbiPairAccumulator_ += nbiPairsPerSecond_ * static_cast<double>(dt_s);
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

        Vec3 spawnPos = nbi_.injectorPos + (nbi_.injectionNormal * config_.minorRadius_m) +
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
        budget_.beamInjected_J += 0.5 * static_cast<double>(constants::kMassDeuterium_kg) * beamSpeedSq;
        budget_.beamInjected_J += 0.5 * static_cast<double>(constants::kMassTritium_kg) * beamSpeedSq;
    }
}

void TokamakEngine::SortGrid() {
    const float gridOffset_m = config_.majorRadius_m + config_.minorRadius_m;
    const uint64_t clampCount = SortParticlesIntoGrid(particles_, *grid_, gridOffset_m);
    counters_.outOfDomainCellClampEvents += clampCount;
    stepCounters_.outOfDomainCellClampEvents += clampCount;
}

void TokamakEngine::SelectCollisionEvents(float dt_s) {
    const auto summary = tokamak::SelectCollisionEvents(runConfig_, particles_, *grid_, rng_, dt_s, pendingFusionEvents_);
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
        &acceptedFusionPositions_);
    stepCounters_.particleCapHitEvents += (counters_.particleCapHitEvents - before.particleCapHitEvents);
    stepCounters_.rejectedFusionAsh += (counters_.rejectedFusionAsh - before.rejectedFusionAsh);
    stepCounters_.fusionAccepted += (counters_.fusionAccepted - before.fusionAccepted);
    stepCounters_.fusionWeightAccepted += (counters_.fusionWeightAccepted - before.fusionWeightAccepted);
    stepCounters_.fuelWeightConsumedD += (counters_.fuelWeightConsumedD - before.fuelWeightConsumedD);
    stepCounters_.fuelWeightConsumedT += (counters_.fuelWeightConsumedT - before.fuelWeightConsumedT);
    stepCounters_.ashWeightProducedHe += (counters_.ashWeightProducedHe - before.ashWeightProducedHe);

    for (const Vec3& position : acceptedFusionPositions_) {
        const std::size_t binIndex = FusionRadiusBinIndex(position);
        if (binIndex < fusionAcceptedByRadiusBins_.size()) {
            ++fusionAcceptedByRadiusBins_[binIndex];
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

SpeciesCounts TokamakEngine::CountSpecies(double* totalKineticEnergy_J, double* totalCharge_C) const {
    SpeciesCounts counts;
    double totalKineticEnergy = 0.0;
    double totalCharge = 0.0;

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
        totalKineticEnergy += 0.5 * static_cast<double>(masses[i]) * speedSq;
        totalCharge += static_cast<double>(charges[i]) * static_cast<double>(weights[i]);
    }

    if (totalKineticEnergy_J != nullptr) {
        *totalKineticEnergy_J = totalKineticEnergy;
    }
    if (totalCharge_C != nullptr) {
        *totalCharge_C = totalCharge;
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
    snapshot.species = CountSpecies(&totalKineticEnergy_J, &totalCharge_C);

    const auto alive = snapshot.species.AliveCount();
    if (alive > 0) {
        snapshot.avgEnergy_keV =
            (totalKineticEnergy_J / static_cast<double>(alive)) /
            (1000.0 * static_cast<double>(constants::kElementaryCharge_C));
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
