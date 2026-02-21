#pragma once

#include <memory>
#include <random>
#include <vector>

#include "tokamak/collision.hpp"
#include "tokamak/electrostatic_field.hpp"
#include "tokamak/electrostatic_models.hpp"
#include "tokamak/magnetic_field.hpp"

namespace tokamak {

class TokamakEngine {
public:
    explicit TokamakEngine(const RunConfig& runConfig);
    TokamakEngine(
        const RunConfig& runConfig,
        std::unique_ptr<IElectrostaticSourceModel> sourceModel,
        std::unique_ptr<IWallInteractionModel> wallInteractionModel);

    void Step(float dt_s);
    void PrintTelemetry(int step) const;

    TelemetrySnapshot Snapshot(int step) const;
    bool HasFiniteState() const;

    const RuntimeCounters& Counters() const { return counters_; }
    const EnergyChargeBudget& Budget() const { return budget_; }

    const ParticleSystem& Particles() const { return particles_; }
    ParticleSystem& MutableParticles() { return particles_; }

    const SpatialGrid& Grid() const { return *grid_; }
    SpatialGrid& MutableGrid() { return *grid_; }

    const TokamakConfig& Config() const { return config_; }
    const NBIConfig& BeamConfig() const { return nbi_; }
    const PlasmaCurrentProfileConfig& PlasmaCurrentProfile() const { return plasmaCurrentProfile_; }
    ElectricFieldMode CurrentElectricFieldMode() const { return runConfig_.electricFieldMode; }
    const std::vector<double>& MagneticFieldRadialMean() const { return magneticFieldRadialMean_T_; }
    const std::vector<uint64_t>& MagneticFieldRadialCounts() const { return magneticFieldRadialCounts_; }
    uint64_t ActiveSeed() const { return activeSeed_; }
    double SimTimeSeconds() const { return time_s_; }
    const std::vector<uint64_t>& FusionAcceptedByRadiusBins() const { return fusionAcceptedByRadiusBins_; }

private:
    Vec3 CalculateBField(const Vec3& position) const;
    Vec3 CalculateEField(const Vec3& position) const;
    Vec3 CalculatePlaceholderEField(const Vec3& position) const;
    Vec3 CalculateElectrostaticEField(const Vec3& position) const;

    void PushParticles(float dt_s);
    void PrepareElectrostaticFieldStep();
    void RunNBI(float dt_s);
    void SortGrid();
    void SelectCollisionEvents(float dt_s);
    void ApplyCollisionEvents();
    SpeciesCounts CountSpecies(double* totalKineticEnergy_J, double* totalCharge_C) const;
    std::size_t FusionRadiusBinIndex(const Vec3& position) const;
    void ResetMagneticFieldDiagnosticsStep();
    void AccumulateMagneticFieldDiagnosticsSample(const MagneticFieldSample& sample);
    void FinalizeMagneticFieldDiagnosticsStep();

    void ValidateStateDebug() const;

    RunConfig runConfig_;
    TokamakConfig config_;
    NBIConfig nbi_;
    PlasmaCurrentProfileConfig plasmaCurrentProfile_;
    std::unique_ptr<IElectrostaticSourceModel> electrostaticSourceModel_;
    std::unique_ptr<IWallInteractionModel> wallInteractionModel_;

    ParticleSystem particles_;
    std::unique_ptr<SpatialGrid> grid_;

    std::mt19937 rng_;
    uint64_t activeSeed_ = 0;
    double time_s_ = 0.0;
    double nbiPairsPerSecond_ = 0.0;
    double nbiPairAccumulator_ = 0.0;

    RuntimeCounters counters_;
    RuntimeCounters stepCounters_;
    EnergyChargeBudget budget_;

    std::vector<PendingFusionEvent> pendingFusionEvents_;
    std::vector<Vec3> acceptedFusionPositions_;
    std::vector<uint64_t> fusionAcceptedByRadiusBins_;

    std::vector<double> magneticFieldRadialMean_T_;
    std::vector<double> magneticFieldRadialSum_T_;
    std::vector<uint64_t> magneticFieldRadialCounts_;
    double magneticFieldStepMax_T_ = 0.0;
    double magneticFieldRecommendedDt_s_ = 0.0;

    ElectrostaticMeshGeometry electrostaticGeometry_;
    ElectrostaticSolveConfig electrostaticSolveConfig_;
    std::vector<double> electrostaticChargeDensity_CPerM3_;
    std::vector<float> electrostaticEffectiveCharges_C_;
    std::vector<double> electrostaticPotential_V_;
    std::vector<Vec3> electrostaticField_VPerM_;
    SolverResidualSnapshot electrostaticSolverResidualStep_;
    double electrostaticMaxField_VPerM_ = 0.0;
    double electrostaticMeanField_VPerM_ = 0.0;
};

}  // namespace tokamak
