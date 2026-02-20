#pragma once

#include <cstdint>
#include <limits>
#include <vector>

namespace tokamak {

struct RuntimeCounters {
    uint64_t particleCapHitEvents = 0;
    uint64_t rejectedInjectionPairs = 0;
    uint64_t rejectedFusionAsh = 0;
    uint64_t outOfDomainCellClampEvents = 0;
    uint64_t fusionAttempts = 0;
    uint64_t fusionAccepted = 0;
    uint32_t maxReactionsInCell = 0;
};

struct EnergyChargeBudget {
    double kinetic_J = 0.0;
    double beamInjected_J = 0.0;
    double fusionAlphaInjected_J = 0.0;
    double totalCharge_C = 0.0;
};

struct SpeciesCounts {
    uint32_t deuterium = 0;
    uint32_t tritium = 0;
    uint32_t helium = 0;

    uint32_t AliveCount() const { return deuterium + tritium + helium; }
};

enum class ResidualSolverKind : uint8_t {
    None = 0,
    Sor = 1,
    Other = 255,
};

enum class ResidualStatus : uint8_t {
    Placeholder = 0,
    Measured = 1,
    Unavailable = 2,
    Failed = 3,
};

struct SolverResidualSnapshot {
    bool residualAvailable = false;
    double residualL2 = std::numeric_limits<double>::quiet_NaN();
    ResidualSolverKind solverKind = ResidualSolverKind::None;
    ResidualStatus status = ResidualStatus::Placeholder;
    uint32_t iterations = 0;
    bool converged = false;
    double tolerance = std::numeric_limits<double>::quiet_NaN();
};

struct MagneticFieldDiagnostics {
    double maxField_T = 0.0;
    double recommendedDt_s = 0.0;
    std::vector<double> radialMeanField_T;
    std::vector<uint64_t> radialSampleCounts;
};

struct ElectrostaticDiagnostics {
    double maxElectricField_VPerM = 0.0;
    double meanElectricField_VPerM = 0.0;
    uint32_t solveIterations = 0;
    bool solveConverged = false;
};

struct TelemetrySnapshot {
    int step = 0;
    double time_s = 0.0;
    uint64_t activeSeed = 0;
    SpeciesCounts species;
    double avgEnergy_keV = 0.0;
    uint64_t fusionEvents = 0;
    RuntimeCounters stepCounters;
    RuntimeCounters counters;
    EnergyChargeBudget budget;
    MagneticFieldDiagnostics magneticField;
    ElectrostaticDiagnostics electrostaticField;
    SolverResidualSnapshot solverResidual;
};

}  // namespace tokamak
