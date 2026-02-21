#include "tokamak/engine.hpp"

#include <iomanip>
#include <iostream>

namespace tokamak {
namespace {

const char* ResidualStatusText(ResidualStatus status) {
    switch (status) {
        case ResidualStatus::Placeholder:
            return "placeholder";
        case ResidualStatus::Measured:
            return "measured";
        case ResidualStatus::Unavailable:
            return "unavailable";
        case ResidualStatus::Failed:
            return "failed";
    }
    return "unknown";
}

}  // namespace

void TokamakEngine::PrintTelemetry(int step) const {
    const TelemetrySnapshot snapshot = Snapshot(step);

    std::cout << std::fixed << std::setprecision(2)
              << "[Step " << std::setw(5) << snapshot.step << " | Time " << (snapshot.time_s * 1000.0) << " ms] "
              << "TotalIons: " << std::setw(6) << snapshot.species.AliveCount() << " | "
              << "D: " << snapshot.species.deuterium << " T: " << snapshot.species.tritium
              << " He: " << snapshot.species.helium << " | "
              << "AvgE: " << snapshot.avgEnergy_keV << " keV | "
              << "FusionEvents: " << snapshot.fusionEvents << " | "
              << "StepCtrs cap-hit: " << snapshot.stepCounters.particleCapHitEvents
              << " rejected-injection: " << snapshot.stepCounters.rejectedInjectionPairs
              << " rejected-ash: " << snapshot.stepCounters.rejectedFusionAsh
              << " out-of-domain-clamp: " << snapshot.stepCounters.outOfDomainCellClampEvents
              << " fusion-attempts: " << snapshot.stepCounters.fusionAttempts
              << " fusion-accepted: " << snapshot.stepCounters.fusionAccepted
              << " fusion-w-attempted: " << snapshot.stepCounters.fusionWeightAttempted
              << " fusion-w-accepted: " << snapshot.stepCounters.fusionWeightAccepted
              << " wall-hits: " << snapshot.stepCounters.wallHitCount
              << " max-cell-reaction: " << snapshot.stepCounters.maxReactionsInCell << " | "
              << "Ctrs cap-hit: " << snapshot.counters.particleCapHitEvents
              << " rejected-injection: " << snapshot.counters.rejectedInjectionPairs
              << " rejected-ash: " << snapshot.counters.rejectedFusionAsh
              << " out-of-domain-clamp: " << snapshot.counters.outOfDomainCellClampEvents
              << " fusion-attempts: " << snapshot.counters.fusionAttempts
              << " fusion-accepted: " << snapshot.counters.fusionAccepted
              << " fusion-w-attempted: " << snapshot.counters.fusionWeightAttempted
              << " fusion-w-accepted: " << snapshot.counters.fusionWeightAccepted
              << " wall-hits: " << snapshot.counters.wallHitCount
              << " max-cell-reaction: " << snapshot.counters.maxReactionsInCell
              << std::scientific << std::setprecision(6)
              << " | Bmax_T: " << snapshot.magneticField.maxField_T
              << " dt_gyro_recommended_s: " << snapshot.magneticField.recommendedDt_s
              << " | Emax_V_per_m: " << snapshot.electrostaticField.maxElectricField_VPerM
              << " Emean_V_per_m: " << snapshot.electrostaticField.meanElectricField_VPerM
              << " solver_iters: " << snapshot.electrostaticField.solveIterations
              << " solver_converged: " << (snapshot.electrostaticField.solveConverged ? "true" : "false")
              << " residual_l2: " << snapshot.solverResidual.residualL2
              << " residual_status: " << ResidualStatusText(snapshot.solverResidual.status)
              << " | Energy kinetic_J: " << snapshot.budget.kinetic_J
              << " beam_J: " << snapshot.budget.beamInjected_J
              << " fusion_alpha_J: " << snapshot.budget.fusionAlphaInjected_J
              << " total_charge_C: " << snapshot.budget.totalCharge_C << '\n'
              << std::defaultfloat;
}

}  // namespace tokamak
