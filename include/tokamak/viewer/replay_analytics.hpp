#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

#include "tokamak/viewer/replay_manifest.hpp"

namespace tokamak::viewer {

struct ReplayRadialProfileBin {
    int step = 0;
    double time_s = 0.0;
    int binIndex = 0;
    double rInner_m = 0.0;
    double rOuter_m = 0.0;
    double rCenter_m = 0.0;
    uint64_t ionCount = 0;
    double macroWeight = 0.0;
    double shellVolume_m3 = 0.0;
    double density_m3 = 0.0;
    double avgIonEnergy_keV = 0.0;
    uint64_t fusionEventsCumulative = 0;
    double fusionRate_m3_s = 0.0;
    bool fusionRatePlaceholder = false;
};

struct ReplayMagneticFieldBin {
    int step = 0;
    double time_s = 0.0;
    int binIndex = 0;
    double rInner_m = 0.0;
    double rOuter_m = 0.0;
    double rCenter_m = 0.0;
    double meanField_T = 0.0;
    uint64_t sampleCount = 0;
    double stepMaxField_T = 0.0;
    double recommendedDt_s = 0.0;
    std::string profileKind;
};

struct ReplayMagneticFieldStepSummary {
    int step = 0;
    double time_s = 0.0;
    double stepMaxField_T = 0.0;
    double recommendedDt_s = 0.0;
    std::string profileKind;
};

struct ReplayElectrostaticPoint {
    int step = 0;
    double time_s = 0.0;
    std::string electricFieldMode;
    std::string boundaryCondition;
    std::string chargeAssignment;
    double maxElectricField_VPerM = 0.0;
    double meanElectricField_VPerM = 0.0;
    uint32_t solverIterations = 0;
    bool solverConverged = false;
    double residualL2 = 0.0;
};

struct ReplaySpeedHistogramBin {
    int step = 0;
    double time_s = 0.0;
    int binIndex = 0;
    double speedMin_mPerS = 0.0;
    double speedMax_mPerS = 0.0;
    uint64_t count = 0;
    uint64_t totalSamples = 0;
};

struct ReplayPitchHistogramBin {
    int step = 0;
    double time_s = 0.0;
    int binIndex = 0;
    double pitchMin_deg = 0.0;
    double pitchMax_deg = 0.0;
    uint64_t count = 0;
    uint64_t totalSamples = 0;
    uint64_t invalidSamples = 0;
};

struct ReplaySolverResidualPoint {
    int step = 0;
    double time_s = 0.0;
    bool residualAvailable = false;
    double residualL2 = 0.0;
    std::string solverName;
    std::string status;
    uint32_t iterations = 0;
    bool converged = false;
    double tolerance = 0.0;
    std::string note;
};

struct ReplayFusionReactivityPoint {
    int step = 0;
    double time_s = 0.0;
    std::string reactivityModel;
    double crossSectionScale = 0.0;
    double probabilityClamp = 0.0;
    double minEnergy_keV = 0.0;
    uint64_t fusionAttemptsStep = 0;
    uint64_t fusionAcceptedStep = 0;
    double fusionWeightAttemptedStep = 0.0;
    double fusionWeightAcceptedStep = 0.0;
    double fuelWeightConsumedDStep = 0.0;
    double fuelWeightConsumedTStep = 0.0;
    double ashWeightProducedHeStep = 0.0;
    double avgSigma_m2Step = 0.0;
    double avgProbabilityStep = 0.0;
    double avgRelativeSpeed_mPerSStep = 0.0;
    uint64_t fusionAttemptsTotal = 0;
    uint64_t fusionAcceptedTotal = 0;
};

struct ReplayWallInteractionPoint {
    int step = 0;
    double time_s = 0.0;
    std::string wallMode;
    double recycleFraction = 0.0;
    uint64_t wallHitCountStep = 0;
    double wallImpactEnergy_JStep = 0.0;
    double wallLossWeightStep = 0.0;
    uint64_t wallHitCountTotal = 0;
    double wallImpactEnergy_JTotal = 0.0;
    double wallLossWeightTotal = 0.0;
};

class ReplayAnalytics {
public:
    bool LoadFromManifest(const ReplayManifest& manifest, std::string* errorOut);
    void Clear();

    bool HasAnyData() const;

    const std::vector<ReplayRadialProfileBin>* RadialProfileForStep(int step) const;
    const std::vector<ReplayMagneticFieldBin>* MagneticFieldBinsForStep(int step) const;
    const ReplayElectrostaticPoint* ElectrostaticForStep(int step) const;
    const std::vector<ReplaySpeedHistogramBin>* SpeedHistogramForStep(int step) const;
    const std::vector<ReplayPitchHistogramBin>* PitchHistogramForStep(int step) const;
    const ReplaySolverResidualPoint* SolverResidualForStep(int step) const;
    const ReplayFusionReactivityPoint* FusionReactivityForStep(int step) const;
    const ReplayWallInteractionPoint* WallInteractionForStep(int step) const;

    const std::vector<ReplayMagneticFieldStepSummary>& MagneticFieldSeries() const { return magneticFieldSeries_; }
    const std::vector<ReplayElectrostaticPoint>& ElectrostaticSeries() const { return electrostaticSeries_; }
    const std::vector<ReplaySolverResidualPoint>& SolverResidualSeries() const { return solverResidualSeries_; }
    const std::vector<ReplayFusionReactivityPoint>& FusionReactivitySeries() const { return fusionReactivitySeries_; }
    const std::vector<ReplayWallInteractionPoint>& WallInteractionSeries() const { return wallInteractionSeries_; }

private:
    template <typename T>
    const T* LookupPoint(const std::unordered_map<int, std::size_t>& byStep, const std::vector<T>& rows, int step) const {
        const auto it = byStep.find(step);
        if (it == byStep.end() || it->second >= rows.size()) {
            return nullptr;
        }
        return &rows[it->second];
    }

    template <typename T>
    const std::vector<T>* LookupBins(const std::unordered_map<int, std::vector<T>>& byStep, int step) const {
        const auto it = byStep.find(step);
        return (it == byStep.end()) ? nullptr : &it->second;
    }

    std::vector<ReplayRadialProfileBin> radialProfiles_;
    std::unordered_map<int, std::vector<ReplayRadialProfileBin>> radialProfilesByStep_;

    std::vector<ReplayMagneticFieldBin> magneticFieldBins_;
    std::unordered_map<int, std::vector<ReplayMagneticFieldBin>> magneticFieldBinsByStep_;
    std::vector<ReplayMagneticFieldStepSummary> magneticFieldSeries_;
    std::unordered_map<int, std::size_t> magneticFieldSeriesByStep_;

    std::vector<ReplayElectrostaticPoint> electrostaticSeries_;
    std::unordered_map<int, std::size_t> electrostaticSeriesByStep_;

    std::vector<ReplaySpeedHistogramBin> speedHistogramBins_;
    std::unordered_map<int, std::vector<ReplaySpeedHistogramBin>> speedHistogramBinsByStep_;

    std::vector<ReplayPitchHistogramBin> pitchHistogramBins_;
    std::unordered_map<int, std::vector<ReplayPitchHistogramBin>> pitchHistogramBinsByStep_;

    std::vector<ReplaySolverResidualPoint> solverResidualSeries_;
    std::unordered_map<int, std::size_t> solverResidualSeriesByStep_;

    std::vector<ReplayFusionReactivityPoint> fusionReactivitySeries_;
    std::unordered_map<int, std::size_t> fusionReactivitySeriesByStep_;

    std::vector<ReplayWallInteractionPoint> wallInteractionSeries_;
    std::unordered_map<int, std::size_t> wallInteractionSeriesByStep_;
};

}  // namespace tokamak::viewer
