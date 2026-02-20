#pragma once

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

#include "tokamak/engine.hpp"

namespace tokamak {

constexpr uint32_t kMilestone7OutputSchemaVersion = 2;

struct ArtifactExportConfig {
    std::string outputRootDirectory = "output/runs";
    int metricsEveryNSteps = 100;
    int particleSnapshotEveryNSteps = 100;
    std::size_t maxParticlesPerSnapshot = 50000;
    std::size_t radialBinCount = 32;
    std::size_t speedHistogramBinCount = 40;
    std::size_t pitchHistogramBinCount = 36;
};

class Milestone7ArtifactExporter {
public:
    Milestone7ArtifactExporter() = default;
    ~Milestone7ArtifactExporter();

    bool Start(
        const RunConfig& runConfig,
        const TokamakEngine& engine,
        const ArtifactExportConfig& config);
    bool WriteStep(const TokamakEngine& engine, const TelemetrySnapshot& telemetry, bool forceWrite);
    bool Finish();

    const std::string& LastError() const { return lastError_; }
    const std::string& RunDirectory() const { return runDirectory_; }
    const std::string& RunId() const { return runId_; }

private:
    bool WriteRunConfigJson(const RunConfig& runConfig, const TokamakEngine& engine);
    bool WriteSummaryCsvRow(const RunConfig& runConfig, const TelemetrySnapshot& telemetry);
    bool WriteRadialProfileRows(const TokamakEngine& engine, const TelemetrySnapshot& telemetry);
    bool WriteMagneticFieldDiagnosticsRows(const TokamakEngine& engine, const TelemetrySnapshot& telemetry);
    bool WriteElectrostaticDiagnosticsRow(const RunConfig& runConfig, const TelemetrySnapshot& telemetry);
    bool WriteSpeedHistogramRows(const TokamakEngine& engine, const TelemetrySnapshot& telemetry);
    bool WritePitchHistogramRows(const TokamakEngine& engine, const TelemetrySnapshot& telemetry);
    bool WriteSolverResidualRow(const TelemetrySnapshot& telemetry);
    bool WriteParticleSnapshotCsv(const TokamakEngine& engine, const TelemetrySnapshot& telemetry);
    bool WriteManifestJson();

    void CloseFiles();
    void SetError(const std::string& message);

    ArtifactExportConfig config_;
    RunConfig runConfig_;

    std::ofstream summaryCsv_;
    std::ofstream radialCsv_;
    std::ofstream magneticFieldCsv_;
    std::ofstream electrostaticDiagnosticsCsv_;
    std::ofstream speedHistogramCsv_;
    std::ofstream pitchHistogramCsv_;
    std::ofstream solverResidualCsv_;

    std::string runDirectory_;
    std::string runId_;
    std::string runConfigRelativePath_;
    std::string summaryRelativePath_;
    std::string radialRelativePath_;
    std::string magneticFieldRelativePath_;
    std::string electrostaticDiagnosticsRelativePath_;
    std::string speedHistogramRelativePath_;
    std::string pitchHistogramRelativePath_;
    std::string solverResidualRelativePath_;
    std::vector<std::string> particleSnapshotRelativePaths_;
    std::string lastError_;

    int lastMetricsStepWritten_ = -1;
    int lastParticleSnapshotStepWritten_ = -1;
    bool started_ = false;
};

}  // namespace tokamak
