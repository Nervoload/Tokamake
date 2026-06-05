#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "tokamak/magnetic_field.hpp"

namespace tokamak::viewer {

struct ReplayManifestFiles {
    std::string runConfigJson;
    std::string summaryCsv;
    std::string radialProfilesCsv;
    std::string magneticFieldDiagnosticsCsv;
    std::string electrostaticDiagnosticsCsv;
    std::string fusionReactivityDiagnosticsCsv;
    std::string wallInteractionBridgeCsv;
    std::string speedHistogramCsv;
    std::string pitchHistogramCsv;
    std::string solverResidualCsv;
    std::string fieldProbeSamplesCsv;
    std::vector<std::string> particleSnapshotCsvFiles;
};

struct ReplayManifest {
    std::filesystem::path manifestPath;
    std::filesystem::path runDirectory;
    uint32_t schemaVersion = 0;
    std::string runId;
    std::string createdUtc;
    ReplayManifestFiles files;
};

struct ReplayRunConfig {
    bool hasTokamakGeometry = false;
    float majorRadius_m = 2.0f;
    float minorRadius_m = 0.5f;
    TokamakConfig tokamakConfig;
    bool hasTokamakConfig = false;
    PlasmaCurrentProfileConfig plasmaCurrentProfile;
    bool hasPlasmaCurrentProfile = false;
    ElectricFieldMode electricFieldMode = ElectricFieldMode::Placeholder;
    bool hasElectricFieldMode = false;
    NBIConfig nbiConfig;
    bool hasNbiConfig = false;
    uint64_t maxParticlesPerSnapshot = 0;
    bool hasMaxParticlesPerSnapshot = false;
    double startupRampDuration_s = 0.0;
    double fusionStartDelay_s = 0.0;
    double fusionRampDuration_s = 0.0;
    std::string scenario;
    uint64_t seed = 0;
    bool hasSeed = false;
};

bool ParseManifestV2File(
    const std::filesystem::path& manifestPath,
    ReplayManifest* outManifest,
    std::string* errorOut);

bool ParseRunConfigV2File(
    const std::filesystem::path& runConfigPath,
    ReplayRunConfig* outRunConfig,
    std::string* errorOut);

}  // namespace tokamak::viewer
