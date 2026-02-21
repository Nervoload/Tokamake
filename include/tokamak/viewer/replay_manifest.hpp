#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace tokamak::viewer {

struct ReplayManifestFiles {
    std::string runConfigJson;
    std::string summaryCsv;
    std::string radialProfilesCsv;
    std::string magneticFieldDiagnosticsCsv;
    std::string electrostaticDiagnosticsCsv;
    std::string speedHistogramCsv;
    std::string pitchHistogramCsv;
    std::string solverResidualCsv;
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
