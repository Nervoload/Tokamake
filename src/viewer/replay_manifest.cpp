#include "tokamak/viewer/replay_manifest.hpp"

#include <fstream>
#include <regex>
#include <sstream>
#include <string>

namespace tokamak::viewer {
namespace {

bool ReadTextFile(const std::filesystem::path& path, std::string* outText, std::string* errorOut) {
    std::ifstream input(path);
    if (!input.is_open()) {
        if (errorOut != nullptr) {
            *errorOut = "Failed to open file: " + path.string();
        }
        return false;
    }

    std::ostringstream buffer;
    buffer << input.rdbuf();
    if (!input.good() && !input.eof()) {
        if (errorOut != nullptr) {
            *errorOut = "Failed to read file: " + path.string();
        }
        return false;
    }

    *outText = buffer.str();
    return true;
}

bool ExtractStringField(const std::string& text, const char* key, std::string* outValue) {
    const std::string pattern = std::string("\"") + key + "\"\\s*:\\s*\"([^\"]*)\"";
    const std::regex re(pattern);
    std::smatch match;
    if (!std::regex_search(text, match, re)) {
        return false;
    }
    *outValue = match[1].str();
    return true;
}

bool ExtractUInt32Field(const std::string& text, const char* key, uint32_t* outValue) {
    const std::string pattern = std::string("\"") + key + "\"\\s*:\\s*(\\d+)";
    const std::regex re(pattern);
    std::smatch match;
    if (!std::regex_search(text, match, re)) {
        return false;
    }
    const unsigned long value = std::stoul(match[1].str());
    *outValue = static_cast<uint32_t>(value);
    return true;
}

bool ExtractUInt64Field(const std::string& text, const char* key, uint64_t* outValue) {
    const std::string pattern = std::string("\"") + key + "\"\\s*:\\s*(\\d+)";
    const std::regex re(pattern);
    std::smatch match;
    if (!std::regex_search(text, match, re)) {
        return false;
    }
    const unsigned long long value = std::stoull(match[1].str());
    *outValue = static_cast<uint64_t>(value);
    return true;
}

bool ExtractFloatField(const std::string& text, const char* key, float* outValue) {
    const std::string pattern = std::string("\"") + key + "\"\\s*:\\s*([-+]?(?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][-+]?\\d+)?)";
    const std::regex re(pattern);
    std::smatch match;
    if (!std::regex_search(text, match, re)) {
        return false;
    }
    *outValue = std::stof(match[1].str());
    return true;
}

bool ExtractStringArrayField(
    const std::string& text,
    const char* key,
    std::vector<std::string>* outValues) {
    const std::string pattern = std::string("\"") + key + "\"\\s*:\\s*\\[([\\s\\S]*?)\\]";
    const std::regex re(pattern, std::regex::ECMAScript);
    std::smatch match;
    if (!std::regex_search(text, match, re)) {
        return false;
    }

    const std::string arrayBody = match[1].str();
    const std::regex itemRe("\"([^\"]*)\"");
    for (std::sregex_iterator it(arrayBody.begin(), arrayBody.end(), itemRe), end; it != end; ++it) {
        outValues->push_back((*it)[1].str());
    }

    return true;
}

std::string JoinMissingKeys(const std::vector<std::string>& keys) {
    std::ostringstream oss;
    for (std::size_t i = 0; i < keys.size(); ++i) {
        if (i > 0) {
            oss << ", ";
        }
        oss << keys[i];
    }
    return oss.str();
}

}  // namespace

bool ParseManifestV2File(
    const std::filesystem::path& manifestPath,
    ReplayManifest* outManifest,
    std::string* errorOut) {
    if (outManifest == nullptr) {
        if (errorOut != nullptr) {
            *errorOut = "ParseManifestV2File: outManifest is null";
        }
        return false;
    }

    std::string text;
    if (!ReadTextFile(manifestPath, &text, errorOut)) {
        return false;
    }

    ReplayManifest parsed;
    parsed.manifestPath = manifestPath;
    parsed.runDirectory = manifestPath.parent_path();

    std::vector<std::string> missingKeys;
    if (!ExtractUInt32Field(text, "schema_version", &parsed.schemaVersion)) {
        missingKeys.emplace_back("schema_version");
    }
    if (!ExtractStringField(text, "run_id", &parsed.runId)) {
        missingKeys.emplace_back("run_id");
    }
    if (!ExtractStringField(text, "created_utc", &parsed.createdUtc)) {
        missingKeys.emplace_back("created_utc");
    }

    if (!ExtractStringField(text, "run_config", &parsed.files.runConfigJson)) {
        missingKeys.emplace_back("files.run_config");
    }
    if (!ExtractStringField(text, "summary_csv", &parsed.files.summaryCsv)) {
        missingKeys.emplace_back("files.summary_csv");
    }
    if (!ExtractStringField(text, "radial_profiles_csv", &parsed.files.radialProfilesCsv)) {
        missingKeys.emplace_back("files.radial_profiles_csv");
    }
    if (!ExtractStringField(text, "magnetic_field_diagnostics_csv", &parsed.files.magneticFieldDiagnosticsCsv)) {
        missingKeys.emplace_back("files.magnetic_field_diagnostics_csv");
    }
    if (!ExtractStringField(text, "electrostatic_diagnostics_csv", &parsed.files.electrostaticDiagnosticsCsv)) {
        missingKeys.emplace_back("files.electrostatic_diagnostics_csv");
    }
    if (!ExtractStringField(text, "speed_histogram_csv", &parsed.files.speedHistogramCsv)) {
        missingKeys.emplace_back("files.speed_histogram_csv");
    }
    if (!ExtractStringField(text, "pitch_histogram_csv", &parsed.files.pitchHistogramCsv)) {
        missingKeys.emplace_back("files.pitch_histogram_csv");
    }
    if (!ExtractStringField(text, "solver_residual_log_csv", &parsed.files.solverResidualCsv)) {
        missingKeys.emplace_back("files.solver_residual_log_csv");
    }

    if (!ExtractStringArrayField(text, "particle_snapshot_csv_files", &parsed.files.particleSnapshotCsvFiles)) {
        missingKeys.emplace_back("files.particle_snapshot_csv_files");
    }

    if (!missingKeys.empty()) {
        if (errorOut != nullptr) {
            *errorOut = "Manifest missing required key(s): " + JoinMissingKeys(missingKeys);
        }
        return false;
    }

    if (parsed.schemaVersion != 2u) {
        if (errorOut != nullptr) {
            *errorOut = "Unsupported manifest schema_version: " + std::to_string(parsed.schemaVersion) +
                        " (expected 2)";
        }
        return false;
    }

    if (parsed.files.particleSnapshotCsvFiles.empty()) {
        if (errorOut != nullptr) {
            *errorOut = "Manifest has empty files.particle_snapshot_csv_files";
        }
        return false;
    }

    *outManifest = std::move(parsed);
    return true;
}

bool ParseRunConfigV2File(
    const std::filesystem::path& runConfigPath,
    ReplayRunConfig* outRunConfig,
    std::string* errorOut) {
    if (outRunConfig == nullptr) {
        if (errorOut != nullptr) {
            *errorOut = "ParseRunConfigV2File: outRunConfig is null";
        }
        return false;
    }

    std::string text;
    if (!ReadTextFile(runConfigPath, &text, errorOut)) {
        return false;
    }

    ReplayRunConfig parsed;
    ExtractStringField(text, "scenario", &parsed.scenario);

    uint64_t seed = 0;
    if (ExtractUInt64Field(text, "seed", &seed)) {
        parsed.seed = seed;
        parsed.hasSeed = true;
    }

    float majorRadius = 0.0f;
    float minorRadius = 0.0f;
    const bool hasMajor = ExtractFloatField(text, "major_radius_m", &majorRadius);
    const bool hasMinor = ExtractFloatField(text, "minor_radius_m", &minorRadius);

    if (hasMajor && hasMinor && majorRadius > 0.0f && minorRadius > 0.0f) {
        parsed.hasTokamakGeometry = true;
        parsed.majorRadius_m = majorRadius;
        parsed.minorRadius_m = minorRadius;
    }

    *outRunConfig = parsed;
    return true;
}

}  // namespace tokamak::viewer
