#include "tokamak/artifact_export.hpp"
#include "tokamak/engine.hpp"

#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using tokamak::ArtifactExportConfig;
using tokamak::ChargeAssignmentScheme;
using tokamak::CurrentProfilePoint;
using tokamak::ElectricFieldMode;
using tokamak::ElectrostaticBoundaryCondition;
using tokamak::PlasmaCurrentProfileKind;
using tokamak::RunConfig;
using tokamak::Scenario;

struct CliOptions {
    RunConfig runConfig;
    ArtifactExportConfig artifactConfig;
    bool exportArtifacts = true;
    bool helpRequested = false;
};

bool ParseUInt32(const char* text, uint32_t* outValue) {
    char* end = nullptr;
    errno = 0;
    const unsigned long parsed = std::strtoul(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || parsed > std::numeric_limits<uint32_t>::max()) {
        return false;
    }
    *outValue = static_cast<uint32_t>(parsed);
    return true;
}

bool ParseSizeT(const char* text, std::size_t* outValue) {
    char* end = nullptr;
    errno = 0;
    const unsigned long long parsed = std::strtoull(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0') {
        return false;
    }
    *outValue = static_cast<std::size_t>(parsed);
    return true;
}

bool ParseInt(const char* text, int* outValue) {
    char* end = nullptr;
    errno = 0;
    const long parsed = std::strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || parsed < std::numeric_limits<int>::min() ||
        parsed > std::numeric_limits<int>::max()) {
        return false;
    }
    *outValue = static_cast<int>(parsed);
    return true;
}

bool ParseFloat(const char* text, float* outValue) {
    char* end = nullptr;
    errno = 0;
    const float parsed = std::strtof(text, &end);
    if (errno != 0 || end == text || *end != '\0') {
        return false;
    }
    *outValue = parsed;
    return true;
}

bool ParseDouble(const char* text, double* outValue) {
    char* end = nullptr;
    errno = 0;
    const double parsed = std::strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0') {
        return false;
    }
    *outValue = parsed;
    return true;
}

std::string TrimAscii(const std::string& text) {
    const std::size_t first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return std::string();
    }
    const std::size_t last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, (last - first) + 1);
}

bool ParseCurrentProfileTableFile(
    const std::string& path,
    std::vector<CurrentProfilePoint>* outPoints,
    std::string* errorOut) {
    std::ifstream input(path);
    if (!input.is_open()) {
        if (errorOut != nullptr) {
            *errorOut = "Unable to open current-profile table file: " + path;
        }
        return false;
    }

    std::vector<CurrentProfilePoint> points;
    std::string line;
    int lineNumber = 0;
    bool parsedAnyData = false;

    while (std::getline(input, line)) {
        ++lineNumber;
        std::string trimmed = TrimAscii(line);
        if (trimmed.empty() || trimmed[0] == '#') {
            continue;
        }

        for (char& c : trimmed) {
            if (c == ',' || c == ';' || c == '\t') {
                c = ' ';
            }
        }

        std::istringstream lineStream(trimmed);
        double normalizedRadius = 0.0;
        double enclosedFraction = 0.0;
        if (!(lineStream >> normalizedRadius >> enclosedFraction)) {
            if (!parsedAnyData) {
                continue;  // Allow a leading header row.
            }
            if (errorOut != nullptr) {
                *errorOut = "Invalid current-profile table row at line " + std::to_string(lineNumber);
            }
            return false;
        }

        std::string trailingToken;
        if (lineStream >> trailingToken) {
            if (errorOut != nullptr) {
                *errorOut = "Too many values in current-profile table row at line " + std::to_string(lineNumber);
            }
            return false;
        }

        if (!std::isfinite(normalizedRadius) || !std::isfinite(enclosedFraction) ||
            normalizedRadius < 0.0 || normalizedRadius > 1.0 ||
            enclosedFraction < 0.0 || enclosedFraction > 1.0) {
            if (errorOut != nullptr) {
                *errorOut =
                    "Current-profile table values must be finite and in [0,1] at line " +
                    std::to_string(lineNumber);
            }
            return false;
        }

        points.push_back(CurrentProfilePoint{
            static_cast<float>(normalizedRadius),
            static_cast<float>(enclosedFraction),
        });
        parsedAnyData = true;
    }

    if (points.empty()) {
        if (errorOut != nullptr) {
            *errorOut = "Current-profile table file has no numeric data rows: " + path;
        }
        return false;
    }

    *outPoints = std::move(points);
    return true;
}

void PrintUsage(const char* argv0) {
    std::cout << "Usage: " << argv0 << " [options]\n"
              << "  --scenario <cold|ignition|failure>\n"
              << "  --seed <uint32>\n"
              << "  --dt <seconds>\n"
              << "  --steps <int>\n"
              << "  --telemetry-every <int>\n"
              << "  --particle-cap <uint>\n"
              << "  --current-profile <uniform|parabolic|custom>\n"
              << "  --current-profile-axis-epsilon <meters>\n"
              << "  --current-profile-axis-blend <meters>\n"
              << "  --current-profile-table <csv_path>\n"
              << "  --mag-field-bins <N>\n"
              << "  --mag-field-dt-safety <fraction>\n"
              << "  --electric-field-mode <placeholder|electrostatic>\n"
              << "  --electrostatic-bc <dirichlet0|neumann0>\n"
              << "  --charge-assignment <ngp|cic>\n"
              << "  --electrostatic-grid-bins <N>\n"
              << "  --electrostatic-tol <value>\n"
              << "  --electrostatic-max-iters <N>\n"
              << "  --electrostatic-omega <value>\n"
              << "  --artifacts-root <path>\n"
              << "  --artifact-every <int>\n"
              << "  --particle-snapshot-every <int>\n"
              << "  --particle-snapshot-max <uint>\n"
              << "  --no-artifacts\n"
              << "  --help\n";
}

bool ParseArgs(int argc, char** argv, CliOptions* options) {
    for (int i = 1; i < argc; ++i) {
        const std::string_view arg(argv[i]);
        auto needValue = [&](const char* optionName) -> const char* {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << optionName << "\n";
                return nullptr;
            }
            return argv[++i];
        };

        if (arg == "--help") {
            PrintUsage(argv[0]);
            options->helpRequested = true;
            return false;
        }
        if (arg == "--scenario") {
            const char* value = needValue("--scenario");
            if (value == nullptr) {
                return false;
            }
            Scenario scenario;
            if (!tokamak::ParseScenario(value, &scenario)) {
                std::cerr << "Invalid scenario: " << value << "\n";
                return false;
            }
            options->runConfig.scenario = scenario;
            continue;
        }
        if (arg == "--seed") {
            const char* value = needValue("--seed");
            if (value == nullptr) {
                return false;
            }
            uint32_t seed = 0;
            if (!ParseUInt32(value, &seed)) {
                std::cerr << "Invalid seed: " << value << "\n";
                return false;
            }
            options->runConfig.seed = seed;
            continue;
        }
        if (arg == "--dt") {
            const char* value = needValue("--dt");
            if (value == nullptr) {
                return false;
            }
            float dt_s = 0.0f;
            if (!ParseFloat(value, &dt_s) || dt_s < 0.0f) {
                std::cerr << "Invalid dt: " << value << "\n";
                return false;
            }
            options->runConfig.timeStep_s = dt_s;
            continue;
        }
        if (arg == "--steps") {
            const char* value = needValue("--steps");
            if (value == nullptr) {
                return false;
            }
            int steps = 0;
            if (!ParseInt(value, &steps) || steps < 0) {
                std::cerr << "Invalid steps: " << value << "\n";
                return false;
            }
            options->runConfig.totalSteps = steps;
            continue;
        }
        if (arg == "--telemetry-every") {
            const char* value = needValue("--telemetry-every");
            if (value == nullptr) {
                return false;
            }
            int cadence = 0;
            if (!ParseInt(value, &cadence) || cadence <= 0) {
                std::cerr << "Invalid telemetry cadence: " << value << "\n";
                return false;
            }
            options->runConfig.telemetryEveryNSteps = cadence;
            continue;
        }
        if (arg == "--particle-cap") {
            const char* value = needValue("--particle-cap");
            if (value == nullptr) {
                return false;
            }
            std::size_t cap = 0;
            if (!ParseSizeT(value, &cap) || cap == 0) {
                std::cerr << "Invalid particle cap: " << value << "\n";
                return false;
            }
            options->runConfig.particleCap = cap;
            continue;
        }
        if (arg == "--current-profile") {
            const char* value = needValue("--current-profile");
            if (value == nullptr) {
                return false;
            }
            PlasmaCurrentProfileKind kind = PlasmaCurrentProfileKind::Uniform;
            if (!tokamak::ParsePlasmaCurrentProfileKind(value, &kind)) {
                std::cerr << "Invalid current-profile value: " << value << "\n";
                return false;
            }
            options->runConfig.plasmaCurrentProfileKind = kind;
            continue;
        }
        if (arg == "--current-profile-axis-epsilon") {
            const char* value = needValue("--current-profile-axis-epsilon");
            if (value == nullptr) {
                return false;
            }
            float epsilon_m = 0.0f;
            if (!ParseFloat(value, &epsilon_m) || epsilon_m <= 0.0f) {
                std::cerr << "Invalid current-profile-axis-epsilon: " << value << "\n";
                return false;
            }
            options->runConfig.plasmaCurrentAxisEpsilon_m = epsilon_m;
            continue;
        }
        if (arg == "--current-profile-axis-blend") {
            const char* value = needValue("--current-profile-axis-blend");
            if (value == nullptr) {
                return false;
            }
            float blendRadius_m = 0.0f;
            if (!ParseFloat(value, &blendRadius_m) || blendRadius_m < 0.0f) {
                std::cerr << "Invalid current-profile-axis-blend: " << value << "\n";
                return false;
            }
            options->runConfig.plasmaCurrentCustomAxisBlendRadius_m = blendRadius_m;
            continue;
        }
        if (arg == "--current-profile-table") {
            const char* value = needValue("--current-profile-table");
            if (value == nullptr) {
                return false;
            }
            std::vector<CurrentProfilePoint> parsedPoints;
            std::string tableError;
            if (!ParseCurrentProfileTableFile(value, &parsedPoints, &tableError)) {
                std::cerr << tableError << "\n";
                return false;
            }
            options->runConfig.plasmaCurrentCustomTable = std::move(parsedPoints);
            continue;
        }
        if (arg == "--mag-field-bins") {
            const char* value = needValue("--mag-field-bins");
            if (value == nullptr) {
                return false;
            }
            std::size_t bins = 0;
            if (!ParseSizeT(value, &bins) || bins == 0) {
                std::cerr << "Invalid mag-field-bins: " << value << "\n";
                return false;
            }
            options->runConfig.magneticFieldRadialBinCount = bins;
            continue;
        }
        if (arg == "--mag-field-dt-safety") {
            const char* value = needValue("--mag-field-dt-safety");
            if (value == nullptr) {
                return false;
            }
            double safetyFraction = 0.0;
            if (!ParseDouble(value, &safetyFraction) || safetyFraction < 0.0) {
                std::cerr << "Invalid mag-field-dt-safety: " << value << "\n";
                return false;
            }
            options->runConfig.magneticFieldDtSafetyFraction = safetyFraction;
            continue;
        }
        if (arg == "--electric-field-mode") {
            const char* value = needValue("--electric-field-mode");
            if (value == nullptr) {
                return false;
            }
            ElectricFieldMode mode = ElectricFieldMode::Placeholder;
            if (!tokamak::ParseElectricFieldMode(value, &mode)) {
                std::cerr << "Invalid electric-field-mode value: " << value << "\n";
                return false;
            }
            options->runConfig.electricFieldMode = mode;
            continue;
        }
        if (arg == "--electrostatic-bc") {
            const char* value = needValue("--electrostatic-bc");
            if (value == nullptr) {
                return false;
            }
            ElectrostaticBoundaryCondition boundaryCondition = ElectrostaticBoundaryCondition::DirichletZero;
            if (!tokamak::ParseElectrostaticBoundaryCondition(value, &boundaryCondition)) {
                std::cerr << "Invalid electrostatic-bc value: " << value << "\n";
                return false;
            }
            options->runConfig.electrostaticBoundaryCondition = boundaryCondition;
            continue;
        }
        if (arg == "--charge-assignment") {
            const char* value = needValue("--charge-assignment");
            if (value == nullptr) {
                return false;
            }
            ChargeAssignmentScheme scheme = ChargeAssignmentScheme::CIC;
            if (!tokamak::ParseChargeAssignmentScheme(value, &scheme)) {
                std::cerr << "Invalid charge-assignment value: " << value << "\n";
                return false;
            }
            options->runConfig.chargeAssignmentScheme = scheme;
            continue;
        }
        if (arg == "--electrostatic-grid-bins") {
            const char* value = needValue("--electrostatic-grid-bins");
            if (value == nullptr) {
                return false;
            }
            std::size_t bins = 0;
            if (!ParseSizeT(value, &bins) || bins < 2) {
                std::cerr << "Invalid electrostatic-grid-bins: " << value << "\n";
                return false;
            }
            options->runConfig.electrostaticGridBinCount = bins;
            continue;
        }
        if (arg == "--electrostatic-tol") {
            const char* value = needValue("--electrostatic-tol");
            if (value == nullptr) {
                return false;
            }
            double tolerance = 0.0;
            if (!ParseDouble(value, &tolerance) || tolerance <= 0.0) {
                std::cerr << "Invalid electrostatic-tol: " << value << "\n";
                return false;
            }
            options->runConfig.electrostaticSolverTolerance = tolerance;
            continue;
        }
        if (arg == "--electrostatic-max-iters") {
            const char* value = needValue("--electrostatic-max-iters");
            if (value == nullptr) {
                return false;
            }
            uint32_t maxIterations = 0;
            if (!ParseUInt32(value, &maxIterations) || maxIterations == 0) {
                std::cerr << "Invalid electrostatic-max-iters: " << value << "\n";
                return false;
            }
            options->runConfig.electrostaticSolverMaxIterations = maxIterations;
            continue;
        }
        if (arg == "--electrostatic-omega") {
            const char* value = needValue("--electrostatic-omega");
            if (value == nullptr) {
                return false;
            }
            double omega = 0.0;
            if (!ParseDouble(value, &omega) || omega <= 0.0 || omega >= 2.0) {
                std::cerr << "Invalid electrostatic-omega: " << value << " (expected 0 < omega < 2)\n";
                return false;
            }
            options->runConfig.electrostaticSorOmega = omega;
            continue;
        }
        if (arg == "--artifacts-root") {
            const char* value = needValue("--artifacts-root");
            if (value == nullptr) {
                return false;
            }
            options->artifactConfig.outputRootDirectory = value;
            continue;
        }
        if (arg == "--artifact-every") {
            const char* value = needValue("--artifact-every");
            if (value == nullptr) {
                return false;
            }
            int cadence = 0;
            if (!ParseInt(value, &cadence) || cadence <= 0) {
                std::cerr << "Invalid artifact cadence: " << value << "\n";
                return false;
            }
            options->artifactConfig.metricsEveryNSteps = cadence;
            continue;
        }
        if (arg == "--particle-snapshot-every") {
            const char* value = needValue("--particle-snapshot-every");
            if (value == nullptr) {
                return false;
            }
            int cadence = 0;
            if (!ParseInt(value, &cadence) || cadence <= 0) {
                std::cerr << "Invalid particle snapshot cadence: " << value << "\n";
                return false;
            }
            options->artifactConfig.particleSnapshotEveryNSteps = cadence;
            continue;
        }
        if (arg == "--particle-snapshot-max") {
            const char* value = needValue("--particle-snapshot-max");
            if (value == nullptr) {
                return false;
            }
            std::size_t count = 0;
            if (!ParseSizeT(value, &count) || count == 0) {
                std::cerr << "Invalid particle snapshot max: " << value << "\n";
                return false;
            }
            options->artifactConfig.maxParticlesPerSnapshot = count;
            continue;
        }
        if (arg == "--no-artifacts") {
            options->exportArtifacts = false;
            continue;
        }

        std::cerr << "Unknown option: " << arg << "\n";
        return false;
    }

    if (options->runConfig.plasmaCurrentProfileKind == PlasmaCurrentProfileKind::CustomTable &&
        options->runConfig.plasmaCurrentCustomTable.empty()) {
        std::cerr << "Custom current profile requires --current-profile-table <csv_path>\n";
        return false;
    }
    if (options->runConfig.plasmaCurrentProfileKind != PlasmaCurrentProfileKind::CustomTable &&
        !options->runConfig.plasmaCurrentCustomTable.empty()) {
        std::cerr << "--current-profile-table requires --current-profile custom\n";
        return false;
    }

    return true;
}

const char* ScenarioDescription(Scenario scenario) {
    switch (scenario) {
        case Scenario::ColdVacuum:
            return "Cold Vacuum. No NBI Heating.";
        case Scenario::NbiIgnition:
            return "NBI Ignition. High power heating to target fusion.";
        case Scenario::MagneticFailure:
            return "Magnetic Failure. Poloidal field collapsing.";
    }
    return "Unknown";
}

}  // namespace

int main(int argc, char** argv) {
    CliOptions options;
    if (!ParseArgs(argc, argv, &options)) {
        return options.helpRequested ? 0 : 1;
    }

    RunConfig& runConfig = options.runConfig;
    tokamak::TokamakEngine engine(runConfig);

    std::cout << "--- INITIALIZING TOKAMAK SIMULATION ---\n";
    std::cout << "SCENARIO: " << ScenarioDescription(runConfig.scenario) << "\n";
    std::cout << "RUN MANIFEST | scenario=" << tokamak::ScenarioName(runConfig.scenario)
              << " seed=" << engine.ActiveSeed() << " dt_s=" << runConfig.timeStep_s
              << " steps=" << runConfig.totalSteps << " telemetry_every=" << runConfig.telemetryEveryNSteps
              << " particle_cap=" << runConfig.particleCap
              << " current_profile=" << tokamak::PlasmaCurrentProfileKindName(runConfig.plasmaCurrentProfileKind)
              << " mag_field_bins=" << runConfig.magneticFieldRadialBinCount
              << " mag_field_dt_safety=" << runConfig.magneticFieldDtSafetyFraction
              << " electric_field_mode=" << tokamak::ElectricFieldModeName(runConfig.electricFieldMode)
              << " electrostatic_bc=" << tokamak::ElectrostaticBoundaryConditionName(runConfig.electrostaticBoundaryCondition)
              << " charge_assignment=" << tokamak::ChargeAssignmentSchemeName(runConfig.chargeAssignmentScheme)
              << " electrostatic_grid_bins=" << runConfig.electrostaticGridBinCount
              << " electrostatic_tol=" << runConfig.electrostaticSolverTolerance
              << " electrostatic_max_iters=" << runConfig.electrostaticSolverMaxIterations
              << " electrostatic_omega=" << runConfig.electrostaticSorOmega << "\n";

    tokamak::Milestone7ArtifactExporter artifactExporter;
    if (options.exportArtifacts) {
        if (!artifactExporter.Start(runConfig, engine, options.artifactConfig)) {
            std::cerr << "Artifact export start failed: " << artifactExporter.LastError() << "\n";
            return 1;
        }

        const tokamak::TelemetrySnapshot initialSnapshot = engine.Snapshot(0);
        if (!artifactExporter.WriteStep(engine, initialSnapshot, true)) {
            std::cerr << "Artifact export write failed at step 0: " << artifactExporter.LastError() << "\n";
            return 1;
        }

        std::cout << "ARTIFACT RUN | id=" << artifactExporter.RunId()
                  << " dir=" << artifactExporter.RunDirectory() << "\n";
    }

    std::cout << "Starting main physics loop...\n";
    std::cout << "------------------------------------------------------------\n";

    engine.PrintTelemetry(0);
    for (int step = 1; step <= runConfig.totalSteps; ++step) {
        engine.Step(runConfig.timeStep_s);
        if (step % runConfig.telemetryEveryNSteps == 0) {
            engine.PrintTelemetry(step);
        }
        if (options.exportArtifacts) {
            const tokamak::TelemetrySnapshot telemetry = engine.Snapshot(step);
            if (!artifactExporter.WriteStep(engine, telemetry, false)) {
                std::cerr << "Artifact export write failed at step " << step << ": "
                          << artifactExporter.LastError() << "\n";
                return 1;
            }
        }
    }

    if (options.exportArtifacts) {
        const tokamak::TelemetrySnapshot finalTelemetry = engine.Snapshot(runConfig.totalSteps);
        if (!artifactExporter.WriteStep(engine, finalTelemetry, true)) {
            std::cerr << "Artifact export final write failed: " << artifactExporter.LastError() << "\n";
            return 1;
        }
        if (!artifactExporter.Finish()) {
            std::cerr << "Artifact export finalize failed: " << artifactExporter.LastError() << "\n";
            return 1;
        }
        std::cout << "ARTIFACTS WRITTEN | " << artifactExporter.RunDirectory() << "\n";
    }

    std::cout << "------------------------------------------------------------\n";
    std::cout << "Simulation Complete.\n";
    return 0;
}
