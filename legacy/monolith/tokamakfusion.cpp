#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include "snapshot_stream.h"
#include "tokamak_engine.h"

namespace {

void PrintUsage() {
    std::cout
        << "Usage: tokamakfusion [options]\n"
        << "  --scenario <cold|ignition|failure>\n"
        << "  --steps <N>\n"
        << "  --dt <seconds>\n"
        << "  --record <file.tksnap>\n";
}

}  // namespace

int main(int argc, char** argv) {
    Scenario scenario = Scenario::NBI_IGNITION;
    int totalSteps = 10000;
    float dt = 1.0e-7f;
    std::string recordPath;

    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        auto nextValue = [&](const char* flag) -> const char* {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << flag << "\n";
                return nullptr;
            }
            return argv[++i];
        };

        if (arg == "--scenario") {
            const char* value = nextValue("--scenario");
            if (value == nullptr) {
                return 1;
            }
            if (std::strcmp(value, "cold") == 0) {
                scenario = Scenario::COLD_VACUUM;
            } else if (std::strcmp(value, "ignition") == 0) {
                scenario = Scenario::NBI_IGNITION;
            } else if (std::strcmp(value, "failure") == 0) {
                scenario = Scenario::MAGNETIC_FAILURE;
            } else {
                std::cerr << "Invalid scenario\n";
                return 1;
            }
        } else if (arg == "--steps") {
            const char* value = nextValue("--steps");
            if (value == nullptr) {
                return 1;
            }
            totalSteps = std::max(1, std::atoi(value));
        } else if (arg == "--dt") {
            const char* value = nextValue("--dt");
            if (value == nullptr) {
                return 1;
            }
            dt = std::max(1.0e-10f, static_cast<float>(std::atof(value)));
        } else if (arg == "--record") {
            const char* value = nextValue("--record");
            if (value == nullptr) {
                return 1;
            }
            recordPath = value;
        } else if (arg == "--help" || arg == "-h") {
            PrintUsage();
            return 0;
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            PrintUsage();
            return 1;
        }
    }

    TokamakEngine engine(scenario);

    SnapshotStreamWriter writer;
    bool shouldRecord = !recordPath.empty();
    if (shouldRecord && !writer.Open(recordPath)) {
        std::cerr << "Warning: snapshot recording disabled: " << writer.LastError() << "\n";
        shouldRecord = false;
    }

    SnapshotRequest request;
    request.maxSampledParticles = 4096;
    request.includeVelocities = true;
    request.includeMassAndCharge = true;

    SimulationSnapshot snapshot;
    snapshot.Reserve(request.maxSampledParticles);

    std::cout << "Starting main physics loop...\n";
    std::cout << "------------------------------------------------------------\n";

    for (int step = 0; step <= totalSteps; ++step) {
        engine.Step(dt);

        if (step % 100 == 0) {
            engine.PrintTelemetry(step);
        }

        if (shouldRecord && (step % 10 == 0)) {
            engine.ExportSnapshot(snapshot, request);
            if (!writer.Write(snapshot)) {
                std::cerr << "Warning: snapshot write failed, disabling recorder: "
                          << writer.LastError() << "\n";
                shouldRecord = false;
                writer.Close();
            }
        }
    }

    std::cout << "------------------------------------------------------------\n";
    std::cout << "Simulation Complete.\n";

    writer.Close();
    return 0;
}
