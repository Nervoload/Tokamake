#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "tokamak/collision.hpp"
#include "tokamak/engine.hpp"

namespace {

std::string RunCommandCapture(const std::string& command) {
    std::array<char, 256> buffer{};
    std::string output;

    FILE* pipe = popen(command.c_str(), "r");
    if (pipe == nullptr) {
        throw std::runtime_error("popen failed");
    }

    while (fgets(buffer.data(), static_cast<int>(buffer.size()), pipe) != nullptr) {
        output += buffer.data();
    }

    const int rc = pclose(pipe);
    if (rc != 0) {
        std::ostringstream oss;
        oss << "command failed with rc=" << rc << ": " << command << "\n" << output;
        throw std::runtime_error(oss.str());
    }

    return output;
}

std::size_t CountOccurrences(const std::string& text, const std::string& needle) {
    std::size_t count = 0;
    std::size_t pos = 0;
    while ((pos = text.find(needle, pos)) != std::string::npos) {
        ++count;
        pos += needle.size();
    }
    return count;
}

tokamak::ParticleSystem BuildCollisionPool(bool groupedOrder) {
    constexpr int kPairsPerSpecies = 1000;

    tokamak::ParticleSystem particles(2200);
    const tokamak::Vec3 position(2.0f, 0.0f, 0.0f);
    const tokamak::Vec3 velocityD(2.0e6f, 0.0f, 0.0f);
    const tokamak::Vec3 velocityT(-2.0e6f, 0.0f, 0.0f);

    auto addD = [&]() {
        return particles.AddParticle(
            position,
            velocityD,
            tokamak::constants::kMassDeuterium_kg,
            tokamak::constants::kElementaryCharge_C,
            tokamak::ParticleType::Deuterium);
    };
    auto addT = [&]() {
        return particles.AddParticle(
            position,
            velocityT,
            tokamak::constants::kMassTritium_kg,
            tokamak::constants::kElementaryCharge_C,
            tokamak::ParticleType::Tritium);
    };

    if (groupedOrder) {
        for (int i = 0; i < kPairsPerSpecies; ++i) {
            EXPECT_TRUE(addD());
        }
        for (int i = 0; i < kPairsPerSpecies; ++i) {
            EXPECT_TRUE(addT());
        }
    } else {
        for (int i = 0; i < kPairsPerSpecies; ++i) {
            EXPECT_TRUE(addD());
            EXPECT_TRUE(addT());
        }
    }

    EXPECT_EQ(particles.Size(), static_cast<std::size_t>(kPairsPerSpecies * 2));
    return particles;
}

double ComputeMedian(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    if (values.empty()) {
        return 0.0;
    }
    const std::size_t mid = values.size() / 2;
    if (values.size() % 2 == 0) {
        return 0.5 * (values[mid - 1] + values[mid]);
    }
    return values[mid];
}

std::string ExtractArtifactDirectory(const std::string& output) {
    const std::string key = "ARTIFACT RUN | id=";
    const std::size_t keyPos = output.find(key);
    if (keyPos == std::string::npos) {
        return std::string();
    }

    const std::string dirKey = " dir=";
    const std::size_t dirPos = output.find(dirKey, keyPos);
    if (dirPos == std::string::npos) {
        return std::string();
    }
    const std::size_t start = dirPos + dirKey.size();
    const std::size_t end = output.find('\n', start);
    return output.substr(start, end == std::string::npos ? std::string::npos : (end - start));
}

std::string ReadFirstLine(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input.is_open()) {
        throw std::runtime_error("failed to open file: " + path.string());
    }
    std::string line;
    std::getline(input, line);
    return line;
}

std::string ReadWholeFile(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input.is_open()) {
        throw std::runtime_error("failed to open file: " + path.string());
    }
    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

std::size_t CountDataLines(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input.is_open()) {
        throw std::runtime_error("failed to open file: " + path.string());
    }
    std::size_t lineCount = 0;
    std::string line;
    while (std::getline(input, line)) {
        ++lineCount;
    }
    return (lineCount == 0) ? 0 : (lineCount - 1);
}

}  // namespace

TEST(CollisionTest, FusionNoCapacityNoFuelConsumption) {
    tokamak::ParticleSystem particles(2);
    EXPECT_TRUE(particles.AddParticle(
        tokamak::Vec3(2.0f, 0.0f, 0.0f),
        tokamak::Vec3(1.0e6f, 0.0f, 0.0f),
        tokamak::constants::kMassDeuterium_kg,
        tokamak::constants::kElementaryCharge_C,
        tokamak::ParticleType::Deuterium));
    EXPECT_TRUE(particles.AddParticle(
        tokamak::Vec3(2.0f, 0.0f, 0.0f),
        tokamak::Vec3(-1.0e6f, 0.0f, 0.0f),
        tokamak::constants::kMassTritium_kg,
        tokamak::constants::kElementaryCharge_C,
        tokamak::ParticleType::Tritium));

    tokamak::RuntimeCounters counters;
    tokamak::EnergyChargeBudget budget;

    std::vector<tokamak::PendingFusionEvent> events;
    events.push_back(tokamak::PendingFusionEvent{
        0,
        1,
        tokamak::Vec3(2.0f, 0.0f, 0.0f),
        tokamak::Vec3(0.0f, 0.0f, 0.0f),
    });

    tokamak::ApplyCollisionEvents(particles, events, counters, budget);

    EXPECT_EQ(particles.Size(), static_cast<std::size_t>(2));
    EXPECT_NE(particles.SpeciesAt(0), tokamak::ParticleType::Dead);
    EXPECT_NE(particles.SpeciesAt(1), tokamak::ParticleType::Dead);
    EXPECT_EQ(counters.rejectedFusionAsh, static_cast<uint64_t>(1));
    EXPECT_EQ(counters.fusionAccepted, static_cast<uint64_t>(0));
}

TEST(CollisionTest, CollisionOrderBiasReduced) {
    constexpr int kSeeds = 30;
    const float dt_s = 1.0e-9f;

    std::vector<double> relativeDeltas;
    relativeDeltas.reserve(kSeeds);

    for (int seed = 1; seed <= kSeeds; ++seed) {
        tokamak::ParticleSystem grouped = BuildCollisionPool(true);
        tokamak::ParticleSystem interleaved = BuildCollisionPool(false);

        tokamak::SpatialGrid groupedGrid(2.5f, 0.2f);
        tokamak::SpatialGrid interleavedGrid(2.5f, 0.2f);
        tokamak::SortParticlesIntoGrid(grouped, groupedGrid, 2.5f);
        tokamak::SortParticlesIntoGrid(interleaved, interleavedGrid, 2.5f);

        std::mt19937 rngA(static_cast<uint32_t>(seed));
        std::mt19937 rngB(static_cast<uint32_t>(seed));

        std::vector<tokamak::PendingFusionEvent> eventsA;
        std::vector<tokamak::PendingFusionEvent> eventsB;

        tokamak::RuntimeCounters countersA;
        tokamak::RuntimeCounters countersB;
        tokamak::EnergyChargeBudget budgetA;
        tokamak::EnergyChargeBudget budgetB;

        tokamak::SelectCollisionEvents(grouped, groupedGrid, rngA, dt_s, eventsA);
        tokamak::SelectCollisionEvents(interleaved, interleavedGrid, rngB, dt_s, eventsB);
        tokamak::ApplyCollisionEvents(grouped, eventsA, countersA, budgetA);
        tokamak::ApplyCollisionEvents(interleaved, eventsB, countersB, budgetB);

        const double acceptedA = static_cast<double>(countersA.fusionAccepted);
        const double acceptedB = static_cast<double>(countersB.fusionAccepted);
        const double denom = std::max(1.0, std::max(acceptedA, acceptedB));
        relativeDeltas.push_back(std::fabs(acceptedA - acceptedB) / denom);
    }

    const double medianDelta = ComputeMedian(relativeDeltas);
    const double maxDelta = *std::max_element(relativeDeltas.begin(), relativeDeltas.end());

    EXPECT_LE(medianDelta, 0.05);
    EXPECT_LE(maxDelta, 0.15);
}

TEST(IntegrationTest, TelemetryContainsRequiredCountersAndSmokeOutput) {
    const std::string command = std::string("\"") + TOKAMAKFUSION_PATH +
                                "\" --scenario cold --seed 21 --steps 100 --telemetry-every 100";
    const std::string output = RunCommandCapture(command);

    EXPECT_GT(output.size(), static_cast<std::size_t>(0));
    EXPECT_NE(output.find("RUN MANIFEST |"), std::string::npos);
    EXPECT_NE(output.find("seed="), std::string::npos);
    EXPECT_NE(output.find("cap-hit:"), std::string::npos);
    EXPECT_NE(output.find("rejected-injection:"), std::string::npos);
    EXPECT_NE(output.find("rejected-ash:"), std::string::npos);
    EXPECT_NE(output.find("out-of-domain-clamp:"), std::string::npos);
    EXPECT_NE(output.find("max-cell-reaction:"), std::string::npos);
    EXPECT_NE(output.find("Bmax_T:"), std::string::npos);
    EXPECT_NE(output.find("dt_gyro_recommended_s:"), std::string::npos);
    EXPECT_NE(output.find("Emax_V_per_m:"), std::string::npos);
    EXPECT_NE(output.find("solver_iters:"), std::string::npos);
    EXPECT_NE(output.find("residual_status:"), std::string::npos);
    EXPECT_NE(output.find("Simulation Complete."), std::string::npos);
    EXPECT_GE(CountOccurrences(output, "[Step"), static_cast<std::size_t>(2));
}

TEST(IntegrationTest, CliSupportsCurrentProfileModesAndCustomTablePath) {
    const std::string uniformCommand = std::string("\"") + TOKAMAKFUSION_PATH +
                                       "\" --scenario cold --seed 301 --steps 20 --telemetry-every 20"
                                       " --no-artifacts --current-profile uniform --mag-field-bins 8"
                                       " --mag-field-dt-safety 0.04";
    const std::string uniformOutput = RunCommandCapture(uniformCommand);
    EXPECT_NE(uniformOutput.find("current_profile=uniform"), std::string::npos);
    EXPECT_NE(uniformOutput.find("Simulation Complete."), std::string::npos);

    const std::string parabolicCommand = std::string("\"") + TOKAMAKFUSION_PATH +
                                         "\" --scenario cold --seed 302 --steps 20 --telemetry-every 20"
                                         " --no-artifacts --current-profile parabolic --mag-field-bins 8";
    const std::string parabolicOutput = RunCommandCapture(parabolicCommand);
    EXPECT_NE(parabolicOutput.find("current_profile=parabolic"), std::string::npos);
    EXPECT_NE(parabolicOutput.find("Simulation Complete."), std::string::npos);

    const auto nonce = std::chrono::steady_clock::now().time_since_epoch().count();
    const std::filesystem::path customTablePath =
        std::filesystem::temp_directory_path() / ("tokamak_profile_table_" + std::to_string(nonce) + ".csv");
    {
        std::ofstream out(customTablePath);
        ASSERT_TRUE(out.is_open());
        out << "rho,enc\n";
        out << "0.9,0.4\n";
        out << "0.1,0.2\n";
        out << "0.5,0.8\n";
    }

    const std::string customCommand = std::string("\"") + TOKAMAKFUSION_PATH +
                                      "\" --scenario cold --seed 303 --steps 20 --telemetry-every 20"
                                      " --no-artifacts --current-profile custom --current-profile-table \"" +
                                      customTablePath.string() + "\"";
    const std::string customOutput = RunCommandCapture(customCommand);
    EXPECT_NE(customOutput.find("current_profile=custom"), std::string::npos);
    EXPECT_NE(customOutput.find("Simulation Complete."), std::string::npos);

    std::filesystem::remove(customTablePath);
}

TEST(IntegrationTest, CliSupportsElectrostaticModeAndFlags) {
    const std::string command = std::string("\"") + TOKAMAKFUSION_PATH +
                                "\" --scenario cold --seed 613 --steps 4 --telemetry-every 2 --no-artifacts"
                                " --electric-field-mode electrostatic --electrostatic-bc dirichlet0"
                                " --charge-assignment cic --electrostatic-grid-bins 8 --electrostatic-tol 1e-5"
                                " --electrostatic-max-iters 100 --electrostatic-omega 1.5";
    const std::string output = RunCommandCapture(command);

    EXPECT_NE(output.find("electric_field_mode=electrostatic"), std::string::npos);
    EXPECT_NE(output.find("charge_assignment=cic"), std::string::npos);
    EXPECT_NE(output.find("residual_status: measured"), std::string::npos);
    EXPECT_NE(output.find("Simulation Complete."), std::string::npos);
}

TEST(IntegrationTest, CliRejectsCurrentProfileTableWithoutCustomProfile) {
    const auto nonce = std::chrono::steady_clock::now().time_since_epoch().count();
    const std::filesystem::path tablePath =
        std::filesystem::temp_directory_path() / ("tokamak_profile_table_invalid_" + std::to_string(nonce) + ".csv");
    {
        std::ofstream out(tablePath);
        ASSERT_TRUE(out.is_open());
        out << "0.0,0.0\n";
        out << "1.0,1.0\n";
    }

    const std::string command =
        std::string("\"") + TOKAMAKFUSION_PATH +
        "\" --scenario cold --seed 990 --steps 1 --telemetry-every 1 --no-artifacts --current-profile uniform"
        " --current-profile-table \"" + tablePath.string() + "\" 2>&1 || true";
    const std::string output = RunCommandCapture(command);
    EXPECT_NE(output.find("--current-profile-table requires --current-profile custom"), std::string::npos);

    std::filesystem::remove(tablePath);
}

TEST(IntegrationTest, ArtifactSidecarExistsAndV2HeadersRemainStable) {
    const auto nonce = std::chrono::steady_clock::now().time_since_epoch().count();
    const std::filesystem::path artifactRoot =
        std::filesystem::temp_directory_path() / ("tokamak_m2_artifacts_" + std::to_string(nonce));
    std::filesystem::create_directories(artifactRoot);

    const std::string command = std::string("\"") + TOKAMAKFUSION_PATH +
                                "\" --scenario cold --seed 404 --steps 12 --telemetry-every 6"
                                " --artifact-every 5 --particle-snapshot-every 50 --mag-field-bins 4"
                                " --current-profile parabolic --artifacts-root \"" +
                                artifactRoot.string() + "\"";
    const std::string output = RunCommandCapture(command);
    const std::string runDir = ExtractArtifactDirectory(output);
    ASSERT_FALSE(runDir.empty());

    const std::filesystem::path runPath(runDir);
    const std::filesystem::path sidecarPath = runPath / "magnetic_field_diagnostics_v2.csv";
    ASSERT_TRUE(std::filesystem::exists(sidecarPath));
    ASSERT_TRUE(std::filesystem::exists(runPath / "electrostatic_diagnostics_v2.csv"));

    const std::string sidecarHeader = ReadFirstLine(sidecarPath);
    EXPECT_EQ(
        sidecarHeader,
        "schema_version,step,time_s,bin_index,r_inner_m,r_outer_m,r_center_m,mean_b_t,sample_count,step_max_b_t,recommended_dt_s,profile_kind");
    EXPECT_EQ(CountDataLines(sidecarPath), static_cast<std::size_t>(16));

    const std::string summaryHeader = ReadFirstLine(runPath / "summary_v2.csv");
    EXPECT_EQ(
        summaryHeader,
        "schema_version,step,time_s,active_seed,scenario,total_ions,deuterium,tritium,helium,avg_energy_kev,fusion_events_total,particle_cap_hit_events,rejected_injection_pairs,rejected_fusion_ash,out_of_domain_cell_clamp_events,fusion_attempts,fusion_accepted,max_reactions_in_cell,kinetic_j,beam_injected_j,fusion_alpha_injected_j,total_charge_c");

    const std::string radialHeader = ReadFirstLine(runPath / "radial_profiles_v2.csv");
    EXPECT_EQ(
        radialHeader,
        "schema_version,step,time_s,bin_index,r_inner_m,r_outer_m,r_center_m,ion_count,macro_weight,shell_volume_m3,density_m3,avg_ion_energy_kev,fusion_events_cumulative,fusion_rate_m3_s,fusion_rate_placeholder");

    const std::string manifestText = ReadWholeFile(runPath / "manifest_v2.json");
    EXPECT_NE(manifestText.find("\"magnetic_field_diagnostics_csv\""), std::string::npos);
    EXPECT_NE(manifestText.find("\"electrostatic_diagnostics_csv\""), std::string::npos);

    const std::string electroHeader = ReadFirstLine(runPath / "electrostatic_diagnostics_v2.csv");
    EXPECT_EQ(
        electroHeader,
        "schema_version,step,time_s,electric_field_mode,boundary_condition,charge_assignment,max_electric_field_v_per_m,mean_electric_field_v_per_m,solver_iterations,solver_converged,residual_l2");

    std::filesystem::remove_all(artifactRoot);
}
