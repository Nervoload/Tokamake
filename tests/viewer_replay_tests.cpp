#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <sys/wait.h>

#include "gtest/gtest.h"
#include "tokamak/viewer/replay_loader.hpp"
#include "tokamak/viewer/replay_manifest.hpp"
#include "tokamak/viewer/replay_snapshot.hpp"

namespace {

std::filesystem::path MakeTempDir(const std::string& prefix) {
    const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / (prefix + "_" + std::to_string(now));
    std::filesystem::create_directories(path);
    return path;
}

void WriteTextFile(const std::filesystem::path& path, const std::string& text) {
    std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path);
    ASSERT_TRUE(out.is_open());
    out << text;
    out.flush();
    ASSERT_TRUE(out.good());
}

std::filesystem::path WriteMinimalRunArtifacts(const std::filesystem::path& rootDir) {
    const std::filesystem::path runDir = rootDir / "run_test";
    std::filesystem::create_directories(runDir / "snapshots");

    WriteTextFile(
        runDir / "run_config_v2.json",
        "{\n"
        "  \"schema\": \"tokamak.milestone7.run_config\",\n"
        "  \"schema_version\": 2,\n"
        "  \"scenario\": \"COLD_VACUUM\",\n"
        "  \"seed\": 20260220,\n"
        "  \"tokamak_config\": {\n"
        "    \"major_radius_m\": 2.0,\n"
        "    \"minor_radius_m\": 0.5\n"
        "  }\n"
        "}\n");

    WriteTextFile(
        runDir / "summary_v2.csv",
        "schema_version,step,time_s,active_seed,scenario,total_ions,deuterium,tritium,helium,avg_energy_kev,fusion_events_total,particle_cap_hit_events,rejected_injection_pairs,rejected_fusion_ash,out_of_domain_cell_clamp_events,fusion_attempts,fusion_accepted,max_reactions_in_cell,kinetic_j,beam_injected_j,fusion_alpha_injected_j,total_charge_c\n"
        "2,0,0.0,20260220,COLD_VACUUM,2,1,1,0,0.10,0,0,0,0,0,0,0,0,0,0,0,0.0\n"
        "2,10,0.000001,20260220,COLD_VACUUM,2,1,1,0,0.12,1,0,0,0,0,0,1,1,0,0,0,0.0\n");

    WriteTextFile(
        runDir / "snapshots/particles_step_00000000.csv",
        "schema_version,step,time_s,total_particles,sampled_particles,sample_stride,particle_index,species,species_name,x_m,y_m,z_m,vx_m_per_s,vy_m_per_s,vz_m_per_s,mass_kg,charge_c,q_over_m\n"
        "2,0,0.0,2,2,1,0,0,Deuterium,1.0,2.0,3.0,0,0,0,1,1,1\n"
        "2,0,0.0,2,2,1,1,1,Tritium,4.0,5.0,6.0,0,0,0,1,1,1\n");

    WriteTextFile(
        runDir / "radial_profiles_v2.csv",
        "schema_version,step,time_s,bin_index,r_inner_m,r_outer_m,r_center_m,ion_count,macro_weight,shell_volume_m3,density_m3,avg_ion_energy_kev,fusion_events_cumulative,fusion_rate_m3_s,fusion_rate_placeholder\n");
    WriteTextFile(
        runDir / "magnetic_field_diagnostics_v2.csv",
        "schema_version,step,time_s,bin_index,r_inner_m,r_outer_m,r_center_m,mean_b_t,sample_count,step_max_b_t,recommended_dt_s,profile_kind\n");
    WriteTextFile(
        runDir / "electrostatic_diagnostics_v2.csv",
        "schema_version,step,time_s,electric_field_mode,boundary_condition,charge_assignment,max_electric_field_v_per_m,mean_electric_field_v_per_m,solver_iterations,solver_converged,residual_l2\n");
    WriteTextFile(
        runDir / "speed_histogram_v2.csv",
        "schema_version,step,time_s,bin_index,speed_min_m_per_s,speed_max_m_per_s,count,total_samples\n");
    WriteTextFile(
        runDir / "pitch_angle_histogram_v2.csv",
        "schema_version,step,time_s,bin_index,pitch_min_deg,pitch_max_deg,count,total_samples,invalid_samples\n");
    WriteTextFile(
        runDir / "solver_residuals_v2.csv",
        "schema_version,step,time_s,residual_available,residual_l2,solver_name,status,iterations,converged,tolerance,note\n");

    WriteTextFile(
        runDir / "manifest_v2.json",
        "{\n"
        "  \"schema\": \"tokamak.milestone7.manifest\",\n"
        "  \"schema_version\": 2,\n"
        "  \"run_id\": \"run_test\",\n"
        "  \"created_utc\": \"2026-02-20T00:00:00Z\",\n"
        "  \"files\": {\n"
        "    \"run_config\": \"run_config_v2.json\",\n"
        "    \"summary_csv\": \"summary_v2.csv\",\n"
        "    \"radial_profiles_csv\": \"radial_profiles_v2.csv\",\n"
        "    \"magnetic_field_diagnostics_csv\": \"magnetic_field_diagnostics_v2.csv\",\n"
        "    \"electrostatic_diagnostics_csv\": \"electrostatic_diagnostics_v2.csv\",\n"
        "    \"speed_histogram_csv\": \"speed_histogram_v2.csv\",\n"
        "    \"pitch_histogram_csv\": \"pitch_angle_histogram_v2.csv\",\n"
        "    \"solver_residual_log_csv\": \"solver_residuals_v2.csv\",\n"
        "    \"particle_snapshot_csv_files\": [\"snapshots/particles_step_00000000.csv\"]\n"
        "  }\n"
        "}\n");

    return runDir;
}

int RunCommandCaptureStatus(const std::string& command, std::string* outputOut) {
    std::array<char, 256> buffer{};
    std::string output;

    FILE* pipe = popen(command.c_str(), "r");
    if (pipe == nullptr) {
        if (outputOut != nullptr) {
            *outputOut = "popen failed";
        }
        return -1;
    }

    while (fgets(buffer.data(), static_cast<int>(buffer.size()), pipe) != nullptr) {
        output += buffer.data();
    }

    const int rawStatus = pclose(pipe);
    if (outputOut != nullptr) {
        *outputOut = output;
    }

    if (rawStatus == -1) {
        return -1;
    }

    if (WIFEXITED(rawStatus)) {
        return WEXITSTATUS(rawStatus);
    }
    return rawStatus;
}

}  // namespace

TEST(ViewerReplayTest, ViewerManifestV2ParsesRequiredKeys) {
    const std::filesystem::path tempDir = MakeTempDir("viewer_manifest_parse");
    const std::filesystem::path runDir = WriteMinimalRunArtifacts(tempDir);

    tokamak::viewer::ReplayManifest manifest;
    std::string error;
    EXPECT_TRUE(tokamak::viewer::ParseManifestV2File(runDir / "manifest_v2.json", &manifest, &error));
    EXPECT_EQ(manifest.schemaVersion, 2u);
    EXPECT_EQ(manifest.runId, "run_test");
    EXPECT_EQ(manifest.files.summaryCsv, "summary_v2.csv");
    EXPECT_EQ(manifest.files.particleSnapshotCsvFiles.size(), static_cast<std::size_t>(1));
}

TEST(ViewerReplayTest, ViewerManifestV2FailsWhenParticleSnapshotsKeyMissing) {
    const std::filesystem::path tempDir = MakeTempDir("viewer_manifest_missing");
    const std::filesystem::path manifestPath = tempDir / "manifest_v2.json";

    WriteTextFile(
        manifestPath,
        "{\n"
        "  \"schema_version\": 2,\n"
        "  \"run_id\": \"run_test\",\n"
        "  \"created_utc\": \"2026-02-20T00:00:00Z\",\n"
        "  \"files\": {\n"
        "    \"run_config\": \"run_config_v2.json\",\n"
        "    \"summary_csv\": \"summary_v2.csv\",\n"
        "    \"radial_profiles_csv\": \"radial_profiles_v2.csv\",\n"
        "    \"magnetic_field_diagnostics_csv\": \"magnetic_field_diagnostics_v2.csv\",\n"
        "    \"electrostatic_diagnostics_csv\": \"electrostatic_diagnostics_v2.csv\",\n"
        "    \"speed_histogram_csv\": \"speed_histogram_v2.csv\",\n"
        "    \"pitch_histogram_csv\": \"pitch_angle_histogram_v2.csv\",\n"
        "    \"solver_residual_log_csv\": \"solver_residuals_v2.csv\"\n"
        "  }\n"
        "}\n");

    tokamak::viewer::ReplayManifest manifest;
    std::string error;
    EXPECT_FALSE(tokamak::viewer::ParseManifestV2File(manifestPath, &manifest, &error));
    EXPECT_NE(error.find("particle_snapshot_csv_files"), std::string::npos);
}

TEST(ViewerReplayTest, ViewerSnapshotCsvParsesSpeciesAndPositions) {
    const std::filesystem::path tempDir = MakeTempDir("viewer_snapshot_ok");
    const std::filesystem::path snapshotPath = tempDir / "particles_step_00000000.csv";

    WriteTextFile(
        snapshotPath,
        "schema_version,step,time_s,total_particles,sampled_particles,sample_stride,particle_index,species,species_name,x_m,y_m,z_m\n"
        "2,0,0.0,2,2,1,0,0,Deuterium,1.0,2.0,3.0\n"
        "2,0,0.0,2,2,1,1,1,Tritium,4.0,5.0,6.0\n");

    tokamak::viewer::ReplayFrame frame;
    std::string error;
    EXPECT_TRUE(tokamak::viewer::ParseParticleSnapshotCsv(snapshotPath, &frame, &error));
    EXPECT_EQ(frame.step, 0);
    EXPECT_EQ(frame.particles.size(), static_cast<std::size_t>(2));
    EXPECT_EQ(frame.particles[0].species, tokamak::viewer::ReplaySpecies::Deuterium);
    EXPECT_NEAR(frame.particles[1].position_m.z, 6.0f, 1.0e-6f);
}

TEST(ViewerReplayTest, ViewerSnapshotCsvFailsOnHeaderMismatch) {
    const std::filesystem::path tempDir = MakeTempDir("viewer_snapshot_header_fail");
    const std::filesystem::path snapshotPath = tempDir / "particles_step_00000000.csv";

    WriteTextFile(
        snapshotPath,
        "schema_version,step,time_s,total_particles,sampled_particles,sample_stride,particle_index,species,x_m,y_m,z_m\n"
        "2,0,0.0,2,2,1,0,0,1.0,2.0,3.0\n");

    tokamak::viewer::ReplayFrame frame;
    std::string error;
    EXPECT_FALSE(tokamak::viewer::ParseParticleSnapshotCsv(snapshotPath, &frame, &error));
    EXPECT_NE(error.find("species_name"), std::string::npos);
}

TEST(ViewerReplayTest, ViewerLoaderResolvesRelativePathsFromManifest) {
    const std::filesystem::path tempDir = MakeTempDir("viewer_loader_resolve");
    const std::filesystem::path runDir = WriteMinimalRunArtifacts(tempDir);

    tokamak::viewer::ReplayLoader loader;
    EXPECT_TRUE(loader.OpenFromManifest(runDir / "manifest_v2.json"));
    EXPECT_TRUE(loader.HasData());
    EXPECT_EQ(loader.FrameCount(), static_cast<std::size_t>(1));

    tokamak::viewer::ReplayFrame frame;
    EXPECT_TRUE(loader.LoadFrameByOrderedIndex(0, &frame));
    EXPECT_EQ(frame.step, 0);
    EXPECT_EQ(frame.particles.size(), static_cast<std::size_t>(2));

    const tokamak::viewer::ReplaySummaryPoint* summary = loader.SummaryForStep(0);
    ASSERT_NE(summary, nullptr);
    EXPECT_NEAR(summary->avgEnergy_keV, 0.10, 1.0e-12);
}

TEST(ViewerReplayTest, ViewerLoaderRejectsMissingSnapshotFile) {
    const std::filesystem::path tempDir = MakeTempDir("viewer_loader_missing_snapshot");
    const std::filesystem::path runDir = tempDir / "run_test";
    std::filesystem::create_directories(runDir);

    WriteTextFile(
        runDir / "run_config_v2.json",
        "{\"schema_version\":2,\"scenario\":\"COLD_VACUUM\",\"tokamak_config\":{\"major_radius_m\":2.0,\"minor_radius_m\":0.5}}\n");
    WriteTextFile(
        runDir / "summary_v2.csv",
        "schema_version,step,time_s,total_ions,avg_energy_kev,fusion_events_total\n"
        "2,0,0.0,0,0.0,0\n");

    WriteTextFile(
        runDir / "manifest_v2.json",
        "{\n"
        "  \"schema\": \"tokamak.milestone7.manifest\",\n"
        "  \"schema_version\": 2,\n"
        "  \"run_id\": \"run_test\",\n"
        "  \"created_utc\": \"2026-02-20T00:00:00Z\",\n"
        "  \"files\": {\n"
        "    \"run_config\": \"run_config_v2.json\",\n"
        "    \"summary_csv\": \"summary_v2.csv\",\n"
        "    \"radial_profiles_csv\": \"radial_profiles_v2.csv\",\n"
        "    \"magnetic_field_diagnostics_csv\": \"magnetic_field_diagnostics_v2.csv\",\n"
        "    \"electrostatic_diagnostics_csv\": \"electrostatic_diagnostics_v2.csv\",\n"
        "    \"speed_histogram_csv\": \"speed_histogram_v2.csv\",\n"
        "    \"pitch_histogram_csv\": \"pitch_angle_histogram_v2.csv\",\n"
        "    \"solver_residual_log_csv\": \"solver_residuals_v2.csv\",\n"
        "    \"particle_snapshot_csv_files\": [\"snapshots/particles_step_00000000.csv\"]\n"
        "  }\n"
        "}\n");

    tokamak::viewer::ReplayLoader loader;
    EXPECT_FALSE(loader.OpenFromManifest(runDir / "manifest_v2.json"));
    EXPECT_NE(loader.LastError().find("does not exist"), std::string::npos);
}

TEST(ViewerReplayTest, CMakeBuildViewerOffDoesNotRequireImGuiGlad) {
    const std::filesystem::path tempDir = MakeTempDir("viewer_cmake_off");
    const std::filesystem::path buildDir = tempDir / "build_off";

    std::ostringstream command;
    command << "cmake -S \"" << TOKAMAK_SOURCE_DIR << "\" -B \"" << buildDir.string()
            << "\" -DBUILD_VIEWER=OFF -DBUILD_TESTING=OFF 2>&1";

    std::string output;
    const int status = RunCommandCaptureStatus(command.str(), &output);
    EXPECT_EQ(status, 0);
}

TEST(ViewerReplayTest, CMakeBuildViewerOnFailsWithClearMessageWhenVendorFilesMissing) {
    const std::filesystem::path sourceDir = TOKAMAK_SOURCE_DIR;
    const std::vector<std::filesystem::path> requiredVendorFiles = {
        sourceDir / "third_party/glad/src/glad.c",
        sourceDir / "third_party/glad/include/glad/glad.h",
        sourceDir / "third_party/glad/include/KHR/khrplatform.h",
        sourceDir / "third_party/imgui/imgui.cpp",
        sourceDir / "third_party/imgui/imgui_draw.cpp",
        sourceDir / "third_party/imgui/imgui_tables.cpp",
        sourceDir / "third_party/imgui/imgui_widgets.cpp",
        sourceDir / "third_party/imgui/backends/imgui_impl_glfw.cpp",
        sourceDir / "third_party/imgui/backends/imgui_impl_opengl3.cpp",
    };

    bool missingVendorFile = false;
    for (const auto& path : requiredVendorFiles) {
        if (!std::filesystem::exists(path)) {
            missingVendorFile = true;
            break;
        }
    }

    if (!missingVendorFile) {
        EXPECT_TRUE(true);
        return;
    }

    const std::filesystem::path tempDir = MakeTempDir("viewer_cmake_on_missing");
    const std::filesystem::path buildDir = tempDir / "build_on";

    std::ostringstream command;
    command << "cmake -S \"" << TOKAMAK_SOURCE_DIR << "\" -B \"" << buildDir.string()
            << "\" -DBUILD_VIEWER=ON -DVIEWER_STRICT_VENDOR_CHECK=ON -DBUILD_TESTING=OFF 2>&1";

    std::string output;
    const int status = RunCommandCaptureStatus(command.str(), &output);
    EXPECT_NE(status, 0);
    EXPECT_NE(output.find("requires vendored ImGui/GLAD files"), std::string::npos);
}
