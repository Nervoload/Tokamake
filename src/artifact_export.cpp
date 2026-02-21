#include "tokamak/artifact_export.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <sstream>

namespace {

constexpr double kPi = 3.14159265358979323846;

std::string JsonEscape(const std::string& text) {
    std::string out;
    out.reserve(text.size() + 8);
    for (const char ch : text) {
        switch (ch) {
            case '\\':
                out += "\\\\";
                break;
            case '\"':
                out += "\\\"";
                break;
            case '\n':
                out += "\\n";
                break;
            case '\r':
                out += "\\r";
                break;
            case '\t':
                out += "\\t";
                break;
            default:
                out.push_back(ch);
                break;
        }
    }
    return out;
}

std::string TimestampForPathUtc() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t nowTime = std::chrono::system_clock::to_time_t(now);
    std::tm utc{};
#if defined(_WIN32)
    gmtime_s(&utc, &nowTime);
#else
    utc = *std::gmtime(&nowTime);
#endif
    std::ostringstream oss;
    oss << std::put_time(&utc, "%Y%m%dT%H%M%SZ");
    return oss.str();
}

uint64_t UnixNowNs() {
    const auto now = std::chrono::system_clock::now().time_since_epoch();
    return static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(now).count());
}

std::string TimestampIsoUtc() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t nowTime = std::chrono::system_clock::to_time_t(now);
    std::tm utc{};
#if defined(_WIN32)
    gmtime_s(&utc, &nowTime);
#else
    utc = *std::gmtime(&nowTime);
#endif
    std::ostringstream oss;
    oss << std::put_time(&utc, "%Y-%m-%dT%H:%M:%SZ");
    return oss.str();
}

double MinorRadiusMeters(const tokamak::TokamakConfig& config, const tokamak::Vec3& position) {
    const double majorR = std::sqrt(
        static_cast<double>(position.x) * static_cast<double>(position.x) +
        static_cast<double>(position.y) * static_cast<double>(position.y));
    const double radial = majorR - static_cast<double>(config.majorRadius_m);
    return std::sqrt(radial * radial + static_cast<double>(position.z) * static_cast<double>(position.z));
}

tokamak::Vec3 ApproximateBField(
    const tokamak::TokamakConfig& config,
    const tokamak::PlasmaCurrentProfileConfig& profileConfig,
    const tokamak::Vec3& position) {
    return tokamak::EvaluateMagneticFieldSample(config, profileConfig, position).totalField_T;
}

const char* SpeciesName(tokamak::ParticleType species) {
    switch (species) {
        case tokamak::ParticleType::Deuterium:
            return "Deuterium";
        case tokamak::ParticleType::Tritium:
            return "Tritium";
        case tokamak::ParticleType::Helium:
            return "Helium";
        case tokamak::ParticleType::Dead:
            return "Dead";
    }
    return "Unknown";
}

const char* ResidualSolverName(tokamak::ResidualSolverKind kind) {
    switch (kind) {
        case tokamak::ResidualSolverKind::None:
            return "none";
        case tokamak::ResidualSolverKind::Sor:
            return "sor";
        case tokamak::ResidualSolverKind::Other:
            return "other";
    }
    return "unknown";
}

const char* ResidualStatusName(tokamak::ResidualStatus status) {
    switch (status) {
        case tokamak::ResidualStatus::Placeholder:
            return "placeholder";
        case tokamak::ResidualStatus::Measured:
            return "measured";
        case tokamak::ResidualStatus::Unavailable:
            return "unavailable";
        case tokamak::ResidualStatus::Failed:
            return "failed";
    }
    return "unknown";
}

const char* ResidualStatusNote(const tokamak::TelemetrySnapshot& telemetry) {
    switch (telemetry.solverResidual.status) {
        case tokamak::ResidualStatus::Placeholder:
            return "Electrostatic solver residual unavailable in current field model";
        case tokamak::ResidualStatus::Measured:
            return "Residual measured by active electrostatic solver";
        case tokamak::ResidualStatus::Unavailable:
            return "Solver exists but did not provide residual for this step";
        case tokamak::ResidualStatus::Failed:
            return "Solver attempted but residual is invalid due to solver failure";
    }
    return "Residual status unknown";
}

}  // namespace

namespace tokamak {

Milestone7ArtifactExporter::~Milestone7ArtifactExporter() {
    CloseFiles();
}

bool Milestone7ArtifactExporter::Start(
    const RunConfig& runConfig,
    const TokamakEngine& engine,
    const ArtifactExportConfig& config) {
    CloseFiles();
    lastError_.clear();
    runConfig_ = runConfig;
    config_ = config;

    config_.metricsEveryNSteps = std::max(1, config_.metricsEveryNSteps);
    config_.particleSnapshotEveryNSteps = std::max(1, config_.particleSnapshotEveryNSteps);
    config_.maxParticlesPerSnapshot = std::max<std::size_t>(1, config_.maxParticlesPerSnapshot);
    config_.radialBinCount = std::max<std::size_t>(1, config_.radialBinCount);
    config_.speedHistogramBinCount = std::max<std::size_t>(1, config_.speedHistogramBinCount);
    config_.pitchHistogramBinCount = std::max<std::size_t>(1, config_.pitchHistogramBinCount);

    const std::string timestamp = TimestampForPathUtc();
    std::ostringstream runId;
    runId << "run_" << timestamp << "_"
          << ScenarioName(runConfig.scenario) << "_seed" << engine.ActiveSeed()
          << "_ns" << UnixNowNs();
    runId_ = runId.str();

    const std::size_t fusionBinCount = engine.FusionAcceptedByRadiusBins().size();
    if (fusionBinCount > 0) {
        config_.radialBinCount = fusionBinCount;
    }

    const std::filesystem::path rootPath(config_.outputRootDirectory.empty() ? "output/runs" : config_.outputRootDirectory);
    const std::filesystem::path runDirectoryPath = rootPath / runId_;
    const std::filesystem::path snapshotDirectoryPath = runDirectoryPath / "snapshots";

    std::error_code ec;
    std::filesystem::create_directories(snapshotDirectoryPath, ec);
    if (ec) {
        SetError("Failed to create artifact output directory: " + runDirectoryPath.string());
        return false;
    }

    runDirectory_ = runDirectoryPath.string();

    runConfigRelativePath_ = "run_config_v2.json";
    summaryRelativePath_ = "summary_v2.csv";
    radialRelativePath_ = "radial_profiles_v2.csv";
    magneticFieldRelativePath_ = "magnetic_field_diagnostics_v2.csv";
    electrostaticDiagnosticsRelativePath_ = "electrostatic_diagnostics_v2.csv";
    fusionReactivityDiagnosticsRelativePath_ = "fusion_reactivity_diagnostics_v2.csv";
    wallInteractionBridgeRelativePath_ = "wall_interaction_bridge_v2.csv";
    speedHistogramRelativePath_ = "speed_histogram_v2.csv";
    pitchHistogramRelativePath_ = "pitch_angle_histogram_v2.csv";
    solverResidualRelativePath_ = "solver_residuals_v2.csv";

    summaryCsv_.open((runDirectoryPath / summaryRelativePath_).string(), std::ios::trunc);
    radialCsv_.open((runDirectoryPath / radialRelativePath_).string(), std::ios::trunc);
    magneticFieldCsv_.open((runDirectoryPath / magneticFieldRelativePath_).string(), std::ios::trunc);
    electrostaticDiagnosticsCsv_.open((runDirectoryPath / electrostaticDiagnosticsRelativePath_).string(), std::ios::trunc);
    fusionReactivityDiagnosticsCsv_.open((runDirectoryPath / fusionReactivityDiagnosticsRelativePath_).string(), std::ios::trunc);
    wallInteractionBridgeCsv_.open((runDirectoryPath / wallInteractionBridgeRelativePath_).string(), std::ios::trunc);
    speedHistogramCsv_.open((runDirectoryPath / speedHistogramRelativePath_).string(), std::ios::trunc);
    pitchHistogramCsv_.open((runDirectoryPath / pitchHistogramRelativePath_).string(), std::ios::trunc);
    solverResidualCsv_.open((runDirectoryPath / solverResidualRelativePath_).string(), std::ios::trunc);

    if (!summaryCsv_.is_open() || !radialCsv_.is_open() || !magneticFieldCsv_.is_open() ||
        !electrostaticDiagnosticsCsv_.is_open() ||
        !fusionReactivityDiagnosticsCsv_.is_open() ||
        !wallInteractionBridgeCsv_.is_open() ||
        !speedHistogramCsv_.is_open() ||
        !pitchHistogramCsv_.is_open() || !solverResidualCsv_.is_open()) {
        SetError("Failed to open one or more artifact CSV files under: " + runDirectory_);
        return false;
    }

    summaryCsv_
        << "schema_version,step,time_s,active_seed,scenario,total_ions,deuterium,tritium,helium,avg_energy_kev,"
        << "fusion_events_total,particle_cap_hit_events,rejected_injection_pairs,rejected_fusion_ash,out_of_domain_cell_clamp_events,fusion_attempts,"
        << "fusion_accepted,max_reactions_in_cell,kinetic_j,beam_injected_j,fusion_alpha_injected_j,total_charge_c\n";

    radialCsv_
        << "schema_version,step,time_s,bin_index,r_inner_m,r_outer_m,r_center_m,ion_count,macro_weight,shell_volume_m3,"
        << "density_m3,avg_ion_energy_kev,fusion_events_cumulative,fusion_rate_m3_s,fusion_rate_placeholder\n";

    magneticFieldCsv_
        << "schema_version,step,time_s,bin_index,r_inner_m,r_outer_m,r_center_m,mean_b_t,sample_count,"
        << "step_max_b_t,recommended_dt_s,profile_kind\n";

    electrostaticDiagnosticsCsv_
        << "schema_version,step,time_s,electric_field_mode,boundary_condition,charge_assignment,max_electric_field_v_per_m,"
        << "mean_electric_field_v_per_m,solver_iterations,solver_converged,residual_l2\n";

    fusionReactivityDiagnosticsCsv_
        << "schema_version,step,time_s,reactivity_model,cross_section_scale,probability_clamp,min_energy_kev,"
        << "fusion_attempts_step,fusion_accepted_step,fusion_weight_attempted_step,fusion_weight_accepted_step,"
        << "fuel_weight_consumed_d_step,fuel_weight_consumed_t_step,ash_weight_produced_he_step,avg_sigma_m2_step,"
        << "avg_probability_step,avg_relative_speed_m_per_s_step,fusion_attempts_total,fusion_accepted_total\n";

    wallInteractionBridgeCsv_
        << "schema_version,step,time_s,wall_mode,recycle_fraction,wall_hit_count_step,wall_impact_energy_j_step,"
        << "wall_loss_weight_step,wall_hit_count_total,wall_impact_energy_j_total,wall_loss_weight_total\n";

    speedHistogramCsv_
        << "schema_version,step,time_s,bin_index,speed_min_m_per_s,speed_max_m_per_s,count,total_samples\n";

    pitchHistogramCsv_
        << "schema_version,step,time_s,bin_index,pitch_min_deg,pitch_max_deg,count,total_samples,invalid_samples\n";

    solverResidualCsv_
        << "schema_version,step,time_s,residual_available,residual_l2,solver_name,status,iterations,converged,tolerance,note\n";

    if (!summaryCsv_.good() || !radialCsv_.good() || !magneticFieldCsv_.good() || !electrostaticDiagnosticsCsv_.good() ||
        !fusionReactivityDiagnosticsCsv_.good() || !wallInteractionBridgeCsv_.good() ||
        !speedHistogramCsv_.good() ||
        !pitchHistogramCsv_.good() || !solverResidualCsv_.good()) {
        SetError("Failed to write CSV headers for artifact files");
        return false;
    }

    if (!WriteRunConfigJson(runConfig, engine)) {
        return false;
    }

    particleSnapshotRelativePaths_.clear();
    lastMetricsStepWritten_ = -1;
    lastParticleSnapshotStepWritten_ = -1;
    started_ = true;
    return true;
}

bool Milestone7ArtifactExporter::WriteStep(
    const TokamakEngine& engine,
    const TelemetrySnapshot& telemetry,
    bool forceWrite) {
    if (!started_) {
        SetError("Artifact exporter was not started");
        return false;
    }

    const bool metricsDue =
        forceWrite || (telemetry.step % config_.metricsEveryNSteps == 0);
    const bool particleSnapshotDue =
        forceWrite || (telemetry.step % config_.particleSnapshotEveryNSteps == 0);

    if (metricsDue && telemetry.step != lastMetricsStepWritten_) {
        if (!WriteSummaryCsvRow(runConfig_, telemetry) ||
            !WriteRadialProfileRows(engine, telemetry) ||
            !WriteMagneticFieldDiagnosticsRows(engine, telemetry) ||
            !WriteElectrostaticDiagnosticsRow(runConfig_, telemetry) ||
            !WriteFusionReactivityDiagnosticsRow(runConfig_, telemetry) ||
            !WriteWallInteractionBridgeRow(runConfig_, telemetry) ||
            !WriteSpeedHistogramRows(engine, telemetry) ||
            !WritePitchHistogramRows(engine, telemetry) ||
            !WriteSolverResidualRow(telemetry)) {
            return false;
        }
        lastMetricsStepWritten_ = telemetry.step;
    }

    if (particleSnapshotDue && telemetry.step != lastParticleSnapshotStepWritten_) {
        if (!WriteParticleSnapshotCsv(engine, telemetry)) {
            return false;
        }
        lastParticleSnapshotStepWritten_ = telemetry.step;
    }

    return true;
}

bool Milestone7ArtifactExporter::Finish() {
    if (!started_) {
        return true;
    }

    summaryCsv_.flush();
    radialCsv_.flush();
    magneticFieldCsv_.flush();
    electrostaticDiagnosticsCsv_.flush();
    fusionReactivityDiagnosticsCsv_.flush();
    wallInteractionBridgeCsv_.flush();
    speedHistogramCsv_.flush();
    pitchHistogramCsv_.flush();
    solverResidualCsv_.flush();

    if (!summaryCsv_.good() || !radialCsv_.good() || !magneticFieldCsv_.good() || !electrostaticDiagnosticsCsv_.good() ||
        !fusionReactivityDiagnosticsCsv_.good() || !wallInteractionBridgeCsv_.good() ||
        !speedHistogramCsv_.good() ||
        !pitchHistogramCsv_.good() || !solverResidualCsv_.good()) {
        SetError("Failed to flush artifact CSV streams");
        return false;
    }

    if (!WriteManifestJson()) {
        return false;
    }

    CloseFiles();
    started_ = false;
    return true;
}

bool Milestone7ArtifactExporter::WriteRunConfigJson(const RunConfig& runConfig, const TokamakEngine& engine) {
    const std::filesystem::path path = std::filesystem::path(runDirectory_) / runConfigRelativePath_;
    std::ofstream out(path.string(), std::ios::trunc);
    if (!out.is_open()) {
        SetError("Failed to open run config JSON path: " + path.string());
        return false;
    }

    const TokamakConfig& tokamakConfig = engine.Config();
    const NBIConfig& nbi = engine.BeamConfig();

    out << "{\n";
    out << "  \"schema\": \"tokamak.milestone7.run_config\",\n";
    out << "  \"schema_version\": " << kMilestone7OutputSchemaVersion << ",\n";
    out << "  \"run_id\": \"" << JsonEscape(runId_) << "\",\n";
    out << "  \"created_utc\": \"" << TimestampIsoUtc() << "\",\n";
    out << "  \"scenario\": \"" << ScenarioName(runConfig.scenario) << "\",\n";
    out << "  \"seed\": " << engine.ActiveSeed() << ",\n";
    out << "  \"requested_seed\": ";
    if (runConfig.seed.has_value()) {
        out << runConfig.seed.value();
    } else {
        out << "null";
    }
    out << ",\n";
    out << "  \"time_step_s\": " << std::setprecision(std::numeric_limits<double>::max_digits10)
        << static_cast<double>(runConfig.timeStep_s) << ",\n";
    out << "  \"total_steps\": " << runConfig.totalSteps << ",\n";
    out << "  \"telemetry_every_n_steps\": " << runConfig.telemetryEveryNSteps << ",\n";
    out << "  \"particle_cap\": " << runConfig.particleCap << ",\n";
    out << "  \"tokamak_config\": {\n";
    out << "    \"major_radius_m\": " << tokamakConfig.majorRadius_m << ",\n";
    out << "    \"minor_radius_m\": " << tokamakConfig.minorRadius_m << ",\n";
    out << "    \"toroidal_current_a\": " << tokamakConfig.toroidalCurrent_A << ",\n";
    out << "    \"toroidal_coil_turns\": " << tokamakConfig.toroidalCoilTurns << ",\n";
    out << "    \"plasma_current_a\": " << tokamakConfig.plasmaCurrent_A << "\n";
    out << "  },\n";
    out << "  \"nbi_config\": {\n";
    out << "    \"is_active\": " << (nbi.isActive ? "true" : "false") << ",\n";
    out << "    \"beam_energy_kev\": " << nbi.beamEnergy_keV << ",\n";
    out << "    \"particles_per_step\": " << nbi.particlesPerStep << ",\n";
    out << "    \"injector_pos_m\": [" << nbi.injectorPos.x << ", " << nbi.injectorPos.y << ", " << nbi.injectorPos.z << "],\n";
    out << "    \"injection_normal\": ["
        << nbi.injectionNormal.x << ", " << nbi.injectionNormal.y << ", " << nbi.injectionNormal.z << "]\n";
    out << "  },\n";
    out << "  \"artifact_settings\": {\n";
    out << "    \"output_root_directory\": \"" << JsonEscape(config_.outputRootDirectory) << "\",\n";
    out << "    \"metrics_every_n_steps\": " << config_.metricsEveryNSteps << ",\n";
    out << "    \"particle_snapshot_every_n_steps\": " << config_.particleSnapshotEveryNSteps << ",\n";
    out << "    \"max_particles_per_snapshot\": " << config_.maxParticlesPerSnapshot << ",\n";
    out << "    \"radial_bin_count\": " << config_.radialBinCount << ",\n";
    out << "    \"speed_histogram_bin_count\": " << config_.speedHistogramBinCount << ",\n";
    out << "    \"pitch_histogram_bin_count\": " << config_.pitchHistogramBinCount << "\n";
    out << "  },\n";
    out << "  \"magnetic_field_settings\": {\n";
    out << "    \"current_profile_kind\": \""
        << PlasmaCurrentProfileKindName(engine.PlasmaCurrentProfile().kind) << "\",\n";
    out << "    \"current_profile_axis_epsilon_m\": " << engine.PlasmaCurrentProfile().axisEpsilon_m << ",\n";
    out << "    \"current_profile_axis_blend_m\": " << engine.PlasmaCurrentProfile().customAxisBlendRadius_m << ",\n";
    out << "    \"configured_radial_bin_count\": " << runConfig.magneticFieldRadialBinCount << ",\n";
    out << "    \"configured_dt_safety_fraction\": " << runConfig.magneticFieldDtSafetyFraction << ",\n";
    out << "    \"custom_profile_point_count\": " << engine.PlasmaCurrentProfile().customTable.size() << "\n";
    out << "  },\n";
    out << "  \"electrostatic_settings\": {\n";
    out << "    \"electric_field_mode\": \"" << ElectricFieldModeName(runConfig.electricFieldMode) << "\",\n";
    out << "    \"boundary_condition\": \""
        << ElectrostaticBoundaryConditionName(runConfig.electrostaticBoundaryCondition) << "\",\n";
    out << "    \"charge_assignment_scheme\": \"" << ChargeAssignmentSchemeName(runConfig.chargeAssignmentScheme)
        << "\",\n";
    out << "    \"grid_bin_count\": " << runConfig.electrostaticGridBinCount << ",\n";
    out << "    \"solver_tolerance\": " << runConfig.electrostaticSolverTolerance << ",\n";
    out << "    \"solver_max_iterations\": " << runConfig.electrostaticSolverMaxIterations << ",\n";
    out << "    \"sor_omega\": " << runConfig.electrostaticSorOmega << ",\n";
    out << "    \"neutralizing_background_fraction\": " << runConfig.electrostaticNeutralizingBackgroundFraction << "\n";
    out << "  },\n";
    out << "  \"fusion_settings\": {\n";
    out << "    \"reactivity_model\": \"" << FusionReactivityModelKindName(runConfig.fusionReactivityModelKind)
        << "\",\n";
    out << "    \"cross_section_scale\": " << runConfig.fusionCrossSectionScale << ",\n";
    out << "    \"probability_clamp\": " << runConfig.fusionProbabilityClamp << ",\n";
    out << "    \"min_energy_kev\": " << runConfig.fusionMinEnergy_keV << ",\n";
    out << "    \"diagnostics_radial_bins\": " << runConfig.fusionDiagnosticsRadialBins << "\n";
    out << "  },\n";
    out << "  \"wall_settings\": {\n";
    out << "    \"boundary_mode\": \"" << WallBoundaryModeName(runConfig.wallBoundaryMode) << "\",\n";
    out << "    \"recycle_fraction\": " << runConfig.recycleFraction << "\n";
    out << "  }\n";
    out << "}\n";

    out.flush();
    if (!out.good()) {
        SetError("Failed to write run config JSON file: " + path.string());
        return false;
    }

    return true;
}

bool Milestone7ArtifactExporter::WriteSummaryCsvRow(const RunConfig& runConfig, const TelemetrySnapshot& telemetry) {
    summaryCsv_ << std::setprecision(std::numeric_limits<double>::max_digits10)
                << kMilestone7OutputSchemaVersion << ','
                << telemetry.step << ','
                << telemetry.time_s << ','
                << telemetry.activeSeed << ','
                << ScenarioName(runConfig.scenario) << ','
                << telemetry.species.AliveCount() << ','
                << telemetry.species.deuterium << ','
                << telemetry.species.tritium << ','
                << telemetry.species.helium << ','
                << telemetry.avgEnergy_keV << ','
                << telemetry.fusionEvents << ','
                << telemetry.counters.particleCapHitEvents << ','
                << telemetry.counters.rejectedInjectionPairs << ','
                << telemetry.counters.rejectedFusionAsh << ','
                << telemetry.counters.outOfDomainCellClampEvents << ','
                << telemetry.counters.fusionAttempts << ','
                << telemetry.counters.fusionAccepted << ','
                << telemetry.counters.maxReactionsInCell << ','
                << telemetry.budget.kinetic_J << ','
                << telemetry.budget.beamInjected_J << ','
                << telemetry.budget.fusionAlphaInjected_J << ','
                << telemetry.budget.totalCharge_C << '\n';

    if (!summaryCsv_.good()) {
        SetError("Failed writing summary CSV row");
        return false;
    }
    return true;
}

bool Milestone7ArtifactExporter::WriteRadialProfileRows(
    const TokamakEngine& engine,
    const TelemetrySnapshot& telemetry) {
    const TokamakConfig& config = engine.Config();
    const ParticleSystem& particles = engine.Particles();
    const std::size_t binCount = config_.radialBinCount;
    const double minorRadius_m = static_cast<double>(config.minorRadius_m);
    const double majorRadius_m = static_cast<double>(config.majorRadius_m);

    std::vector<uint64_t> ionCounts(binCount, 0);
    std::vector<double> energySumKeV(binCount, 0.0);

    const auto& positions = particles.Positions();
    const auto& velocities = particles.Velocities();
    const auto& masses = particles.Masses();
    const auto& species = particles.Species();
    for (std::size_t i = 0; i < particles.Size(); ++i) {
        if (species[i] == ParticleType::Dead) {
            continue;
        }

        const double minorRadius = MinorRadiusMeters(config, positions[i]);
        const double clampedRadius = std::max(0.0, std::min(minorRadius, minorRadius_m));
        const double normalized = (minorRadius_m > 0.0) ? (clampedRadius / minorRadius_m) : 0.0;
        std::size_t bin = static_cast<std::size_t>(normalized * static_cast<double>(binCount));
        if (bin >= binCount) {
            bin = binCount - 1;
        }

        const double speedSquared = static_cast<double>(Vec3::Dot(velocities[i], velocities[i]));
        const double kineticEnergyJ = 0.5 * static_cast<double>(masses[i]) * speedSquared;
        const double kineticEnergyKeV = kineticEnergyJ / (1000.0 * static_cast<double>(constants::kElementaryCharge_C));

        ++ionCounts[bin];
        energySumKeV[bin] += kineticEnergyKeV;
    }

    const std::vector<uint64_t>& fusionBins = engine.FusionAcceptedByRadiusBins();
    const bool fusionBinCompatible = (fusionBins.size() == binCount);
    const double macroWeight = static_cast<double>(particles.macroWeight);
    const double elapsedSeconds = std::max(telemetry.time_s, 0.0);

    for (std::size_t bin = 0; bin < binCount; ++bin) {
        const double normLower = static_cast<double>(bin) / static_cast<double>(binCount);
        const double normUpper = static_cast<double>(bin + 1) / static_cast<double>(binCount);
        const double rInner = minorRadius_m * normLower;
        const double rOuter = minorRadius_m * normUpper;
        const double rCenter = 0.5 * (rInner + rOuter);

        const double shellVolume = 2.0 * kPi * kPi * majorRadius_m * (rOuter * rOuter - rInner * rInner);
        const double weightedCount = static_cast<double>(ionCounts[bin]) * macroWeight;
        const double density = (shellVolume > 0.0) ? (weightedCount / shellVolume) : 0.0;
        const double avgEnergyKeV = (ionCounts[bin] > 0)
            ? (energySumKeV[bin] / static_cast<double>(ionCounts[bin]))
            : 0.0;

        uint64_t fusionCount = 0;
        bool fusionRatePlaceholder = false;
        double fusionRate = 0.0;
        if (fusionBinCompatible) {
            fusionCount = fusionBins[bin];
            if (elapsedSeconds > 0.0 && shellVolume > 0.0) {
                fusionRate =
                    (static_cast<double>(fusionCount) * macroWeight) /
                    (elapsedSeconds * shellVolume);
            } else {
                fusionRate = 0.0;
            }
        } else {
            fusionRatePlaceholder = true;
            fusionRate = std::numeric_limits<double>::quiet_NaN();
        }

        radialCsv_ << std::setprecision(std::numeric_limits<double>::max_digits10)
                   << kMilestone7OutputSchemaVersion << ','
                   << telemetry.step << ','
                   << telemetry.time_s << ','
                   << bin << ','
                   << rInner << ','
                   << rOuter << ','
                   << rCenter << ','
                   << ionCounts[bin] << ','
                   << macroWeight << ','
                   << shellVolume << ','
                   << density << ','
                   << avgEnergyKeV << ','
                   << fusionCount << ','
                   << fusionRate << ','
                   << (fusionRatePlaceholder ? 1 : 0) << '\n';
    }

    if (!radialCsv_.good()) {
        SetError("Failed writing radial profile CSV rows");
        return false;
    }
    return true;
}

bool Milestone7ArtifactExporter::WriteMagneticFieldDiagnosticsRows(
    const TokamakEngine& engine,
    const TelemetrySnapshot& telemetry) {
    const auto& radialMean = telemetry.magneticField.radialMeanField_T;
    const auto& radialCounts = telemetry.magneticField.radialSampleCounts;
    const std::size_t binCount = std::max(radialMean.size(), radialCounts.size());
    if (binCount == 0) {
        return true;
    }

    const double minorRadius_m = std::max(0.0, static_cast<double>(engine.Config().minorRadius_m));
    const char* profileKind = PlasmaCurrentProfileKindName(engine.PlasmaCurrentProfile().kind);

    for (std::size_t bin = 0; bin < binCount; ++bin) {
        const double normLower = static_cast<double>(bin) / static_cast<double>(binCount);
        const double normUpper = static_cast<double>(bin + 1) / static_cast<double>(binCount);
        const double rInner = minorRadius_m * normLower;
        const double rOuter = minorRadius_m * normUpper;
        const double rCenter = 0.5 * (rInner + rOuter);

        const double meanB = (bin < radialMean.size()) ? radialMean[bin] : 0.0;
        const uint64_t sampleCount = (bin < radialCounts.size()) ? radialCounts[bin] : 0U;

        magneticFieldCsv_ << std::setprecision(std::numeric_limits<double>::max_digits10)
                          << kMilestone7OutputSchemaVersion << ','
                          << telemetry.step << ','
                          << telemetry.time_s << ','
                          << bin << ','
                          << rInner << ','
                          << rOuter << ','
                          << rCenter << ','
                          << meanB << ','
                          << sampleCount << ','
                          << telemetry.magneticField.maxField_T << ','
                          << telemetry.magneticField.recommendedDt_s << ','
                          << profileKind << '\n';
    }

    if (!magneticFieldCsv_.good()) {
        SetError("Failed writing magnetic field diagnostics CSV rows");
        return false;
    }
    return true;
}

bool Milestone7ArtifactExporter::WriteElectrostaticDiagnosticsRow(
    const RunConfig& runConfig,
    const TelemetrySnapshot& telemetry) {
    electrostaticDiagnosticsCsv_ << std::setprecision(std::numeric_limits<double>::max_digits10)
                                 << kMilestone7OutputSchemaVersion << ','
                                 << telemetry.step << ','
                                 << telemetry.time_s << ','
                                 << ElectricFieldModeName(runConfig.electricFieldMode) << ','
                                 << ElectrostaticBoundaryConditionName(runConfig.electrostaticBoundaryCondition) << ','
                                 << ChargeAssignmentSchemeName(runConfig.chargeAssignmentScheme) << ','
                                 << telemetry.electrostaticField.maxElectricField_VPerM << ','
                                 << telemetry.electrostaticField.meanElectricField_VPerM << ','
                                 << telemetry.electrostaticField.solveIterations << ','
                                 << (telemetry.electrostaticField.solveConverged ? "true" : "false") << ','
                                 << telemetry.solverResidual.residualL2 << '\n';

    if (!electrostaticDiagnosticsCsv_.good()) {
        SetError("Failed writing electrostatic diagnostics sidecar row");
        return false;
    }
    return true;
}

bool Milestone7ArtifactExporter::WriteFusionReactivityDiagnosticsRow(
    const RunConfig& runConfig,
    const TelemetrySnapshot& telemetry) {
    double avgSigma_m2 = std::numeric_limits<double>::quiet_NaN();
    double avgProbability = std::numeric_limits<double>::quiet_NaN();
    double avgRelativeSpeed_mPerS = std::numeric_limits<double>::quiet_NaN();
    if (telemetry.stepCounters.fusionKineticsSamples > 0) {
        const double invSamples = 1.0 / static_cast<double>(telemetry.stepCounters.fusionKineticsSamples);
        avgSigma_m2 = telemetry.stepCounters.fusionSigmaSum_m2 * invSamples;
        avgProbability = telemetry.stepCounters.fusionProbabilitySum * invSamples;
        avgRelativeSpeed_mPerS = telemetry.stepCounters.fusionRelativeSpeedSum_mPerS * invSamples;
    }

    fusionReactivityDiagnosticsCsv_ << std::setprecision(std::numeric_limits<double>::max_digits10)
                                    << kMilestone7OutputSchemaVersion << ','
                                    << telemetry.step << ','
                                    << telemetry.time_s << ','
                                    << FusionReactivityModelKindName(runConfig.fusionReactivityModelKind) << ','
                                    << runConfig.fusionCrossSectionScale << ','
                                    << runConfig.fusionProbabilityClamp << ','
                                    << runConfig.fusionMinEnergy_keV << ','
                                    << telemetry.stepCounters.fusionAttempts << ','
                                    << telemetry.stepCounters.fusionAccepted << ','
                                    << telemetry.stepCounters.fusionWeightAttempted << ','
                                    << telemetry.stepCounters.fusionWeightAccepted << ','
                                    << telemetry.stepCounters.fuelWeightConsumedD << ','
                                    << telemetry.stepCounters.fuelWeightConsumedT << ','
                                    << telemetry.stepCounters.ashWeightProducedHe << ','
                                    << avgSigma_m2 << ','
                                    << avgProbability << ','
                                    << avgRelativeSpeed_mPerS << ','
                                    << telemetry.counters.fusionAttempts << ','
                                    << telemetry.counters.fusionAccepted << '\n';

    if (!fusionReactivityDiagnosticsCsv_.good()) {
        SetError("Failed writing fusion reactivity diagnostics sidecar row");
        return false;
    }
    return true;
}

bool Milestone7ArtifactExporter::WriteWallInteractionBridgeRow(
    const RunConfig& runConfig,
    const TelemetrySnapshot& telemetry) {
    wallInteractionBridgeCsv_ << std::setprecision(std::numeric_limits<double>::max_digits10)
                              << kMilestone7OutputSchemaVersion << ','
                              << telemetry.step << ','
                              << telemetry.time_s << ','
                              << WallBoundaryModeName(runConfig.wallBoundaryMode) << ','
                              << runConfig.recycleFraction << ','
                              << telemetry.stepCounters.wallHitCount << ','
                              << telemetry.stepCounters.wallImpactEnergy_J << ','
                              << telemetry.stepCounters.wallLossWeight << ','
                              << telemetry.counters.wallHitCount << ','
                              << telemetry.counters.wallImpactEnergy_J << ','
                              << telemetry.counters.wallLossWeight << '\n';

    if (!wallInteractionBridgeCsv_.good()) {
        SetError("Failed writing wall interaction bridge sidecar row");
        return false;
    }
    return true;
}

bool Milestone7ArtifactExporter::WriteSpeedHistogramRows(
    const TokamakEngine& engine,
    const TelemetrySnapshot& telemetry) {
    const ParticleSystem& particles = engine.Particles();
    const auto& velocities = particles.Velocities();
    const auto& species = particles.Species();

    const std::size_t binCount = config_.speedHistogramBinCount;
    std::vector<uint64_t> counts(binCount, 0);

    double maxSpeed = 0.0;
    uint64_t totalSamples = 0;
    for (std::size_t i = 0; i < particles.Size(); ++i) {
        if (species[i] == ParticleType::Dead) {
            continue;
        }
        const double speed = static_cast<double>(velocities[i].Magnitude());
        maxSpeed = std::max(maxSpeed, speed);
        ++totalSamples;
    }
    maxSpeed = std::max(maxSpeed, 1.0);
    const double width = maxSpeed / static_cast<double>(binCount);

    for (std::size_t i = 0; i < particles.Size(); ++i) {
        if (species[i] == ParticleType::Dead) {
            continue;
        }
        const double speed = static_cast<double>(velocities[i].Magnitude());
        std::size_t bin = static_cast<std::size_t>(speed / width);
        if (bin >= binCount) {
            bin = binCount - 1;
        }
        ++counts[bin];
    }

    for (std::size_t bin = 0; bin < binCount; ++bin) {
        const double lower = width * static_cast<double>(bin);
        const double upper = width * static_cast<double>(bin + 1);
        speedHistogramCsv_ << std::setprecision(std::numeric_limits<double>::max_digits10)
                           << kMilestone7OutputSchemaVersion << ','
                           << telemetry.step << ','
                           << telemetry.time_s << ','
                           << bin << ','
                           << lower << ','
                           << upper << ','
                           << counts[bin] << ','
                           << totalSamples << '\n';
    }

    if (!speedHistogramCsv_.good()) {
        SetError("Failed writing speed histogram CSV rows");
        return false;
    }
    return true;
}

bool Milestone7ArtifactExporter::WritePitchHistogramRows(
    const TokamakEngine& engine,
    const TelemetrySnapshot& telemetry) {
    const ParticleSystem& particles = engine.Particles();
    const TokamakConfig& config = engine.Config();
    const PlasmaCurrentProfileConfig& profileConfig = engine.PlasmaCurrentProfile();
    const auto& positions = particles.Positions();
    const auto& velocities = particles.Velocities();
    const auto& species = particles.Species();

    const std::size_t binCount = config_.pitchHistogramBinCount;
    std::vector<uint64_t> counts(binCount, 0);
    uint64_t totalSamples = 0;
    uint64_t invalidSamples = 0;

    const double width = 180.0 / static_cast<double>(binCount);
    for (std::size_t i = 0; i < particles.Size(); ++i) {
        if (species[i] == ParticleType::Dead) {
            continue;
        }
        ++totalSamples;

        const Vec3 bField = ApproximateBField(config, profileConfig, positions[i]);
        const double vMag = static_cast<double>(velocities[i].Magnitude());
        const double bMag = static_cast<double>(bField.Magnitude());
        if (vMag <= 1.0e-12 || bMag <= 1.0e-12) {
            ++invalidSamples;
            continue;
        }

        double cosTheta =
            static_cast<double>(Vec3::Dot(velocities[i], bField)) /
            (vMag * bMag);
        cosTheta = std::max(-1.0, std::min(1.0, cosTheta));
        const double pitchDeg = std::acos(cosTheta) * (180.0 / kPi);
        std::size_t bin = static_cast<std::size_t>(pitchDeg / width);
        if (bin >= binCount) {
            bin = binCount - 1;
        }
        ++counts[bin];
    }

    for (std::size_t bin = 0; bin < binCount; ++bin) {
        const double lower = width * static_cast<double>(bin);
        const double upper = width * static_cast<double>(bin + 1);
        pitchHistogramCsv_ << std::setprecision(std::numeric_limits<double>::max_digits10)
                           << kMilestone7OutputSchemaVersion << ','
                           << telemetry.step << ','
                           << telemetry.time_s << ','
                           << bin << ','
                           << lower << ','
                           << upper << ','
                           << counts[bin] << ','
                           << totalSamples << ','
                           << invalidSamples << '\n';
    }

    if (!pitchHistogramCsv_.good()) {
        SetError("Failed writing pitch-angle histogram CSV rows");
        return false;
    }
    return true;
}

bool Milestone7ArtifactExporter::WriteSolverResidualRow(const TelemetrySnapshot& telemetry) {
    const SolverResidualSnapshot& residual = telemetry.solverResidual;

    solverResidualCsv_ << std::setprecision(std::numeric_limits<double>::max_digits10)
                       << kMilestone7OutputSchemaVersion << ','
                       << telemetry.step << ','
                       << telemetry.time_s << ','
                       << (residual.residualAvailable ? "true" : "false") << ','
                       << residual.residualL2 << ','
                       << ResidualSolverName(residual.solverKind) << ','
                       << ResidualStatusName(residual.status) << ','
                       << residual.iterations << ','
                       << (residual.converged ? "true" : "false") << ','
                       << residual.tolerance << ','
                       << "\"" << ResidualStatusNote(telemetry) << "\"\n";
    if (!solverResidualCsv_.good()) {
        SetError("Failed writing solver residual contract row");
        return false;
    }
    return true;
}

bool Milestone7ArtifactExporter::WriteParticleSnapshotCsv(
    const TokamakEngine& engine,
    const TelemetrySnapshot& telemetry) {
    const ParticleSystem& particles = engine.Particles();
    const auto& positions = particles.Positions();
    const auto& velocities = particles.Velocities();
    const auto& masses = particles.Masses();
    const auto& charges = particles.Charges();
    const auto& qOverM = particles.ChargeToMass();
    const auto& species = particles.Species();

    const std::size_t totalParticles = particles.Size();
    const std::size_t target = std::min(totalParticles, config_.maxParticlesPerSnapshot);
    const std::size_t sampleStride = (target == 0)
        ? 1
        : std::max<std::size_t>(1, (totalParticles + target - 1) / target);

    std::ostringstream filename;
    filename << "particles_step_" << std::setfill('0') << std::setw(8) << telemetry.step << ".csv";
    const std::string relativePath = std::string("snapshots/") + filename.str();
    const std::filesystem::path path = std::filesystem::path(runDirectory_) / relativePath;

    std::ofstream out(path.string(), std::ios::trunc);
    if (!out.is_open()) {
        SetError("Failed to open particle snapshot CSV path: " + path.string());
        return false;
    }

    out << "schema_version,step,time_s,total_particles,sampled_particles,sample_stride,particle_index,species,species_name,"
        << "x_m,y_m,z_m,vx_m_per_s,vy_m_per_s,vz_m_per_s,mass_kg,charge_c,q_over_m\n";

    std::size_t sampledCount = 0;
    for (std::size_t i = 0; i < totalParticles; i += sampleStride) {
        if (sampledCount >= target) {
            break;
        }
        out << std::setprecision(std::numeric_limits<double>::max_digits10)
            << kMilestone7OutputSchemaVersion << ','
            << telemetry.step << ','
            << telemetry.time_s << ','
            << totalParticles << ','
            << target << ','
            << sampleStride << ','
            << i << ','
            << static_cast<uint32_t>(species[i]) << ','
            << SpeciesName(species[i]) << ','
            << positions[i].x << ','
            << positions[i].y << ','
            << positions[i].z << ','
            << velocities[i].x << ','
            << velocities[i].y << ','
            << velocities[i].z << ','
            << masses[i] << ','
            << charges[i] << ','
            << qOverM[i] << '\n';
        ++sampledCount;
    }

    out.flush();
    if (!out.good()) {
        SetError("Failed writing particle snapshot CSV path: " + path.string());
        return false;
    }

    particleSnapshotRelativePaths_.push_back(relativePath);
    return true;
}

bool Milestone7ArtifactExporter::WriteManifestJson() {
    const std::filesystem::path manifestPath = std::filesystem::path(runDirectory_) / "manifest_v2.json";
    std::ofstream out(manifestPath.string(), std::ios::trunc);
    if (!out.is_open()) {
        SetError("Failed to open manifest JSON path: " + manifestPath.string());
        return false;
    }

    out << "{\n";
    out << "  \"schema\": \"tokamak.milestone7.manifest\",\n";
    out << "  \"schema_version\": " << kMilestone7OutputSchemaVersion << ",\n";
    out << "  \"run_id\": \"" << JsonEscape(runId_) << "\",\n";
    out << "  \"created_utc\": \"" << TimestampIsoUtc() << "\",\n";
    out << "  \"files\": {\n";
    out << "    \"run_config\": \"" << runConfigRelativePath_ << "\",\n";
    out << "    \"summary_csv\": \"" << summaryRelativePath_ << "\",\n";
    out << "    \"radial_profiles_csv\": \"" << radialRelativePath_ << "\",\n";
    out << "    \"magnetic_field_diagnostics_csv\": \"" << magneticFieldRelativePath_ << "\",\n";
    out << "    \"electrostatic_diagnostics_csv\": \"" << electrostaticDiagnosticsRelativePath_ << "\",\n";
    out << "    \"speed_histogram_csv\": \"" << speedHistogramRelativePath_ << "\",\n";
    out << "    \"pitch_histogram_csv\": \"" << pitchHistogramRelativePath_ << "\",\n";
    out << "    \"solver_residual_log_csv\": \"" << solverResidualRelativePath_ << "\",\n";
    out << "    \"particle_snapshot_csv_files\": [";
    for (std::size_t i = 0; i < particleSnapshotRelativePaths_.size(); ++i) {
        if (i > 0) {
            out << ", ";
        }
        out << "\"" << JsonEscape(particleSnapshotRelativePaths_[i]) << "\"";
    }
    out << "]\n";
    out << "  }\n";
    out << "}\n";

    out.flush();
    if (!out.good()) {
        SetError("Failed to write manifest JSON file: " + manifestPath.string());
        return false;
    }
    return true;
}

void Milestone7ArtifactExporter::CloseFiles() {
    if (summaryCsv_.is_open()) {
        summaryCsv_.close();
    }
    if (radialCsv_.is_open()) {
        radialCsv_.close();
    }
    if (magneticFieldCsv_.is_open()) {
        magneticFieldCsv_.close();
    }
    if (electrostaticDiagnosticsCsv_.is_open()) {
        electrostaticDiagnosticsCsv_.close();
    }
    if (fusionReactivityDiagnosticsCsv_.is_open()) {
        fusionReactivityDiagnosticsCsv_.close();
    }
    if (wallInteractionBridgeCsv_.is_open()) {
        wallInteractionBridgeCsv_.close();
    }
    if (speedHistogramCsv_.is_open()) {
        speedHistogramCsv_.close();
    }
    if (pitchHistogramCsv_.is_open()) {
        pitchHistogramCsv_.close();
    }
    if (solverResidualCsv_.is_open()) {
        solverResidualCsv_.close();
    }
}

void Milestone7ArtifactExporter::SetError(const std::string& message) {
    lastError_ = message;
}

}  // namespace tokamak
