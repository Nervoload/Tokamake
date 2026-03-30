#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${BUILD_DIR:-${repo_root}/build/viewer-debug}"
artifact_root="${ARTIFACT_ROOT:-${repo_root}/output/runs}"

profile="${PROFILE:-fusion}"
scenario="${SCENARIO:-}"
seed="${SEED:-20260220}"
dt="${DT:-1e-7}"
sim_duration_ms="${SIM_DURATION_MS:-}"
sim_duration_s="${SIM_DURATION_S:-}"
steps_override="${STEPS:-}"
snapshots_per_ms="${SNAPSHOTS_PER_MS:-1000}"
telemetry_per_ms="${TELEMETRY_PER_MS:-100}"
particle_snapshot_max="${PARTICLE_SNAPSHOT_MAX:-6000}"
particle_cap="${PARTICLE_CAP:-6000}"
fusion_cross_section_scale="${FUSION_CROSS_SECTION_SCALE:-}"
fusion_probability_clamp="${FUSION_PROBABILITY_CLAMP:-}"
fusion_min_energy_kev="${FUSION_MIN_ENERGY_KEV:-}"
startup_ramp_ms="${STARTUP_RAMP_MS:-}"
fusion_start_delay_ms="${FUSION_START_DELAY_MS:-}"
fusion_ramp_ms="${FUSION_RAMP_MS:-}"

point_size="${POINT_SIZE:-3.5}"
playback_rate="${PLAYBACK_RATE:-1.0}"
start_step="${START_STEP:-}"
launch_viewer="${LAUNCH_VIEWER:-1}"

round_to_int() {
    awk -v value="$1" 'BEGIN {
        if (value < 1.0) {
            value = 1.0;
        }
        printf "%.0f", value;
    }'
}

ceil_div_int() {
    awk -v numerator="$1" -v denominator="$2" 'BEGIN {
        value = int((numerator + denominator - 1) / denominator);
        if (value < 1) {
            value = 1;
        }
        print value;
    }'
}

print_usage() {
    cat <<EOF
Usage: ./scripts/run_full_system.sh [--no-viewer] [--help]

Builds the simulator + viewer, runs a replay-producing simulation, resolves the
artifact run directory, and launches the viewer on that run.

Environment overrides:
  BUILD_DIR                  default: ${repo_root}/build/viewer-debug
  ARTIFACT_ROOT              default: ${repo_root}/output/runs
  PROFILE                    default: fusion
  SCENARIO                   default: ignition (for fusion/baseline profiles)
  SEED                       default: 20260220
  DT                         default: 1e-7
  SIM_DURATION_MS            default: 0.12
  SIM_DURATION_S             default: unset
  STEPS                      default: unset
  SNAPSHOTS_PER_MS           default: 1000
  TELEMETRY_PER_MS           default: 100
  TELEMETRY_EVERY            default: derived from TELEMETRY_PER_MS
  ARTIFACT_EVERY             default: derived from SNAPSHOTS_PER_MS
  PARTICLE_SNAPSHOT_EVERY    default: derived from SNAPSHOTS_PER_MS
  PARTICLE_SNAPSHOT_MAX      default: 6000
  PARTICLE_CAP               default: 6000
  FUSION_CROSS_SECTION_SCALE default: 1e12 for PROFILE=fusion, 1 for PROFILE=baseline
  FUSION_PROBABILITY_CLAMP   default: engine default
  FUSION_MIN_ENERGY_KEV      default: engine default
  STARTUP_RAMP_MS            default: 0.04 for fusion/baseline
  FUSION_START_DELAY_MS      default: 0.015 for fusion/baseline
  FUSION_RAMP_MS             default: 0.02 for fusion/baseline
  POINT_SIZE                 default: 3.5
  PLAYBACK_RATE              default: 1.0
  START_STEP                 default: unset

Examples:
  ./scripts/run_full_system.sh
  PROFILE=baseline SIM_DURATION_MS=0.24 ./scripts/run_full_system.sh
  STARTUP_RAMP_MS=0.08 FUSION_START_DELAY_MS=0.03 FUSION_RAMP_MS=0.05 ./scripts/run_full_system.sh
  SIM_DURATION_MS=0.50 SNAPSHOTS_PER_MS=600 ./scripts/run_full_system.sh
  SCENARIO=failure LAUNCH_VIEWER=0 ./scripts/run_full_system.sh --no-viewer
EOF
}

while (($# > 0)); do
    case "$1" in
        --help|-h)
            print_usage
            exit 0
            ;;
        --no-viewer)
            launch_viewer=0
            shift
            ;;
        *)
            echo "Unknown option: $1" >&2
            print_usage >&2
            exit 1
            ;;
    esac
done

if ! command -v cmake >/dev/null 2>&1; then
    echo "Error: cmake is required to build the simulator and viewer." >&2
    exit 1
fi

case "${profile}" in
    fusion)
        if [[ -z "${scenario}" ]]; then
            scenario="ignition"
        fi
        if [[ -z "${fusion_cross_section_scale}" ]]; then
            fusion_cross_section_scale="1000000000000"
        fi
        if [[ -z "${startup_ramp_ms}" ]]; then
            startup_ramp_ms="0.04"
        fi
        if [[ -z "${fusion_start_delay_ms}" ]]; then
            fusion_start_delay_ms="0.015"
        fi
        if [[ -z "${fusion_ramp_ms}" ]]; then
            fusion_ramp_ms="0.02"
        fi
        default_duration_ms="0.18"
        ;;
    baseline)
        if [[ -z "${scenario}" ]]; then
            scenario="ignition"
        fi
        if [[ -z "${fusion_cross_section_scale}" ]]; then
            fusion_cross_section_scale="1"
        fi
        if [[ -z "${startup_ramp_ms}" ]]; then
            startup_ramp_ms="0.04"
        fi
        if [[ -z "${fusion_start_delay_ms}" ]]; then
            fusion_start_delay_ms="0.015"
        fi
        if [[ -z "${fusion_ramp_ms}" ]]; then
            fusion_ramp_ms="0.02"
        fi
        default_duration_ms="0.18"
        ;;
    *)
        echo "Error: unsupported PROFILE=${profile}. Supported values: fusion, baseline" >&2
        exit 1
        ;;
esac

if [[ -n "${sim_duration_ms}" && -n "${sim_duration_s}" ]]; then
    echo "Error: specify only one of SIM_DURATION_MS or SIM_DURATION_S." >&2
    exit 1
fi

if [[ -n "${sim_duration_ms}" ]]; then
    sim_duration_s="$(awk -v ms="${sim_duration_ms}" 'BEGIN { printf "%.12g", ms / 1000.0; }')"
elif [[ -n "${sim_duration_s}" ]]; then
    sim_duration_ms="$(awk -v seconds="${sim_duration_s}" 'BEGIN { printf "%.12g", seconds * 1000.0; }')"
elif [[ -n "${steps_override}" ]]; then
    sim_duration_s="$(awk -v steps="${steps_override}" -v dt="${dt}" 'BEGIN { printf "%.12g", steps * dt; }')"
    sim_duration_ms="$(awk -v seconds="${sim_duration_s}" 'BEGIN { printf "%.12g", seconds * 1000.0; }')"
else
    sim_duration_ms="${default_duration_ms}"
    sim_duration_s="$(awk -v ms="${sim_duration_ms}" 'BEGIN { printf "%.12g", ms / 1000.0; }')"
fi

if [[ -n "${steps_override}" && -z "${SIM_DURATION_MS:-}" && -z "${SIM_DURATION_S:-}" ]]; then
    steps="${steps_override}"
else
    steps="$(round_to_int "$(awk -v seconds="${sim_duration_s}" -v dt="${dt}" 'BEGIN { print seconds / dt; }')")"
fi

target_snapshot_count="$(round_to_int "$(awk -v ms="${sim_duration_ms}" -v rate="${snapshots_per_ms}" 'BEGIN { print ms * rate; }')")"
derived_snapshot_every="$(ceil_div_int "${steps}" "${target_snapshot_count}")"
particle_snapshot_every="${PARTICLE_SNAPSHOT_EVERY:-${derived_snapshot_every}}"
artifact_every="${ARTIFACT_EVERY:-${derived_snapshot_every}}"

target_telemetry_count="$(round_to_int "$(awk -v ms="${sim_duration_ms}" -v rate="${telemetry_per_ms}" 'BEGIN { print ms * rate; }')")"
derived_telemetry_every="$(ceil_div_int "${steps}" "${target_telemetry_count}")"
telemetry_every="${TELEMETRY_EVERY:-${derived_telemetry_every}}"

echo "Profile: ${profile}"
echo "Scenario: ${scenario}"
echo "Sim duration: ${sim_duration_ms} ms"
echo "Steps: ${steps}"
echo "Snapshot cadence: every ${particle_snapshot_every} steps"
echo "Fusion cross-section scale: ${fusion_cross_section_scale}"
echo "Startup ramp: ${startup_ramp_ms:-0} ms"
echo "Fusion start delay: ${fusion_start_delay_ms:-0} ms"
echo "Fusion ramp: ${fusion_ramp_ms:-0} ms"

echo "Configuring viewer build in ${build_dir}"
cmake -S "${repo_root}" -B "${build_dir}" -DCMAKE_BUILD_TYPE=Debug -DBUILD_VIEWER=ON -DBUILD_TESTING=ON

echo "Building tokamakfusion + tokamak_viewer"
cmake --build "${build_dir}" --config Debug --target tokamakfusion tokamak_viewer

binary_path="${build_dir}/tokamakfusion"
viewer_path="${build_dir}/tokamak_viewer"

if [[ ! -x "${binary_path}" ]]; then
    echo "Error: simulator binary not found at ${binary_path}" >&2
    exit 1
fi
if [[ ! -x "${viewer_path}" ]]; then
    echo "Error: viewer binary not found at ${viewer_path}" >&2
    exit 1
fi

mkdir -p "${artifact_root}"

log_path="$(mktemp "${TMPDIR:-/tmp}/tokamak_run_full.XXXXXX")"
trap 'rm -f "${log_path}"' EXIT

sim_cmd=(
    "${binary_path}"
    --scenario "${scenario}"
    --seed "${seed}"
    --dt "${dt}"
    --steps "${steps}"
    --telemetry-every "${telemetry_every}"
    --particle-cap "${particle_cap}"
    --fusion-cross-section-scale "${fusion_cross_section_scale}"
    --artifact-every "${artifact_every}"
    --particle-snapshot-every "${particle_snapshot_every}"
    --particle-snapshot-max "${particle_snapshot_max}"
    --artifacts-root "${artifact_root}"
)

if [[ -n "${fusion_probability_clamp}" ]]; then
    sim_cmd+=(--fusion-probability-clamp "${fusion_probability_clamp}")
fi

if [[ -n "${fusion_min_energy_kev}" ]]; then
    sim_cmd+=(--fusion-min-energy-kev "${fusion_min_energy_kev}")
fi

if [[ -n "${startup_ramp_ms}" ]]; then
    sim_cmd+=(--startup-ramp-ms "${startup_ramp_ms}")
fi

if [[ -n "${fusion_start_delay_ms}" ]]; then
    sim_cmd+=(--fusion-start-delay-ms "${fusion_start_delay_ms}")
fi

if [[ -n "${fusion_ramp_ms}" ]]; then
    sim_cmd+=(--fusion-ramp-ms "${fusion_ramp_ms}")
fi

echo "Running simulation"
"${sim_cmd[@]}" | tee "${log_path}"

run_dir="$(sed -n 's/^ARTIFACTS WRITTEN | //p' "${log_path}" | tail -n1)"
if [[ -z "${run_dir}" ]]; then
    run_dir="$(sed -nE 's/^ARTIFACT RUN \| id=[^ ]+ dir=(.*)$/\1/p' "${log_path}" | tail -n1)"
fi

if [[ -z "${run_dir}" || ! -d "${run_dir}" ]]; then
    echo "Error: failed to resolve artifact run directory from simulator output." >&2
    exit 1
fi

echo "Run directory: ${run_dir}"

if [[ "${launch_viewer}" != "1" ]]; then
    exit 0
fi

viewer_cmd=(
    "${viewer_path}"
    --run-dir "${run_dir}"
    --point-size "${point_size}"
    --playback-rate "${playback_rate}"
)

if [[ -n "${start_step}" ]]; then
    viewer_cmd+=(--start-step "${start_step}")
fi

echo "Launching viewer"
exec "${viewer_cmd[@]}"
