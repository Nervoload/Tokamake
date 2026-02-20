#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${BUILD_DIR:-${repo_root}/build/benchmark}"
results_dir="${RESULTS_DIR:-${repo_root}/benchmarks/results}"

scenario="${SCENARIO:-ignition}"
seed="${SEED:-424242}"
steps="${STEPS:-5000}"
dt="${DT:-1e-7}"
telemetry_every="${TELEMETRY_EVERY:-500}"
particle_cap="${PARTICLE_CAP:-2500000}"

mkdir -p "${results_dir}"
timestamp="$(date +%Y%m%d_%H%M%S)"
sim_log="${results_dir}/benchmark_${timestamp}.log"
time_log="${results_dir}/benchmark_${timestamp}.time"

binary_path="${build_dir}/tokamakfusion"
if command -v cmake >/dev/null 2>&1; then
    cmake -S "${repo_root}" -B "${build_dir}" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
    cmake --build "${build_dir}" --config Release --target tokamakfusion
elif [[ -x "${repo_root}/build_manual/tokamakfusion" ]]; then
    binary_path="${repo_root}/build_manual/tokamakfusion"
    echo "cmake not found; using prebuilt binary at ${binary_path}"
else
    echo "Error: cmake not found and no prebuilt binary at build_manual/tokamakfusion" >&2
    exit 1
fi

command=(
    "${binary_path}"
    "--scenario" "${scenario}"
    "--seed" "${seed}"
    "--dt" "${dt}"
    "--steps" "${steps}"
    "--telemetry-every" "${telemetry_every}"
    "--particle-cap" "${particle_cap}"
)

/usr/bin/time -p -o "${time_log}" "${command[@]}" > "${sim_log}"

real_seconds="$(awk '$1=="real" {print $2}' "${time_log}")"
steps_per_second="$(awk -v steps="${steps}" -v seconds="${real_seconds}" 'BEGIN { if (seconds > 0) printf "%.2f", steps / seconds; else printf "inf"; }')"

echo "Benchmark complete"
echo "  command: ${command[*]}"
echo "  wall_time_s: ${real_seconds}"
echo "  steps_per_s: ${steps_per_second}"
echo "  simulation_log: ${sim_log}"
echo "  timing_log: ${time_log}"
