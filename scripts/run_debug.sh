#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${BUILD_DIR:-${repo_root}/build/debug}"

scenario="${SCENARIO:-cold}"
seed="${SEED:-20260220}"
steps="${STEPS:-200}"
dt="${DT:-1e-7}"
telemetry_every="${TELEMETRY_EVERY:-50}"
particle_cap="${PARTICLE_CAP:-2500000}"

binary_path="${build_dir}/tokamakfusion"
if command -v cmake >/dev/null 2>&1; then
    cmake -S "${repo_root}" -B "${build_dir}" -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON
    cmake --build "${build_dir}" --config Debug
elif [[ -x "${repo_root}/build_manual/tokamakfusion" ]]; then
    binary_path="${repo_root}/build_manual/tokamakfusion"
    echo "cmake not found; using prebuilt binary at ${binary_path}"
else
    echo "Error: cmake not found and no prebuilt binary at build_manual/tokamakfusion" >&2
    exit 1
fi

exec "${binary_path}" \
    --scenario "${scenario}" \
    --seed "${seed}" \
    --dt "${dt}" \
    --steps "${steps}" \
    --telemetry-every "${telemetry_every}" \
    --particle-cap "${particle_cap}" \
    "$@"
