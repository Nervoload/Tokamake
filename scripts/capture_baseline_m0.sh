#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${BUILD_DIR:-${repo_root}/build}"
artifact_root="${repo_root}/baselines/m0/generated_runs"
reference_dir="${repo_root}/baselines/m0/reference"
raw_log="${artifact_root}/console.raw.log"

normalize_console() {
  local in_path="$1"
  local out_path="$2"
  sed -E \
    -e 's/id=run_[^[:space:]]+/id=<RUN_ID>/g' \
    -e 's#dir=[^[:space:]]+#dir=<ARTIFACT_DIR>#g' \
    -e 's#ARTIFACTS WRITTEN \| .*#ARTIFACTS WRITTEN | <ARTIFACT_DIR>#g' \
    "$in_path" > "$out_path"
}

cmake -S "${repo_root}" -B "${build_dir}" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build "${build_dir}" --config Release --target tokamakfusion

rm -rf "${artifact_root}"
mkdir -p "${artifact_root}" "${reference_dir}"

binary_path="${build_dir}/tokamakfusion"
"${binary_path}" \
  --scenario cold \
  --seed 20260220 \
  --dt 1e-7 \
  --steps 200 \
  --telemetry-every 50 \
  --artifacts-root "${artifact_root}" > "${raw_log}"

run_dir="$(grep -m1 'ARTIFACT RUN | id=' "${raw_log}" | sed -E 's/.* dir=([^[:space:]]+).*/\1/')"
if [[ -z "${run_dir}" || ! -d "${run_dir}" ]]; then
  echo "Failed to resolve artifact run directory from baseline capture output" >&2
  exit 1
fi

normalize_console "${raw_log}" "${reference_dir}/console.log"
cp "${run_dir}/summary_v2.csv" "${reference_dir}/summary_v2.csv"
cp "${run_dir}/radial_profiles_v2.csv" "${reference_dir}/radial_profiles_v2.csv"
cp "${run_dir}/magnetic_field_diagnostics_v2.csv" "${reference_dir}/magnetic_field_diagnostics_v2.csv"
cp "${run_dir}/solver_residuals_v2.csv" "${reference_dir}/solver_residuals_v2.csv"

(
  cd "${reference_dir}"
  shasum -a 256 \
    console.log \
    summary_v2.csv \
    radial_profiles_v2.csv \
    magnetic_field_diagnostics_v2.csv \
    solver_residuals_v2.csv > checksums.sha256
)

echo "Milestone 0 baseline captured at ${reference_dir}"
