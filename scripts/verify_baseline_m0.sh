#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${BUILD_DIR:-${repo_root}/build}"
reference_dir="${repo_root}/baselines/m0/reference"

normalize_console() {
  local in_path="$1"
  local out_path="$2"
  sed -E \
    -e 's/id=run_[^[:space:]]+/id=<RUN_ID>/g' \
    -e 's#dir=[^[:space:]]+#dir=<ARTIFACT_DIR>#g' \
    -e 's#ARTIFACTS WRITTEN \| .*#ARTIFACTS WRITTEN | <ARTIFACT_DIR>#g' \
    "$in_path" > "$out_path"
}

for required in \
  "${reference_dir}/console.log" \
  "${reference_dir}/summary_v2.csv" \
  "${reference_dir}/radial_profiles_v2.csv" \
  "${reference_dir}/magnetic_field_diagnostics_v2.csv" \
  "${reference_dir}/solver_residuals_v2.csv" \
  "${reference_dir}/checksums.sha256"; do
  if [[ ! -f "${required}" ]]; then
    echo "Missing baseline reference file: ${required}" >&2
    exit 1
  fi
done

cmake -S "${repo_root}" -B "${build_dir}" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build "${build_dir}" --config Release --target tokamakfusion

work_dir="$(mktemp -d /tmp/tokamak_m0_verify.XXXXXX)"
trap 'rm -rf "${work_dir}"' EXIT

artifact_root="${work_dir}/generated_runs"
raw_log="${work_dir}/console.raw.log"
normalized_log="${work_dir}/console.log"
mkdir -p "${artifact_root}"

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
  echo "Failed to resolve artifact run directory from verification run output" >&2
  exit 1
fi

normalize_console "${raw_log}" "${normalized_log}"
cp "${run_dir}/summary_v2.csv" "${work_dir}/summary_v2.csv"
cp "${run_dir}/radial_profiles_v2.csv" "${work_dir}/radial_profiles_v2.csv"
cp "${run_dir}/magnetic_field_diagnostics_v2.csv" "${work_dir}/magnetic_field_diagnostics_v2.csv"
cp "${run_dir}/solver_residuals_v2.csv" "${work_dir}/solver_residuals_v2.csv"

for file in \
  console.log \
  summary_v2.csv \
  radial_profiles_v2.csv \
  magnetic_field_diagnostics_v2.csv \
  solver_residuals_v2.csv; do
  diff -u "${reference_dir}/${file}" "${work_dir}/${file}"
done

(
  cd "${work_dir}"
  shasum -a 256 \
    console.log \
    summary_v2.csv \
    radial_profiles_v2.csv \
    magnetic_field_diagnostics_v2.csv \
    solver_residuals_v2.csv > checksums.sha256
)

diff -u "${reference_dir}/checksums.sha256" "${work_dir}/checksums.sha256"

echo "Milestone 0 baseline verification passed"
