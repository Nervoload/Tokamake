#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

if command -v cmake >/dev/null 2>&1; then
  BUILD_DIR="${BUILD_DIR:-${ROOT_DIR}/build/debug}"
  cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON
  cmake --build "${BUILD_DIR}" --config Debug --target tokamakfusion tokamak_tests
  ctest --test-dir "${BUILD_DIR}" --output-on-failure
  exit 0
fi

echo "cmake not found; using manual fallback build in build_manual/" >&2
mkdir -p build_manual

CXX_BIN="${CXX:-clang++}"

"${CXX_BIN}" -std=c++17 -Iinclude \
  src/artifact_export.cpp \
  src/collision.cpp \
  src/electrostatic_field.cpp \
  src/electrostatic_models.cpp \
  src/engine.cpp \
  src/magnetic_field.cpp \
  src/particle_push.cpp \
  src/particle_system.cpp \
  src/spatial_grid.cpp \
  src/telemetry.cpp \
  src/main.cpp \
  -o build_manual/tokamakfusion

"${CXX_BIN}" -std=c++17 -Iinclude -Ithird_party/googletest -Ithird_party/googletest/include \
  -DTOKAMAKFUSION_PATH=\"./build_manual/tokamakfusion\" \
  src/artifact_export.cpp \
  src/collision.cpp \
  src/electrostatic_field.cpp \
  src/electrostatic_models.cpp \
  src/engine.cpp \
  src/magnetic_field.cpp \
  src/particle_push.cpp \
  src/particle_system.cpp \
  src/spatial_grid.cpp \
  src/telemetry.cpp \
  tests/engine_tests.cpp \
  tests/integration_tests.cpp \
  tests/milestone2_magnetic_field_tests.cpp \
  tests/milestone3_electrostatic_field_tests.cpp \
  tests/milestone6_validation_tests.cpp \
  tests/particle_system_tests.cpp \
  third_party/googletest/src/gtest.cc \
  third_party/googletest/src/gtest_main.cc \
  -pthread \
  -o build_manual/tokamak_tests

./build_manual/tokamak_tests
