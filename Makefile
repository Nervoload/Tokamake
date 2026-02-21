CMAKE ?= cmake

BUILD_DIR ?= build/release
TEST_BUILD_DIR ?= build/debug
BENCH_BUILD_DIR ?= build/benchmark

.PHONY: all sim test benchmark debug clean viewer

all: sim

sim:
	$(CMAKE) -S . -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
	$(CMAKE) --build $(BUILD_DIR) --target tokamakfusion

test:
	./tests/run_validation.sh

benchmark:
	./benchmarks/run_benchmark.sh

debug:
	./scripts/run_debug.sh

viewer:
	@echo "Build the canonical replay viewer with CMake:"
	@echo "  cmake -S . -B build -DBUILD_VIEWER=ON"
	@echo "  cmake --build build --target tokamak_viewer"
	@echo "See docs/VIEWER_SETUP.md for dependency setup."

clean:
	rm -rf $(BUILD_DIR) $(TEST_BUILD_DIR) $(BENCH_BUILD_DIR) build_manual tokamakfusion tokamakfusion_manifest
