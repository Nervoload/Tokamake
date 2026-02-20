# Benchmarks

This directory provides a deterministic benchmark entry point for the CMake tokamak executable.

## Entry point

Run from the repository root:

```bash
./benchmarks/run_benchmark.sh
```

The script will:
1. configure/build a Release binary in `build/benchmark`,
2. run a fixed-seed simulation profile, and
3. write logs in `benchmarks/results/` with wall time and computed steps/second.

If `cmake` is not installed, the script falls back to `build_manual/tokamakfusion` when available.

## Default benchmark profile

- Scenario: `ignition`
- Seed: `424242`
- Steps: `5000`
- Time step: `1e-7`
- Telemetry cadence: `500`

## Overrides

You can override settings via environment variables:

```bash
SEED=7 STEPS=20000 TELEMETRY_EVERY=1000 ./benchmarks/run_benchmark.sh
```

Supported variables: `BUILD_DIR`, `RESULTS_DIR`, `SCENARIO`, `SEED`, `STEPS`, `DT`, `TELEMETRY_EVERY`, `PARTICLE_CAP`.

## Use with contribution workflow

When touching hot paths, include benchmark output (wall time + steps/s) in your change notes alongside validation output from `./tests/run_validation.sh`.
