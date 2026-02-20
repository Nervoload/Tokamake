# Contributing to Tokamak

This repository uses a single canonical path:
- CMake modular runtime (`src/main.cpp`, `include/tokamak/*`, `tests/*`).

Deprecated monolith files are archived under `legacy/` and are not part of contribution checks.
Legacy viewer/snapshot tooling is also archived under `legacy/` and excluded from canonical targets.

## Read before editing

1. `AGENT.md` for architecture, guardrails, and patch workflow.
2. `PHYSICS_MANIFEST.md` for modeled vs placeholder assumptions.
3. `MILESTONE_CHECKLIST.md` for acceptance criteria and milestone status.

## Standard entry points

From the repository root:

1. Debug run (deterministic smoke run):
```bash
./scripts/run_debug.sh
```

2. Validation run (tests):
```bash
./tests/run_validation.sh
```

3. Benchmark run (deterministic release profile):
```bash
./benchmarks/run_benchmark.sh
```

4. Milestone 0 baseline capture:
```bash
./scripts/capture_baseline_m0.sh
```

5. Milestone 0 baseline verification:
```bash
./scripts/verify_baseline_m0.sh
```

## Magnetic profile CLI examples

Examples from repository root:

1. Uniform profile (default behavior, explicit):
```bash
./scripts/run_debug.sh --current-profile uniform
```

2. Parabolic profile with custom diagnostics bins:
```bash
./scripts/run_debug.sh --current-profile parabolic --mag-field-bins 64 --mag-field-dt-safety 0.04
```

3. Custom profile table:
```bash
./scripts/run_debug.sh --current-profile custom --current-profile-table ./path/to/profile.csv
```

If `cmake` is unavailable, fallback behavior is script-specific:
- `tests/run_validation.sh` performs a manual fallback compile in `build_manual/`.
- `scripts/run_debug.sh` and `benchmarks/run_benchmark.sh` use prebuilt fallback binaries in `build_manual/` when present.

## Minimal checks required before merge

1. `./scripts/run_debug.sh` succeeds with a fixed seed.
2. `./tests/run_validation.sh` passes.
3. If performance-sensitive code changed, run `./benchmarks/run_benchmark.sh` and include the result in your change notes.
4. If physics assumptions changed, update `PHYSICS_MANIFEST.md` and relevant acceptance notes in `MILESTONE_CHECKLIST.md`.
5. If baseline artifacts are intentionally updated, create or move the annotated baseline tag using:
   - `baseline/m0-YYYY-MM-DD`

## Deterministic seed policy

1. Any stochastic test must use explicit fixed seeds (`RunConfig.seed` or CLI `--seed`).
2. Do not rely on `std::random_device` behavior in tests or reproducibility checks.
3. Multi-seed tests must use a fixed seed list/range committed to the repository.
4. Include the seed(s) used in validation notes when behavior is stochastic.

## Updating physics assumptions and docs

Update `PHYSICS_MANIFEST.md` when any change affects:
- what is physically modeled,
- placeholder approximations,
- intentionally non-physical stabilizations,
- required diagnostics/validation for that component.

When this happens, also:
1. add or update at least one validation or telemetry check proving the new behavior, and
2. update milestone criteria/status notes in `MILESTONE_CHECKLIST.md` if acceptance conditions changed.
