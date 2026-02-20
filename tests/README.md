# Tests and Validation

This directory contains the CMake/CTest validation suite for the modular tokamak core.

## Entry point

Run validation from the repository root:

```bash
./tests/run_validation.sh
```

This script configures a Debug build in `build/debug`, builds `tokamakfusion` + `tokamak_tests`, and executes `ctest --output-on-failure`.
If `cmake` is not installed, it compiles a manual fallback build in `build_manual/`.

## Scope of current tests

- `particle_system_tests.cpp`: capacity and array-consistency checks.
- `engine_tests.cpp`: deterministic seeded behavior, insertion safety, long-run finite-state checks.
- `integration_tests.cpp`: collision behavior and telemetry/run-manifest smoke coverage.

## Deterministic stochastic tests

1. Use explicit seeds in tests (`RunConfig.seed` or local RNG seed values).
2. Keep seed choices stable and committed.
3. When adding a stochastic test, include a short comment describing why the selected seed set is representative.

## Related workflow docs

- `CONTRIBUTING.md` for required pre-merge checks.
- `AGENT.md` for agent implementation workflow.
- `PHYSICS_MANIFEST.md` for assumptions tied to validation expectations.
