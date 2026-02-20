# Milestone 6 Validation Matrix

## Run command

```bash
./tests/run_validation.sh
```

Expected high-level output:
- GoogleTest executes all `Milestone6Validation*` tests.
- Final line reports: `[  PASSED  ] All tests passed`.

## Validation status

| Test | Status | Tolerance | Rationale |
|---|---|---|---|
| Uniform B-field gyrofrequency and gyroradius | Pass | `1%` freq, `1.5%` radius (relative) | Captures Boris finite-`dt` phase error while still detecting integration regressions. |
| E×B drift velocity | Pass | `2%` relative + `5 m/s` absolute | Averaged drift converges with finite-step gyro motion; this bound is tight enough to catch sign/magnitude mistakes. |
| Wall reflection speed conservation | Pass | `1e-6` relative + `1e-3 m/s` absolute | Reflection should be elastic; tolerance only absorbs float roundoff. |
| Collision reproducibility with fixed seed | Pass | Exact event/counter equality | Deterministic RNG path should be bitwise repeatable for identical setup and seed. |
| Grid consistency/convergence | Pass | Exact rebinned count equality | Same particles sorted on 2x finer grid should aggregate back to the coarse histogram exactly. |
| Electrostatic Poisson known-charge distribution | Pass | Potential L2 `<0.3%`, E-field L2 `<1.5%`, residual `<=1e-6` | 1D Dirichlet Poisson solve with sinusoidal charge density is compared against analytic potential/field. |

## Determinism notes

- All stochastic tests use explicit `std::mt19937` seeds.
- No test depends on wall-clock time or non-deterministic seeding.
- Shared helpers live in `tests/validation_harness.hpp` to keep seed/tolerance policy centralized.
- Current Poisson test validates deterministic electrostatic numerics with an analytic benchmark while core runtime solver APIs are still evolving.
