# Milestone 0 Baseline Specification

This baseline freezes a deterministic reference run for regression checks.

## Canonical command

```bash
./build/tokamakfusion \
  --scenario cold \
  --seed 20260220 \
  --dt 1e-7 \
  --steps 200 \
  --telemetry-every 50 \
  --artifacts-root ./baselines/m0/generated_runs
```

## Reference files

- `baselines/m0/reference/console.log`
- `baselines/m0/reference/summary_v2.csv`
- `baselines/m0/reference/radial_profiles_v2.csv`
- `baselines/m0/reference/magnetic_field_diagnostics_v2.csv`
- `baselines/m0/reference/solver_residuals_v2.csv`
- `baselines/m0/reference/checksums.sha256`

## Workflow

1. Refresh baseline reference artifacts:
   - `./scripts/capture_baseline_m0.sh`
2. Verify current code reproduces baseline:
   - `./scripts/verify_baseline_m0.sh`
3. If baseline intentionally changes, refresh references and move the baseline tag:
   - `baseline/m0-2026-02-20`
