# Milestone 7 Output Schema (v2)

Milestone 7 artifacts are emitted per run under:

`output/runs/run_<UTC_TIMESTAMP>_<SCENARIO>_seed<ACTIVE_SEED>_ns<UNIX_NS>/`

All artifact files include schema version `2` either in filename (`*_v2.*`) or first CSV column (`schema_version`).

## File naming conventions

- `manifest_v2.json`: artifact index for the run
- `run_config_v2.json`: machine-readable run configuration used for the run
- `summary_v2.csv`: run-level telemetry time series
- `radial_profiles_v2.csv`: radial diagnostics (density, average ion energy, fusion rate)
- `magnetic_field_diagnostics_v2.csv`: radial magnetic diagnostics per sampled step
- `electrostatic_diagnostics_v2.csv`: electric-field diagnostics and solver status per sampled step
- `speed_histogram_v2.csv`: speed distribution histogram per sampled step
- `pitch_angle_histogram_v2.csv`: pitch-angle distribution histogram per sampled step
- `solver_residuals_v2.csv`: solver residual log contract rows (placeholder-compatible)
- `snapshots/particles_step_<8-digit-step>.csv`: periodic particle snapshots

## Field semantics

### `summary_v2.csv`
- `schema_version`: artifact schema version
- `step`, `time_s`, `active_seed`, `scenario`
- Species/energy/fusion counters: `total_ions`, `deuterium`, `tritium`, `helium`, `avg_energy_kev`, `fusion_events_total`
- Safety/runtime counters: `particle_cap_hit_events`, `rejected_injection_pairs`, `rejected_fusion_ash`, `out_of_domain_cell_clamp_events`, `fusion_attempts`, `fusion_accepted`, `max_reactions_in_cell`
- Energy/charge budget: `kinetic_j`, `beam_injected_j`, `fusion_alpha_injected_j`, `total_charge_c`

### `radial_profiles_v2.csv`
- Bin geometry: `bin_index`, `r_inner_m`, `r_outer_m`, `r_center_m`
- Plasma profile values: `ion_count`, `macro_weight`, `shell_volume_m3`, `density_m3`, `avg_ion_energy_kev`
- Fusion profile values: `fusion_events_cumulative`, `fusion_rate_m3_s`, `fusion_rate_placeholder`
- If `fusion_rate_placeholder=1`, `fusion_rate_m3_s` is placeholder-compatible and not physically populated.

### `magnetic_field_diagnostics_v2.csv`
- Bin geometry: `bin_index`, `r_inner_m`, `r_outer_m`, `r_center_m`
- Field diagnostics: `mean_b_t`, `sample_count`, `step_max_b_t`, `recommended_dt_s`, `profile_kind`

### `electrostatic_diagnostics_v2.csv`
- Per-step solve diagnostics: `electric_field_mode`, `boundary_condition`, `charge_assignment`
- Field magnitude diagnostics: `max_electric_field_v_per_m`, `mean_electric_field_v_per_m`
- Solver diagnostics: `solver_iterations`, `solver_converged`, `residual_l2`

### `speed_histogram_v2.csv`
- Histogram geometry: `bin_index`, `speed_min_m_per_s`, `speed_max_m_per_s`
- Counts: `count`, `total_samples`

### `pitch_angle_histogram_v2.csv`
- Histogram geometry: `bin_index`, `pitch_min_deg`, `pitch_max_deg`
- Counts/quality: `count`, `total_samples`, `invalid_samples`
- `invalid_samples` are particles where speed or local B-field magnitude is near zero.

### `solver_residuals_v2.csv`
- Contract fields:
  - `residual_available`, `residual_l2`, `solver_name`, `status`
  - `iterations`, `converged`, `tolerance`, `note`
- In the current placeholder electric field path, rows are written with:
  - `residual_available=false`
  - `residual_l2=nan`
  - `solver_name=none`
  - `status=placeholder`
  - `iterations=0`
  - `converged=false`
  - `tolerance=nan`
- This keeps downstream workflows stable until a residual-producing electrostatic solver is integrated.

### `run_config_v2.json`
Contains complete run and model configuration used for reproduction:
- run id/time metadata
- scenario/seed/time-step/steps/cadence/particle cap
- tokamak and NBI parameters
- artifact exporter settings

### `manifest_v2.json`
Top-level index for all files created in a run directory. Downstream ingestion should use this as the canonical entrypoint.
The `files` object includes `magnetic_field_diagnostics_csv` in v2.
The `files` object includes `electrostatic_diagnostics_csv` in v2.
