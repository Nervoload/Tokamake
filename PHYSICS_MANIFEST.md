# PHYSICS_MANIFEST.md — Component Assumptions and Validation Map

## Purpose

This file documents, for each major component of the tokamak simulation, what is physically modeled, what is currently a placeholder, what is intentionally non-physical for stability or debugging, and which diagnostics/tests should validate the component.

This is the scientific honesty layer for the project.

---

## 1) Magnetic field model

### Physically modeled
- A tokamak-like magnetic geometry with:
  - toroidal field component from coil current (Ampère-based approximation)
  - poloidal field component from plasma current (approximate enclosed-current formulation)
- Position-dependent magnetic field used by the particle pusher.
- Configurable enclosed-current profiles (`uniform`, `parabolic`, `custom table`) through run configuration / CLI.

### Placeholder / approximation
- The field is not derived from a full axisymmetric equilibrium solve (no Grad-Shafranov [axisymmetric magnetostatic equilibrium equation]).
- Plasma shaping (elongation, triangularity, divertor geometry) is not modeled unless explicitly added.

### Intentionally non-physical for stability
- Any smoothing, clamping, or bounded-axis regularization used to avoid singular fields near `r = 0`.
- Any field-strength caps used for time-step safety during debug runs.

### Diagnostics / validation
- Radial magnetic field magnitude profile (`|B|` vs minor radius).
- Max/mean `|B|` telemetry per step.
- Gyro-based recommended timestep diagnostic from per-step max `|B|`.
- Sidecar artifact export: `magnetic_field_diagnostics_v2.csv` (included in Milestone 7 v2 manifest contract).
- Uniform-field gyro test (compare gyrofrequency and gyroradius to analytic values).
- Time-step recommendation diagnostic from max gyrofrequency.

---

## 2) Electrostatic field / charge solve

### Physically modeled
- A charge-density-driven electrostatic potential solve (Poisson-like or screened Poisson form).
- Electric field sampled by ions through particle-mesh coupling.
- Selectable electric-field runtime mode (`placeholder` vs `electrostatic`).
- Particle-mesh coupling supports NGP [nearest-grid-point] and CIC [cloud-in-cell] assignment/interpolation.
- Residual-based SOR [successive over-relaxation] stopping in electrostatic mode.

### Placeholder / approximation
- Simplified geometry terms (e.g., omission of full cylindrical Laplacian terms) if not yet implemented.
- Simplified electron closure or neutralizing background profile.
- Fixed Debye-like screening or heuristic screening length.
- Uniform Cartesian mesh is used as first integration target before geometry-aware electrostatic operators.

### Intentionally non-physical for stability
- Placeholder restoring electric field mode that acts like a confining spring.
- Artificial neutralization to suppress runaway charge separation in early prototypes.
- Solver damping or under-relaxation added purely for convergence stability.

### Diagnostics / validation
- Solver residual norm vs iteration (for SOR [successive over-relaxation] or equivalent solver).
- Electric field magnitude histograms.
- Potential and field snapshots on the grid.
- Sidecar artifact export: `electrostatic_diagnostics_v2.csv` (Milestone 7 v2 artifact contract).
- Known-charge-distribution standalone finite-difference harness test case.
- Engine-path runtime electrostatic solver validation test with manufactured 3D Dirichlet solution.
- Grid refinement test for field convergence.

---

## 3) Particle push (Boris integrator)

### Physically modeled
- Lorentz-force motion of charged ions under electric and magnetic fields.
- Time integration with Boris method (good geometric behavior for gyro motion).

### Placeholder / approximation
- Non-relativistic dynamics unless explicitly upgraded.
- Ion-only kinetics if electrons are not explicit.

### Intentionally non-physical for stability
- Time-step restrictions selected conservatively for numerical stability.
- Optional debug clamps on extreme velocity or field values (if enabled).

### Diagnostics / validation
- Uniform B-field gyro orbit test.
- E×B [electric-cross-magnetic] drift test.
- Speed conservation in pure magnetic field (no electric field, no collisions).
- NaN [not a number] / Inf [infinity] guards on positions and velocities.

---

## 4) Particle system and weights

### Physically modeled
- Macro-particle representation of ion populations.
- Species-specific mass and charge.
- Optional per-particle weight for population scaling.

### Placeholder / approximation
- Global uniform macro-weight if per-particle weights are not yet implemented.
- Simplified ash insertion and depletion bookkeeping in early fusion stages.

### Intentionally non-physical for stability
- Hard particle cap (`MAX_PARTICLES`) to prevent memory blow-up.
- Temporary compaction strategies chosen for safety before performance optimization.

### Diagnostics / validation
- Total macro-weighted particle count by species.
- Total macro-weighted charge.
- Particle insertion failure counters.
- Compaction consistency checks (array lengths remain synchronized).

---

## 5) Spatial grid and particle sorting

### Physically modeled
- Spatial locality for collisions and particle-mesh coupling via a 3D Cartesian [three-dimensional rectangular] grid.
- Histogram + prefix sum + scatter organization for efficient cell access.

### Placeholder / approximation
- Cartesian grid in a toroidal domain (many cells are geometrically unused).
- Uniform cell size even where plasma occupancy is non-uniform.

### Intentionally non-physical for stability
- Cell-index clamping at domain boundaries to avoid out-of-bounds access.
- Conservative grid resolution choices to protect memory.

### Diagnostics / validation
- Per-cell occupancy statistics (mean, max, occupancy histogram).
- Sort consistency test: every particle appears exactly once after scatter.
- Bounds checks for cell indices.
- Per-stage timing (histogram, prefix sum, scatter).

---

## 6) Collisions and fusion reactions

### Physically modeled
- Cell-local stochastic pairing and reaction handling.
- Per-cell in-place shuffle before pairing to reduce sorted-order bias.
- Two-phase processing: event selection first, deferred mutation second.
- D-T [deuterium-tritium] fusion event pathway with alpha [helium nucleus] ash creation.
- Bounded event probability (`1 - exp(-lambda * dt)` style).

### Placeholder / approximation
- Simplified reactivity model (`lambda` or toy cross-section) before a validated `sigma(E)` [cross-section as a function of energy] model is added.
- Pairwise local collisions instead of a full Coulomb collision operator.
- Simplified neutron handling (usually omitted from ion state).

### Intentionally non-physical for stability
- Collision pairing simplifications (e.g., one-pass pairing, capped per-cell events) used to avoid index hazards.
- Deferred mutation and compaction after collision passes.
- Hard threshold fusion trigger (`E_k > 15 keV`) remains a placeholder in this phase.

### Diagnostics / validation
- Fusion events per step and per radius bin.
- Per-cell max reaction counter.
- Fuel depletion and ash accumulation curves.
- Fixed-seed reproducibility checks.
- Sensitivity checks against pair-order randomization.
- Energy accounting for alpha energy deposition.

---

## 7) Neutral beam injection (NBI)

### Physically modeled
- Directed injection of energetic ions into the tokamak domain.
- Timestep-invariant source rate using an accumulator.

### Placeholder / approximation
- Beam ionization, neutral penetration, and deposition profiles are simplified.
- Injection geometry spread is heuristic unless calibrated.

### Intentionally non-physical for stability
- Injection clamps at capacity (`MAX_PARTICLES`) to protect memory.
- Simplified paired D/T injection for balanced fueling in prototype runs.

### Safety behavior implemented
- NBI insertion is atomic at pair granularity (`CanInsert(2)` pre-check).
- Rejected pair insertions are surfaced in telemetry counters (no silent drops).

### Diagnostics / validation
- Injected particles per step and cumulative total.
- Rejected injection counter.
- Beam power estimate (from injected kinetic energy).
- Sensitivity test to time-step changes (source rate should remain invariant).

---

## 8) Boundary and wall interactions

### Physically modeled
- Tokamak-like toroidal wall geometry.
- Particle-wall intersection detection.

### Placeholder / approximation
- Perfectly elastic reflection mode for debugging.
- Simplified absorb/recycle modes if implemented later.

### Intentionally non-physical for stability
- Reflective wall mode is intentionally non-physical and used to preserve particle count while debugging confinement behavior.

### Diagnostics / validation
- Wall-hit count by species.
- Impact energy distribution.
- Speed conservation check in reflective mode.
- Wall loss power (for absorb mode).

---

## 9) Telemetry and analysis outputs

### Physically modeled
- Observable summaries of simulated system state (counts, energies, fusion events, profiles).

### Placeholder / approximation
- “Average ion energy” used as a proxy while true temperature decomposition is not implemented.
- Limited profile set in early versions.

### Intentionally non-physical for stability
- Reduced telemetry cadence for performance (sampling instead of every-step output).

### Diagnostics / validation
- Consistency checks across console output and exported files.
- Schema checks for CSV [comma-separated values] / JSON [JavaScript Object Notation] output.
- Regression comparisons against baseline runs.

---

## 10) Validation and benchmarking framework

### Physically modeled
- Not a physics component directly; this is the trust framework used to verify physics and numerics.

### Placeholder / approximation
- Tolerances may start broad and be tightened over time.
- Benchmark scenarios may be simplified before full-physics tests exist.

### Intentionally non-physical for stability
- Controlled synthetic tests (uniform fields, fixed distributions) are intentionally simplified to isolate subsystem behavior.

### Diagnostics / validation
- Required test suite:
  - uniform B-field gyro test
  - E×B drift test
  - electrostatic solver test
  - collision reproducibility test
  - grid consistency test
- Benchmark suite:
  - particles per second
  - stage timing
  - memory footprint

---

## How to use this file

When changing any component:
1. Update the relevant section if assumptions or approximations changed.
2. Add or revise diagnostics if the component behavior changes.
3. Keep placeholder and intentionally non-physical behaviors explicitly documented.

This prevents accidental “physics inflation” in code comments or reports.

---

## Acronyms used

SOR [successive over-relaxation], E×B [electric-cross-magnetic], NaN [not a number], Inf [infinity], D-T [deuterium-tritium], NBI [neutral beam injection], CSV [comma-separated values], JSON [JavaScript Object Notation].
