# Tokamak Fusion Simulation — Milestone Checklist

## Project intent

This project is a kinetic plasma simulation prototype for a tokamak-like reactor geometry in C++. The current engine already supports charged-particle pushing with a Boris integrator, toroidal/poloidal magnetic field approximations, an electrostatic placeholder/closure path, spatial sorting, Monte Carlo-style fusion events, and telemetry.

This checklist defines milestones that move the simulator from a stable prototype into a numerically trustworthy and physically interpretable research tool.

## Current status snapshot

- Milestone 0: baseline run manifest support is available in CLI output.
- Milestone 1: implemented.
- Milestone 2: implemented in canonical runtime path (bounded profile model + CLI profile controls + runtime diagnostics + sidecar artifact export).
- Milestone 3: implemented in canonical runtime path (selectable electrostatic mode + residual-based SOR + CIC particle-mesh coupling + runtime diagnostics).
- Milestone 4+: pending.

---

## Milestone 0 — Baseline stability snapshot

### Goal
Freeze a known-good baseline before changing core physics.

### Tasks
- Tag the current version as a baseline checkpoint.
- Record a reference run configuration (seed, time step, total steps, scenario).
- Save baseline telemetry outputs for regression comparison.
- Add a run manifest to telemetry (time step, seed, particle cap, grid size, beam settings).

### Exit criteria
- You can re-run the same configuration and reproduce the same outputs (within stochastic variation if random seed changes).

### Implementation status
- Done: baseline capture workflow is defined in `scripts/capture_baseline_m0.sh`.
- Done: baseline verification workflow is defined in `scripts/verify_baseline_m0.sh`.
- Done: canonical deterministic baseline command is:
  - `./build/tokamakfusion --scenario cold --seed 20260220 --dt 1e-7 --steps 200 --telemetry-every 50 --artifacts-root ./baselines/m0/generated_runs`
- Done: baseline references are tracked under `baselines/m0/reference/`.

---

## Milestone 1 — Numerical correctness and safety

### Goal
Make the engine numerically trustworthy and remove silent failure modes.

### Tasks
- Randomize collision pairing within each occupied cell (do not pair in sorted-neighbor order only).
- Add deterministic random number generator (RNG) [random number generator] seeding mode for reproducible debug runs.
- Make `AddParticle(...)` return `bool` and handle insertion failures explicitly.
- Add telemetry counters for:
  - rejected injections
  - rejected fusion-product insertions
  - particle-cap hits
- Add debug assertions / guards for:
  - NaN [not a number] and Inf [infinity] values
  - out-of-bounds grid accesses
  - invalid particle state (negative mass, zero mass, invalid charge)
- Add energy and charge bookkeeping:
  - total kinetic energy
  - total macro-weighted charge
  - fusion energy deposited
  - beam energy injected

### Implementation status
- Done: randomized per-cell pairing via shuffle before pairing.
- Done: deterministic seed mode (`--seed`) with run-manifest seed output.
- Done: `AddParticle(...)` returns `bool`; multi-particle insertion pre-checks are explicit.
- Done: telemetry counters for rejected injection pairs, rejected fusion ash, and cap-hit events.
- Done: debug-time finite-state and bounds assertions are present.
- Done: energy/charge diagnostics are included in telemetry output.

### Exit criteria
- No silent particle drops.
- No NaN/Inf propagation in long runs.
- Same seed produces identical results.
- Reordering insertions does not materially change fusion trendlines.

---

## Milestone 2 — Magnetic field model cleanup

### Goal
Replace singular tokamak-current approximations with bounded profiles.

### Tasks
- Replace the singular poloidal field model with an enclosed-current profile:
  - Example: `I_enc(r) = I_p * (r/a)^2` for `r < a`
  - `I_enc(r) = I_p` for `r >= a`
- Implement an axis-safe evaluation for `B_pol` near `r = 0`.
- Add configurable plasma-current profiles (uniform, parabolic, custom function).
- Add field diagnostics:
  - radial magnetic field magnitude profile
  - maximum field magnitude per step
  - recommended time-step estimate from gyrofrequency

### Implementation status
- Done: singular poloidal behavior replaced by enclosed-current profile evaluation in `src/magnetic_field.cpp`.
- Done: axis-safe evaluation uses configurable epsilon near the magnetic axis.
- Done: profile selection is configurable from CLI/RunConfig (`uniform`, `parabolic`, `custom` with table input).
- Done: runtime telemetry includes per-step `Bmax_T` and gyro-based recommended time-step output.
- Done: radial magnetic-field diagnostics are exported to sidecar `magnetic_field_diagnostics_v2.csv`.
- Done: Milestone 7 artifact contract is versioned to `v2` summary/radial/manifest files.

### Exit criteria
- No artificial field singularity near the magnetic axis.
- Time-step constraints are no longer dominated by a non-physical `1/r` blow-up.

---

## Milestone 3 — Electrostatic solver and particle-mesh quality

### Goal
Upgrade the electric-field path from a placeholder closure to a controlled electrostatic particle-mesh loop.

### Tasks
- Keep the restoring-field placeholder as a selectable mode (`Placeholder`, `Electrostatic`).
- In electrostatic mode:
  - Use a residual-based stopping condition for SOR [successive over-relaxation] or the selected solver.
  - Log solver residual per step (or per solve call).
  - Make boundary conditions explicit and configurable.
- Upgrade charge deposition from NGP [nearest-grid-point] to CIC [cloud-in-cell].
- Upgrade field interpolation to matching CIC interpolation.
- Add grid-resolution and convergence tests for the electric field.

### Implementation status
- Done: CLI/RunConfig supports electric-field mode selection (`placeholder`, `electrostatic`), boundary condition selection, charge assignment selection, and SOR controls.
- Done: runtime electrostatic path uses particle-mesh coupling with selectable NGP/CIC deposition and matching interpolation.
- Done: SOR solve exposes per-step residual availability, convergence, iteration count, and tolerance through telemetry/artifacts.
- Done: sidecar artifact `electrostatic_diagnostics_v2.csv` is exported without changing existing v2 summary/radial contracts.
- Done: validation includes standalone electrostatic benchmark plus runtime manufactured-solution and refinement tests.

### Exit criteria
- Electric-field solver convergence is measured, not assumed.
- Field noise and numerical heating are reduced after CIC upgrade.
- Results are less sensitive to grid resolution.

---

## Milestone 4 — Collision and fusion realism v1

### Goal
Replace toy reaction logic with a physically interpretable reaction model.

### Tasks
- Replace hard fusion energy threshold with a continuous energy-dependent reactivity model:
  - `sigma(E)` [cross-section as a function of energy] or
  - simplified `⟨σv⟩(T)` [reactivity] table/fit
- Normalize reaction probability using:
  - relative speed
  - cell volume
  - species weights / number densities
  - time step
- Keep probability bounded with `P = 1 - exp(-lambda * dt)`.
- Add a weighted particle policy:
  - per-particle weight
  - partial depletion support
  - ash production with weight conservation
- Add diagnostics:
  - per-cell fusion event counts
  - fusion rate vs radius
  - fuel depletion vs ash accumulation

### Exit criteria
- Fusion rate changes smoothly with energy and density.
- Species and weight conservation are explicit and internally consistent.

---

## Milestone 5 — Boundary and wall interaction modes

### Goal
Move from debug-friendly reflections to interpretable tokamak wall behavior.

### Tasks
- Implement boundary modes:
  - `Reflect` (debug)
  - `Absorb` (particle loss)
  - `Recycle` (simplified plasma-wall recycling model)
- Add wall-hit diagnostics:
  - hit count by species
  - impact energy distribution
  - estimated wall power load
- Add optional particle loss accounting for confinement analysis.

### Exit criteria
- Boundary behavior is configurable by experiment.
- Wall losses are measurable and visible in telemetry.

---

## Milestone 6 — Validation suite (must-pass physics tests)

### Goal
Validate subsystems against analytic or controlled reference behaviors.

### Required tests
- **Uniform B-field gyro test**  
  Verify gyrofrequency and gyroradius against analytic values.
- **E×B [electric-cross-magnetic] drift test**  
  Verify drift velocity matches the expected `E × B / |B|^2`.
- **Wall reflection test**  
  Verify speed is conserved under reflective boundary mode.
- **Poisson / electrostatic solve test**  
  Compare against a known charge distribution case.
- **Collision reproducibility test**  
  Verify deterministic behavior under fixed seed.
- **Grid convergence test**  
  Check a chosen diagnostic against increasing grid resolution.

### Exit criteria
- All required tests pass within defined tolerances.
- Tolerances are documented.

### Implementation notes
- Poisson validation includes a standalone finite-difference numeric harness test.
- Engine-path electrostatic validation is covered by manufactured-solution and grid-refinement tests.

---

## Milestone 7 — Telemetry and analysis outputs

### Goal
Make the simulation outputs useful for scientific diagnosis and machine learning (ML) [machine learning] workflows.

### Tasks
- Add radial profile diagnostics:
  - density vs minor radius
  - average ion energy vs radius
  - fusion rate vs radius
- Add distribution diagnostics:
  - speed histogram
  - pitch-angle histogram
- Export structured outputs:
  - CSV [comma-separated values] summaries
  - periodic particle snapshots
  - solver residual logs
- Save a machine-readable run configuration file (JSON [JavaScript Object Notation] or YAML [YAML Ain't Markup Language]).

### Exit criteria
- A single run produces enough artifacts to diagnose confinement, heating, and losses without code inspection.

---

## Milestone 8 — Performance and scaling

### Goal
Scale particle count and run length after correctness/fidelity is established.

### Tasks
- Remove per-step allocations in hot paths (reuse buffers).
- Reserve all particle arrays to `MAX_PARTICLES`.
- Benchmark `Vec3` layout vs split scalar arrays (`x[]`, `y[]`, `z[]`).
- Parallelize:
  - particle push
  - grid histogram and scatter
  - collision processing by cell (with deferred mutation queues)
- Add a benchmark harness:
  - particles/sec
  - steps/sec
  - per-stage timing (push, sort, field solve, collisions)

### Exit criteria
- Performance regressions are detectable.
- Scaling behavior is quantified, not guessed.

---

## Milestone 9 — Codex-agent-ready development workflow

### Goal
Make the repository easy for a coding agent to modify safely.

### Tasks
- Add `AGENT.md` with project architecture, assumptions, coding rules, and task workflow.
- Add `PHYSICS_MANIFEST.md` with component-by-component modeling assumptions.
- Add `tests/` and `benchmarks/` folders with documented entry points.
- Add standard run targets (debug, validation, benchmark).
- Add a `CONTRIBUTING.md` with style and validation requirements.

### Exit criteria
- A coding agent can implement a component change and run the correct validations without human guesswork.

---

## Suggested execution order

1. Milestone 1 (correctness and safety)  
2. Milestone 2 (bounded magnetic field model)  
3. Milestone 3 (electrostatic + CIC)  
4. Milestone 4 (fusion realism and weights)  
5. Milestone 6 (validation suite)  
6. Milestone 7 (telemetry outputs)  
7. Milestone 5 (wall modes)  
8. Milestone 8 (performance scaling)  
9. Milestone 9 (workflow hardening)

This order prioritizes scientific trust before optimization.

---

## Acronyms used

RNG [random number generator], NaN [not a number], Inf [infinity], SOR [successive over-relaxation], NGP [nearest-grid-point], CIC [cloud-in-cell], ML [machine learning], CSV [comma-separated values], JSON [JavaScript Object Notation], YAML [YAML Ain't Markup Language].
