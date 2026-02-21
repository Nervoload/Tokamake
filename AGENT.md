# AGENT.md — Codex Agent Context for Tokamak Fusion Simulation

## Purpose of this repository

This repository contains a C++ kinetic plasma simulation prototype for a tokamak-like fusion reactor geometry. The code is not yet a full self-consistent tokamak solver. It is a staged engineering project intended to evolve from a numerically stable charged-particle engine into a more physically credible plasma simulation.

The primary goals are:

1. Build a trustworthy charged-particle simulation core (Boris integrator, particle-grid sorting, bounded fields, safe memory behavior).
2. Incrementally improve physics fidelity (electrostatics, collisions, fusion reactivity, wall interactions).
3. Maintain strong telemetry and diagnostics for scientific debugging and future machine learning (ML) [machine learning] analysis.
4. Keep the implementation modular and testable so an autonomous coding agent can make safe changes.

---

## Milestone implementation status (current repo)

- Milestone 1 (numerical correctness and safety) has been implemented.
- Collision pairing is randomized per occupied cell with deterministic fixed-seed support.
- Particle insertion is explicit and failable (`AddParticle(...) -> bool`) with rejection counters.
- Collision selection and mutation are split into separate phases.
- Grid sorting reuses persistent write-head buffers (no per-step temporary cursor allocation).
- Runtime counters include out-of-domain grid clamp tracking.
- Artifact contract is versioned as Milestone 7 schema `v2` (`manifest_v2.json`, `summary_v2.csv`, sidecar CSVs).

---

## Repository layout (current)

- `CMakeLists.txt`
- `include/tokamak/types.hpp`
- `include/tokamak/config.hpp`
- `include/tokamak/particle_system.hpp`
- `include/tokamak/spatial_grid.hpp`
- `include/tokamak/collision.hpp`
- `include/tokamak/diagnostics.hpp`
- `include/tokamak/engine.hpp`
- `include/tokamak/viewer/replay_manifest.hpp`
- `include/tokamak/viewer/replay_snapshot.hpp`
- `include/tokamak/viewer/replay_loader.hpp`
- `include/tokamak/viewer/viewer_app.hpp`
- `src/main.cpp`
- `src/viewer/main.cpp`
- `src/viewer/replay_loader.cpp`
- `src/viewer/replay_manifest.cpp`
- `src/viewer/replay_snapshot.cpp`
- `src/viewer/viewer_app.cpp`
- `src/viewer/gl_renderer.cpp`
- `src/viewer/camera.cpp`
- `src/engine.cpp`
- `src/particle_system.cpp`
- `src/spatial_grid.cpp`
- `src/collision.cpp`
- `src/telemetry.cpp`
- `tests/particle_system_tests.cpp`
- `tests/engine_tests.cpp`
- `tests/integration_tests.cpp`
- `third_party/googletest/` (local offline test harness compatible with GoogleTest-style APIs)

Legacy monolith and viewer/snapshot files are archived under `legacy/` and are not part of the active build graph.
Canonical visualization is `tokamak_viewer` (artifact replay, schema v2).

---

## Current simulation scope

### What exists (high level)
- Toroidal and poloidal magnetic field approximation in tokamak geometry.
- Particle push via Boris integrator.
- Particle storage in a data-oriented style (struct-of-arrays-inspired layout).
- Spatial binning / counting-sort style grid for cell-local operations.
- Monte Carlo-style fusion event scaffold.
- Electrostatic / restoring-field placeholder and evolving electrostatic solve path.
- Telemetry for ion counts, energy, and fusion events.

### What does **not** exist yet (important)
- Full Particle-in-Cell (PIC) [particle-in-cell] self-consistent electromagnetic field solve.
- Explicit electron kinetics.
- Full Coulomb collision operator (Landau / Fokker-Planck [kinetic diffusion-drift equation] treatment).
- Real tokamak equilibrium solver (e.g., Grad-Shafranov [axisymmetric magnetostatic equilibrium equation]).
- Edge / divertor / neutral recycling / radiation models.

A coding agent should preserve this distinction and avoid implying full tokamak fidelity where placeholders are still in use.

---

## Guiding engineering principles

### 1) Correctness before speed
Do not optimize a component until:
- it has diagnostics,
- it has a basic validation test,
- and its assumptions are documented.

### 2) Scientific honesty
When a component is a placeholder, label it clearly in code comments and telemetry names. Avoid names that imply a physically complete model if the implementation is a surrogate.

### 3) Determinism for debugging
All stochastic systems (collision pairing, Monte Carlo fusion, injections with randomness) should support a fixed-seed mode. Reproducible debugging is mandatory.

### 4) No silent failures
Do not drop particles, skip products, or clamp values silently without telemetry counters or explicit logs. Silent behavior is difficult to debug and corrupts scientific interpretation.

### 5) Deferred mutation in cell loops
When processing cell-local collisions or reactions, avoid mutating global particle storage in ways that invalidate indices during iteration. Prefer:
- marking for deletion,
- deferred insertion queues,
- compaction after the loop.

---

## Repository-level implementation expectations for a coding agent

## Coding style and safety
- Prefer `double` for physics calculations in strict SI [International System of Units] units unless there is a documented nondimensionalization.
- Prefer explicit units in variable names where practical (e.g., `_m`, `_s`, `_J`, `_keV`).
- Prefer `std::unique_ptr` or direct member ownership over raw `new` / `delete`.
- Return status values (`bool`, error codes, or structured results) for operations that can fail (e.g., particle insertion at capacity).
- Add comments for any non-obvious numerical or physical approximation.

## Performance style
- Avoid per-step heap allocations in hot loops.
- Reuse buffers for grid sorting and temporary arrays.
- Keep data layout decisions benchmarkable (do not over-assume auto-vectorization).
- Add timing instrumentation around major phases:
  - field solve
  - particle push
  - grid sort
  - collisions

## Testing style
Every nontrivial component change should include at least one of:
- a validation test,
- a regression test,
- or a benchmark delta with explanation.

If a change affects physics behavior, update the relevant diagnostics or validation target.

---

## Component map (for agent navigation)

The project is organized conceptually into these components. File names may evolve, but the logical split should remain.

1. **Math and constants**
   - Physical constants and vector math utilities.
   - Common source of unit consistency issues.

2. **Configuration and scenarios**
   - Tokamak geometry, magnetic field settings, injection rates, run mode toggles.
   - Should remain serializable for reproducibility.

3. **Particle system**
   - Particle arrays, species, weights, insertion/deletion/compaction.
   - Critical for memory safety and conservation bookkeeping.

4. **Spatial grid / particle sorting**
   - Histogram + prefix sum + scatter (counting sort).
   - Foundation for collisions and particle-mesh coupling.

5. **Field models**
   - Magnetic field model (toroidal + poloidal).
   - Electrostatic closure / Poisson path.
   - Placeholder terms must be clearly separated from validated physics.

6. **Particle push**
   - Boris integrator and boundary handling.
   - Performance-critical and numerically sensitive.

7. **Collisions and fusion**
   - Pairing, reaction probabilities, species conversion, weighted depletion.
   - Must support deterministic debug mode.

8. **Telemetry and outputs**
   - Console diagnostics, exported profiles, snapshots, validation summaries.
   - Includes out-of-domain clamp counters and artifact schema `v2` manifests.
   - Must be informative enough to diagnose behavior without reading code.

9. **Viewer replay tooling**
   - `tokamak_viewer` replays `manifest_v2.json` + snapshot CSV artifacts.
   - Legacy `.tksnap` viewer path is archived only.

10. **Validation and benchmarks**
   - Analytic tests and runtime benchmarks.
   - Required for safe iteration.

---

## Development workflow for a coding agent

When implementing a change, follow this workflow:

Canonical command entry points for this workflow are documented in:
- `CONTRIBUTING.md`
- `tests/README.md`
- `benchmarks/README.md`

### A) Read context
- Read `PHYSICS_MANIFEST.md` first to understand what is modeled vs placeholder.
- Read `MILESTONE_CHECKLIST.md` to identify priority and acceptance criteria.
- Read the relevant component code and existing diagnostics.

### B) Make a minimal scoped change
- Change one subsystem at a time (e.g., collision pairing randomization, CIC deposition).
- Avoid mixing physics changes and performance refactors in the same patch unless necessary.

### C) Add or update diagnostics
- If the change affects physics, add telemetry or validation outputs that prove the effect.
- If the change affects performance, add timing instrumentation or benchmark comparison.

### D) Validate
- Run the smallest relevant test scenario.
- Use `./scripts/run_debug.sh` for a deterministic debug smoke run.
- Use `./tests/run_validation.sh` for the validation suite.
- Use `./benchmarks/run_benchmark.sh` when performance-sensitive paths are changed.
- If stochastic, run fixed-seed and compare against baseline behavior.
- Note expected changes vs suspicious regressions.

### E) Document assumptions
- Update `PHYSICS_MANIFEST.md` if modeling assumptions changed.
- Update milestone progress if acceptance criteria are met.

---

## Common pitfalls (do not repeat)

1. **Index invalidation in collision loops**
   Do not remove particles directly while iterating over sorted cell particle IDs.

2. **Timestep-dependent source terms**
   Injection and source terms must scale with `dt` through rates and accumulators.

3. **Misleading labels**
   Do not label average kinetic energy as temperature unless thermalization and bulk-flow subtraction are handled correctly.

4. **Singular field models near axes**
   Avoid `1/r` singular behavior when a bounded enclosed-current profile is physically appropriate.

5. **Silent capacity clamps**
   Do not silently drop inserted particles at `MAX_PARTICLES`; expose counters and explicit policy.

6. **Unverified vectorization claims**
   Do not assume compiler SIMD [single instruction, multiple data] behavior without profiling evidence.

---

## Priority queue for upcoming improvements

The recommended next tasks for a coding agent are:

1. Randomize collision pairing within cells.
2. Make all particle insertions explicit/failable with telemetry.
3. Replace singular poloidal field with enclosed-current profile.
4. Improve electrostatic solver convergence checks.
5. Upgrade deposition/interpolation to CIC [cloud-in-cell].
6. Add validation tests for gyro motion and E×B [electric-cross-magnetic] drift.
7. Add benchmark harness and hot-path timing.

These are the highest-leverage changes for scientific trust and development speed.

---

## Output expectations for agent-generated patches

A good patch should include:
- Code changes
- A short explanation of what changed and why
- Validation notes (what was run, what passed)
- Any known limitations left in place
- Updates to docs if assumptions changed

The goal is to make every patch auditable by a human scientist/engineer with minimal guesswork.

---

## Acronyms used

ML [machine learning], PIC [particle-in-cell], Fokker-Planck [kinetic diffusion-drift equation], SI [International System of Units], SIMD [single instruction, multiple data], CIC [cloud-in-cell], E×B [electric-cross-magnetic].
