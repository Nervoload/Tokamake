# ⚛ Tokamake: Tokamak Fusion Physics Simulator

**To learn about the physics underlying fusion reactions, and to play with synthetic data for ML models, I made Tokamake--a physics engine for simulating fusion dynamics 😎**

### Currently supporting 3D OpenGL visualization, stepwise telemetry, and producing approximate particle statistics.

## Overview

Tokamake at it's core is a kinetic plasma simulation that models the behavior of charged particles in a tokamak-like magnetic geometry. The engine supports particle pushing, magnetic field approximations, electrostatics, and fusion events.
## Core Features

- **Particle Dynamics**: Charged-particle pushing with a Boris integrator.
- **Field Models**: Toroidal and poloidal magnetic field approximations with bounded, axis-safe profiles.
- **Electrostatics**: A selectable electrostatic solver (Successive Over-Relaxation) with Cloud-in-Cell (CIC) particle-mesh coupling.
- **Fusion Events**: Monte Carlo-style fusion event modeling.
- **Diagnostics**: Detailed telemetry and diagnostics with structured sidecar artifact exports in CSV and JSON formats.
- **Visualization**: A replay viewer (`tokamak_viewer`) for visualizing simulation artifacts from run manifests.

## Project Goals & Roadmap

The development is organized into milestones to incrementally build a trustworthy research tool.

- **Milestone 0: Baseline Stability**
  - **Goal**: Freeze a known-good baseline for regression testing.
  - **Status**: ✅ Implemented

- **Milestone 1: Numerical Correctness & Safety**
  - **Goal**: Ensure the engine is numerically trustworthy and free of silent failure modes.
  - **Status**: ✅ Implemented

- **Milestone 2: Magnetic Field Model Cleanup**
  - **Goal**: Replace singular field approximations with bounded, configurable plasma current profiles.
  - **Status**: ✅ Implemented

- **Milestone 3: Electrostatic Solver & Particle-Mesh Quality**
  - **Goal**: Upgrade the electric field path to a controlled particle-mesh loop with CIC coupling.
  - **Status**: ✅ Implemented

- **Milestone 4: Collision and Fusion Realism**
  - **Goal**: Implement a continuous, energy-dependent fusion reactivity model and weighted particle policies.
  - **Status**: ⏳ Pending

- **Milestone 5: Boundary and Wall Interactions**
  - **Goal**: Add configurable wall behaviors like absorption and simplified plasma-wall recycling.
  - **Status**: ⏳ Pending

- **Milestone 6: Validation Suite**
  - **Goal**: Validate subsystems against analytic solutions (e.g., gyro-motion, E×B drift).
  - **Status**: ⏳ Pending

- **Milestone 7: Telemetry and Analysis Outputs**
  - **Goal**: Enhance outputs for scientific diagnosis and ML workflows.
  - **Status**: ✅ Partially Implemented (Schema v2)

- **Milestone 8: Performance and Scaling**
  - **Goal**: Optimize hot paths and parallelize core components to scale particle count and run length.
  - **Status**: ⏳ Pending

- **Milestone 9: Agent-Ready Workflow**
  - **Goal**: Improve documentation and repository structure for safe modification by autonomous coding agents.
  - **Status**: ⏳ Pending

# Getting Started

### Dependencies

- A C++ compiler with C++17 support.
- CMake (for building the project).
- **GLFW** (for the optional viewer). On macOS, you can install it via Homebrew: `brew install glfw`.

### Building the Simulation

1.  **Configure the project** (from the root directory):
    ```bash
    cmake -S . -B build
    ```

2.  **Build the simulation engine**:
    ```bash
    cmake --build build --target tokamakfusion
    ```

### Running a Simulation

You can run a simulation from the command line. Artifacts will be saved to the specified directory.

```bash
# Example deterministic run from Milestone 0
./build/tokamakfusion --scenario cold --seed 20260220 --dt 1e-7 --steps 200 --telemetry-every 50 --artifacts-root ./output/runs

# Example run to generate replay artifacts for the viewer
./build/tokamakfusion \
  --scenario <scenario_name> \
  --steps 400 \
  --telemetry-every 50 \
  --artifact-every 50 \
  --particle-snapshot-every 50 \
  --artifacts-root ./output/runs
```

## Visualization (`tokamak_viewer`)

The project includes a replay viewer to visualize simulation artifacts.

### Building the Viewer

1.  **Vendor Dependencies**: The viewer requires ImGui and GLAD. Please follow the instructions in `docs/VIEWER_SETUP.md` to place the required files in the `third_party/` directory.

2.  **Configure with viewer enabled**:
    ```bash
    cmake -S . -B build -DBUILD_VIEWER=ON
    ```

3.  **Build the viewer target**:
    ```bash
    cmake --build build --target tokamak_viewer
    ```

### Launching the Viewer

You can launch the viewer by pointing it to a run's manifest file or its directory.

```bash
# By manifest path
./build/tokamak_viewer --manifest ./output/runs/<run_id>/manifest_v2.json

# By run directory
./build/tokamak_viewer --run-dir ./output/runs/<run_id>
```

Thank you for checking out this project!
