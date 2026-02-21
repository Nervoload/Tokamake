# Viewer Setup (`tokamak_viewer`)

## Overview
`tokamak_viewer` replays artifact schema v2 runs from `manifest_v2.json` + snapshot CSV files.

## Build Flags
- `-DBUILD_VIEWER=ON`: enables viewer target.
- `-DVIEWER_STRICT_VENDOR_CHECK=ON`: fail configure if required ImGui/GLAD vendor files are missing.

Default CI path keeps viewer off (`BUILD_VIEWER=OFF`).

## How to Vendor Dependencies

The viewer requires several third-party libraries (ImGui, GLAD) that must be "vendored" — manually placed into the `third_party/` directory before building. Here are the instructions to get the required files.

### 1. Vendor GLAD (OpenGL Loader)

1.  Go to the GLAD web generator.
2.  Set the options as follows:
    -   **Language**: `C/C++`
    -   **Specification**: `OpenGL`
    -   **API gl**: Version `3.3`
    -   **Profile**: `Core`
    -   **Options**: Ensure `Generate a loader` is checked.
3.  Click **Generate**.
4.  On the next page, download the generated `glad.zip` file.
5.  Extract the archive. You will find `include` and `src` folders.
6.  Copy the contents of the extracted `include` folder (which are the `glad` and `KHR` directories) into this project's `third_party/glad/include/`.
7.  Create a `src` directory inside `third_party/glad/` and copy the `src/glad.c` file into it.

### 2. Vendor ImGui (GUI Library)

1.  Go to the Dear ImGui GitHub repository and download the source code, usually via the "Code" -> "Download ZIP" button.
2.  Extract the downloaded archive.
3.  From the extracted folder, copy the required files into the `third_party/imgui/` directory of this project. Create subdirectories as needed.

    **Core files (copy to `third_party/imgui/`):**
    -   `imgui.cpp`
    -   `imgui_draw.cpp`
    -   `imgui_tables.cpp`
    -   `imgui_widgets.cpp`
    -   `imgui.h`

    **Backend files (copy to `third_party/imgui/backends/`):**
    -   `backends/imgui_impl_glfw.cpp`
    -   `backends/imgui_impl_glfw.h`
    -   `backends/imgui_impl_opengl3.cpp`
    -   `backends/imgui_impl_opengl3.h`

After copying, your `third_party/imgui` directory should contain the core files, and a `backends` subdirectory with the platform-specific files.

---

## Required Dependencies (macOS first)

This section lists the dependencies you need to have installed or vendored.

1. GLFW (system/Homebrew):
- `brew install glfw`

2. Vendored GLAD files under `third_party/glad`:
- `third_party/glad/src/glad.c`
- `third_party/glad/include/glad/glad.h`
- `third_party/glad/include/KHR/khrplatform.h`

3. Vendored ImGui files under `third_party/imgui`:
- `third_party/imgui/imgui.cpp`
- `third_party/imgui/imgui_draw.cpp`
- `third_party/imgui/imgui_tables.cpp`
- `third_party/imgui/imgui_widgets.cpp`
- `third_party/imgui/imgui.h`
- `third_party/imgui/backends/imgui_impl_glfw.cpp`
- `third_party/imgui/backends/imgui_impl_glfw.h`
- `third_party/imgui/backends/imgui_impl_opengl3.cpp`
- `third_party/imgui/backends/imgui_impl_opengl3.h`

## Build
```bash
cmake -S /Users/johnsurette/Documents/Codespace🪐/Tokamak -B /Users/johnsurette/Documents/Codespace🪐/Tokamak/build -DBUILD_VIEWER=ON
cmake --build /Users/johnsurette/Documents/Codespace🪐/Tokamak/build --target tokamak_viewer
```

## Generate A Replay Run
```bash
/Users/johnsurette/Documents/Codespace🪐/Tokamak/build/tokamakfusion \
  --scenario <scenario_name> \
  --seed 20260220 \
  --steps 400 \
  --telemetry-every 50 \
  --artifact-every 50 \
  --particle-snapshot-every 50 \
  --artifacts-root /Users/johnsurette/Documents/Codespace🪐/Tokamak/output/runs
```

Find latest run:
```bash
ls -1dt /Users/johnsurette/Documents/Codespace🪐/Tokamak/output/runs/run_* | head -n1
```

## Launch Viewer
By manifest path:
```bash
/Users/johnsurette/Documents/Codespace🪐/Tokamak/build/tokamak_viewer \
  --manifest /Users/johnsurette/Documents/Codespace🪐/Tokamak/output/runs/<run_id>/manifest_v2.json
```

By run directory:
```bash
/Users/johnsurette/Documents/Codespace🪐/Tokamak/build/tokamak_viewer \
  --run-dir /Users/johnsurette/Documents/Codespace🪐/Tokamak/output/runs/<run_id>
```

Optional flags:
- `--point-size 3.0`
- `--playback-rate 1.0`
- `--start-step 0`

## Troubleshooting
- Configure fails with missing vendor files: add the required files listed above.
- Viewer opens but no frames: verify `manifest_v2.json` contains non-empty `files.particle_snapshot_csv_files` and referenced files exist.
- Viewer target not found: confirm `-DBUILD_VIEWER=ON` during configure.
