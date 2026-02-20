# Legacy Archive

This directory stores deprecated simulation/viewer code preserved for reference.

These files are intentionally **not** part of the canonical CMake build/test workflow:

- `legacy/monolith/tokamakfusion.cpp`
- `legacy/monolith/include/tokamak_engine.h`
- `legacy/monolith/src/tokamak_engine.cpp`
- `legacy/viewer/src/viewer_main.cpp`
- `legacy/viewer/src/render/*.cpp`
- `legacy/viewer/src/ui/telemetry_panels.cpp`
- `legacy/shared/src/snapshot_stream.cpp`
- `legacy/shared/include/*`

Use the modular runtime in `src/main.cpp` + `include/tokamak/*` via CMake targets (`tokamak_core`, `tokamakfusion`, `tokamak_tests`).

If legacy code must be compiled for archival debugging, include both header roots:
- `-Ilegacy/monolith/include`
- `-Ilegacy/shared/include`
