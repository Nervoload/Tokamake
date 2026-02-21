# Legacy Quarantine Notice

`legacy/` is preserved for historical reference only.

These files are not part of the canonical build, test, or release pipeline:
- `legacy/monolith/*`
- `legacy/viewer/*`
- `legacy/shared/*`

Canonical runtime and tests live under:
- `include/tokamak/*`
- `src/*`
- `tests/*`

Canonical visualization target is `tokamak_viewer` (replay-first, manifest/CSV v2 based).
Legacy `.tksnap` viewer/replay code remains archived and unreferenced.

No build guarantees are provided for files under `legacy/`.
