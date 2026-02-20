#pragma once

#include "sim_snapshot.h"

struct ViewerMetrics {
    double renderMs = 0.0;
    double fps = 0.0;
    uint64_t producedSnapshots = 0;
    uint64_t consumedSnapshots = 0;
    uint64_t droppedSnapshots = 0;
    size_t pendingSnapshots = 0;
    bool replayMode = false;
    bool paused = false;
    float playbackRate = 1.0f;
};

struct TelemetryUiState {
    bool showKeyTelemetry = true;
    bool showPerformance = true;
    bool showParameters = true;

    int rowsPerPage = 32;
    int page = 0;
    int speciesFilter = -1;
};

void DrawTelemetryPanels(
    const SimulationSnapshot& snapshot,
    const ViewerMetrics& metrics,
    TelemetryUiState& uiState);
