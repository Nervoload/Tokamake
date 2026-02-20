#include "ui/telemetry_panels.h"

#include <algorithm>
#include <cstdint>

#include "imgui.h"

namespace {

const char* SpeciesName(uint8_t species) {
    switch (species) {
        case 0:
            return "Deuterium";
        case 1:
            return "Tritium";
        case 2:
            return "Helium";
        default:
            return "Dead";
    }
}

bool SpeciesMatchesFilter(uint8_t species, int filter) {
    if (filter < 0) {
        return true;
    }
    return static_cast<int>(species) == filter;
}

}  // namespace

void DrawTelemetryPanels(
    const SimulationSnapshot& snapshot,
    const ViewerMetrics& metrics,
    TelemetryUiState& uiState) {
    if (ImGui::Begin("Telemetry Layers")) {
        ImGui::Checkbox("Key telemetry", &uiState.showKeyTelemetry);
        ImGui::Checkbox("Simulation/performance", &uiState.showPerformance);
        ImGui::Checkbox("Parameter values", &uiState.showParameters);
        ImGui::Separator();
        ImGui::Text("Step %llu | Sim time %.6fs", snapshot.stepIndex, snapshot.simTimeSeconds);
    }
    ImGui::End();

    if (uiState.showKeyTelemetry && ImGui::Begin("Key Telemetry")) {
        ImGui::Text("Total ions: %u", snapshot.telemetry.totalIons);
        ImGui::Text(
            "Species counts | D: %u, T: %u, He: %u",
            snapshot.telemetry.deuteriumCount,
            snapshot.telemetry.tritiumCount,
            snapshot.telemetry.heliumCount);
        ImGui::Text("Avg temp: %.3f keV", snapshot.telemetry.avgTempKeV);
        ImGui::Text("Fusion events: %llu", snapshot.telemetry.fusionEvents);
        ImGui::Text(
            "Sampled %zu / %u particles (stride %u)",
            snapshot.sampledParticles.size(),
            snapshot.totalParticleCount,
            snapshot.sampleStride);
    }
    ImGui::End();

    if (uiState.showPerformance && ImGui::Begin("Simulation / Performance")) {
        ImGui::Text("Sim step: %.3f ms", snapshot.performance.simStepMs);
        ImGui::Text("Snapshot export: %.3f ms", snapshot.performance.exportMs);
        ImGui::Text("Render frame: %.3f ms", metrics.renderMs);
        ImGui::Text("Renderer FPS: %.1f", metrics.fps);
        ImGui::Separator();
        ImGui::Text(
            "Ring buffer | produced: %llu consumed: %llu dropped: %llu pending: %zu",
            metrics.producedSnapshots,
            metrics.consumedSnapshots,
            metrics.droppedSnapshots,
            metrics.pendingSnapshots);
        if (metrics.replayMode) {
            ImGui::Separator();
            ImGui::Text("Replay mode");
            ImGui::Text("Playback rate: %.2fx", metrics.playbackRate);
            ImGui::Text("Playback state: %s", metrics.paused ? "Paused" : "Playing");
        }
    }
    ImGui::End();

    if (uiState.showParameters && ImGui::Begin("Parameter Values (Sampled)")) {
        ImGui::SliderInt("Rows per page", &uiState.rowsPerPage, 8, 256);

        const int maxPage = static_cast<int>(
            snapshot.sampledParticles.empty()
                ? 0
                : (snapshot.sampledParticles.size() - 1) / static_cast<size_t>(uiState.rowsPerPage));
        uiState.page = std::clamp(uiState.page, 0, maxPage);
        ImGui::SliderInt("Page", &uiState.page, 0, maxPage);

        const char* speciesItems[] = {"All", "Deuterium", "Tritium", "Helium"};
        int selected = uiState.speciesFilter + 1;
        if (ImGui::Combo("Species filter", &selected, speciesItems, IM_ARRAYSIZE(speciesItems))) {
            uiState.speciesFilter = selected - 1;
            uiState.page = 0;
        }

        const size_t start = static_cast<size_t>(uiState.page) * static_cast<size_t>(uiState.rowsPerPage);
        const size_t end = std::min(
            snapshot.sampledParticles.size(),
            start + static_cast<size_t>(uiState.rowsPerPage));

        if (ImGui::BeginTable("sample_table", 10, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY, ImVec2(0.0f, 320.0f))) {
            ImGui::TableSetupColumn("i");
            ImGui::TableSetupColumn("species");
            ImGui::TableSetupColumn("x");
            ImGui::TableSetupColumn("y");
            ImGui::TableSetupColumn("z");
            ImGui::TableSetupColumn("vx");
            ImGui::TableSetupColumn("vy");
            ImGui::TableSetupColumn("vz");
            ImGui::TableSetupColumn("m");
            ImGui::TableSetupColumn("q/m");
            ImGui::TableHeadersRow();

            for (size_t i = start; i < end; ++i) {
                const ParticleSample& p = snapshot.sampledParticles[i];
                if (!SpeciesMatchesFilter(p.species, uiState.speciesFilter)) {
                    continue;
                }

                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                ImGui::Text("%zu", i);
                ImGui::TableSetColumnIndex(1);
                ImGui::TextUnformatted(SpeciesName(p.species));
                ImGui::TableSetColumnIndex(2);
                ImGui::Text("%.3f", p.x);
                ImGui::TableSetColumnIndex(3);
                ImGui::Text("%.3f", p.y);
                ImGui::TableSetColumnIndex(4);
                ImGui::Text("%.3f", p.z);
                ImGui::TableSetColumnIndex(5);
                ImGui::Text("%.3e", p.vx);
                ImGui::TableSetColumnIndex(6);
                ImGui::Text("%.3e", p.vy);
                ImGui::TableSetColumnIndex(7);
                ImGui::Text("%.3e", p.vz);
                ImGui::TableSetColumnIndex(8);
                ImGui::Text("%.3e", p.mass);
                ImGui::TableSetColumnIndex(9);
                ImGui::Text("%.3e", p.qOverM);
            }

            ImGui::EndTable();
        }
    }
    ImGui::End();
}
