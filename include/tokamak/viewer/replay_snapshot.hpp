#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "tokamak/types.hpp"

namespace tokamak::viewer {

enum class ReplaySpecies : uint8_t {
    Deuterium = 0,
    Tritium = 1,
    Helium = 2,
    Unknown = 255,
};

struct ReplayParticle {
    Vec3 position_m;
    ReplaySpecies species = ReplaySpecies::Unknown;
    std::string speciesName;
};

struct ReplayFrame {
    int step = 0;
    double time_s = 0.0;
    uint64_t totalParticles = 0;
    uint64_t sampledParticles = 0;
    uint64_t sampleStride = 1;
    std::vector<ReplayParticle> particles;
};

struct ReplaySummaryPoint {
    int step = 0;
    double time_s = 0.0;
    uint64_t totalIons = 0;
    double avgEnergy_keV = 0.0;
    uint64_t fusionEventsTotal = 0;
};

bool ParseParticleSnapshotCsv(
    const std::filesystem::path& snapshotCsvPath,
    ReplayFrame* outFrame,
    std::string* errorOut);

bool ParseSummaryCsv(
    const std::filesystem::path& summaryCsvPath,
    std::vector<ReplaySummaryPoint>* outSummary,
    std::string* errorOut);

ReplaySpecies ParseReplaySpeciesName(const std::string& speciesName);

}  // namespace tokamak::viewer
