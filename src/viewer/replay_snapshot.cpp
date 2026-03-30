#include "tokamak/viewer/replay_snapshot.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace tokamak::viewer {
namespace {

std::string TrimAscii(const std::string& text) {
    const std::size_t first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return std::string();
    }
    const std::size_t last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, (last - first) + 1);
}

std::vector<std::string> SplitCsvLine(const std::string& line) {
    std::vector<std::string> fields;
    std::string field;
    bool inQuotes = false;

    for (std::size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '"') {
            if (inQuotes && i + 1 < line.size() && line[i + 1] == '"') {
                field.push_back('"');
                ++i;
            } else {
                inQuotes = !inQuotes;
            }
            continue;
        }

        if (c == ',' && !inQuotes) {
            fields.push_back(field);
            field.clear();
            continue;
        }

        field.push_back(c);
    }

    fields.push_back(field);
    return fields;
}

bool ParseInt(const std::string& text, int* outValue) {
    try {
        *outValue = std::stoi(text);
        return true;
    } catch (...) {
        return false;
    }
}

bool ParseUInt64(const std::string& text, uint64_t* outValue) {
    try {
        *outValue = static_cast<uint64_t>(std::stoull(text));
        return true;
    } catch (...) {
        return false;
    }
}

bool ParseDouble(const std::string& text, double* outValue) {
    try {
        *outValue = std::stod(text);
        return true;
    } catch (...) {
        return false;
    }
}

bool LookupColumn(
    const std::unordered_map<std::string, std::size_t>& headerToIndex,
    const char* name,
    std::size_t* outIndex,
    std::string* errorOut) {
    const auto it = headerToIndex.find(name);
    if (it == headerToIndex.end()) {
        if (errorOut != nullptr) {
            *errorOut = std::string("Missing required CSV column: ") + name;
        }
        return false;
    }
    *outIndex = it->second;
    return true;
}

bool FindOptionalColumn(
    const std::unordered_map<std::string, std::size_t>& headerToIndex,
    const char* name,
    std::size_t* outIndex) {
    const auto it = headerToIndex.find(name);
    if (it == headerToIndex.end()) {
        return false;
    }
    *outIndex = it->second;
    return true;
}

}  // namespace

ReplaySpecies ParseReplaySpeciesName(const std::string& speciesName) {
    if (speciesName == "Deuterium") {
        return ReplaySpecies::Deuterium;
    }
    if (speciesName == "Tritium") {
        return ReplaySpecies::Tritium;
    }
    if (speciesName == "Helium") {
        return ReplaySpecies::Helium;
    }
    return ReplaySpecies::Unknown;
}

bool ParseParticleSnapshotCsv(
    const std::filesystem::path& snapshotCsvPath,
    ReplayFrame* outFrame,
    std::string* errorOut) {
    if (outFrame == nullptr) {
        if (errorOut != nullptr) {
            *errorOut = "ParseParticleSnapshotCsv: outFrame is null";
        }
        return false;
    }

    std::ifstream input(snapshotCsvPath);
    if (!input.is_open()) {
        if (errorOut != nullptr) {
            *errorOut = "Failed to open snapshot CSV: " + snapshotCsvPath.string();
        }
        return false;
    }

    std::string headerLine;
    if (!std::getline(input, headerLine)) {
        if (errorOut != nullptr) {
            *errorOut = "Snapshot CSV missing header row: " + snapshotCsvPath.string();
        }
        return false;
    }

    const std::vector<std::string> headerFields = SplitCsvLine(headerLine);
    std::unordered_map<std::string, std::size_t> headerToIndex;
    for (std::size_t i = 0; i < headerFields.size(); ++i) {
        headerToIndex[TrimAscii(headerFields[i])] = i;
    }

    std::size_t stepColumn = 0;
    std::size_t timeColumn = 0;
    std::size_t totalParticlesColumn = 0;
    std::size_t sampledParticlesColumn = 0;
    std::size_t sampleStrideColumn = 0;
    std::size_t particleIndexColumn = 0;
    std::size_t speciesNameColumn = 0;
    std::size_t xColumn = 0;
    std::size_t yColumn = 0;
    std::size_t zColumn = 0;
    std::size_t vxColumn = 0;
    std::size_t vyColumn = 0;
    std::size_t vzColumn = 0;
    std::size_t weightColumn = 0;
    std::size_t speedColumn = 0;
    std::size_t kineticEnergyColumn = 0;
    std::size_t pitchAngleColumn = 0;

    if (!LookupColumn(headerToIndex, "step", &stepColumn, errorOut) ||
        !LookupColumn(headerToIndex, "time_s", &timeColumn, errorOut) ||
        !LookupColumn(headerToIndex, "total_particles", &totalParticlesColumn, errorOut) ||
        !LookupColumn(headerToIndex, "sampled_particles", &sampledParticlesColumn, errorOut) ||
        !LookupColumn(headerToIndex, "sample_stride", &sampleStrideColumn, errorOut) ||
        !LookupColumn(headerToIndex, "species_name", &speciesNameColumn, errorOut) ||
        !LookupColumn(headerToIndex, "x_m", &xColumn, errorOut) ||
        !LookupColumn(headerToIndex, "y_m", &yColumn, errorOut) ||
        !LookupColumn(headerToIndex, "z_m", &zColumn, errorOut)) {
        return false;
    }

    const bool hasParticleIndex = FindOptionalColumn(headerToIndex, "particle_index", &particleIndexColumn);
    const bool hasVelocity =
        FindOptionalColumn(headerToIndex, "vx_m_per_s", &vxColumn) &&
        FindOptionalColumn(headerToIndex, "vy_m_per_s", &vyColumn) &&
        FindOptionalColumn(headerToIndex, "vz_m_per_s", &vzColumn);
    const bool hasWeight = FindOptionalColumn(headerToIndex, "weight", &weightColumn);
    const bool hasSpeed = FindOptionalColumn(headerToIndex, "speed_m_per_s", &speedColumn);
    const bool hasKineticEnergy = FindOptionalColumn(headerToIndex, "kinetic_energy_kev", &kineticEnergyColumn);
    const bool hasPitchAngle = FindOptionalColumn(headerToIndex, "pitch_angle_deg", &pitchAngleColumn);

    std::size_t maxRequiredColumn = std::max(
        {stepColumn, timeColumn, totalParticlesColumn, sampledParticlesColumn, sampleStrideColumn,
         speciesNameColumn, xColumn, yColumn, zColumn});
    if (hasParticleIndex) {
        maxRequiredColumn = std::max(maxRequiredColumn, particleIndexColumn);
    }
    if (hasVelocity) {
        maxRequiredColumn = std::max(maxRequiredColumn, std::max(vxColumn, std::max(vyColumn, vzColumn)));
    }
    if (hasWeight) {
        maxRequiredColumn = std::max(maxRequiredColumn, weightColumn);
    }
    if (hasSpeed) {
        maxRequiredColumn = std::max(maxRequiredColumn, speedColumn);
    }
    if (hasKineticEnergy) {
        maxRequiredColumn = std::max(maxRequiredColumn, kineticEnergyColumn);
    }
    if (hasPitchAngle) {
        maxRequiredColumn = std::max(maxRequiredColumn, pitchAngleColumn);
    }

    ReplayFrame parsed;
    bool sawAnyRow = false;

    std::string line;
    int lineNumber = 1;
    while (std::getline(input, line)) {
        ++lineNumber;
        if (TrimAscii(line).empty()) {
            continue;
        }

        const std::vector<std::string> fields = SplitCsvLine(line);
        if (fields.size() <= maxRequiredColumn) {
            if (errorOut != nullptr) {
                *errorOut = "Snapshot CSV row has too few columns at line " + std::to_string(lineNumber);
            }
            return false;
        }

        int step = 0;
        double time_s = 0.0;
        uint64_t totalParticles = 0;
        uint64_t sampledParticles = 0;
        uint64_t sampleStride = 0;
        uint64_t particleIndex = 0;
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        double vx = 0.0;
        double vy = 0.0;
        double vz = 0.0;
        double weight = 0.0;
        double speed = 0.0;
        double kineticEnergy = 0.0;
        double pitchAngle = std::numeric_limits<double>::quiet_NaN();

        if (!ParseInt(TrimAscii(fields[stepColumn]), &step) ||
            !ParseDouble(TrimAscii(fields[timeColumn]), &time_s) ||
            !ParseUInt64(TrimAscii(fields[totalParticlesColumn]), &totalParticles) ||
            !ParseUInt64(TrimAscii(fields[sampledParticlesColumn]), &sampledParticles) ||
            !ParseUInt64(TrimAscii(fields[sampleStrideColumn]), &sampleStride) ||
            !ParseDouble(TrimAscii(fields[xColumn]), &x) ||
            !ParseDouble(TrimAscii(fields[yColumn]), &y) ||
            !ParseDouble(TrimAscii(fields[zColumn]), &z)) {
            if (errorOut != nullptr) {
                *errorOut = "Snapshot CSV parse error at line " + std::to_string(lineNumber);
            }
            return false;
        }

        if (hasParticleIndex && !ParseUInt64(TrimAscii(fields[particleIndexColumn]), &particleIndex)) {
            if (errorOut != nullptr) {
                *errorOut = "Snapshot CSV particle_index parse error at line " + std::to_string(lineNumber);
            }
            return false;
        }
        if (hasVelocity &&
            (!ParseDouble(TrimAscii(fields[vxColumn]), &vx) ||
             !ParseDouble(TrimAscii(fields[vyColumn]), &vy) ||
             !ParseDouble(TrimAscii(fields[vzColumn]), &vz))) {
            if (errorOut != nullptr) {
                *errorOut = "Snapshot CSV velocity parse error at line " + std::to_string(lineNumber);
            }
            return false;
        }
        if (hasWeight && !ParseDouble(TrimAscii(fields[weightColumn]), &weight)) {
            if (errorOut != nullptr) {
                *errorOut = "Snapshot CSV weight parse error at line " + std::to_string(lineNumber);
            }
            return false;
        }
        if (hasSpeed && !ParseDouble(TrimAscii(fields[speedColumn]), &speed)) {
            if (errorOut != nullptr) {
                *errorOut = "Snapshot CSV speed parse error at line " + std::to_string(lineNumber);
            }
            return false;
        }
        if (hasKineticEnergy && !ParseDouble(TrimAscii(fields[kineticEnergyColumn]), &kineticEnergy)) {
            if (errorOut != nullptr) {
                *errorOut = "Snapshot CSV kinetic energy parse error at line " + std::to_string(lineNumber);
            }
            return false;
        }
        if (hasPitchAngle && !ParseDouble(TrimAscii(fields[pitchAngleColumn]), &pitchAngle)) {
            if (errorOut != nullptr) {
                *errorOut = "Snapshot CSV pitch angle parse error at line " + std::to_string(lineNumber);
            }
            return false;
        }

        if (!sawAnyRow) {
            parsed.step = step;
            parsed.time_s = time_s;
            parsed.totalParticles = totalParticles;
            parsed.sampledParticles = sampledParticles;
            parsed.sampleStride = std::max<uint64_t>(1u, sampleStride);
            sawAnyRow = true;
        }

        ReplayParticle particle;
        particle.particleIndex = hasParticleIndex ? particleIndex : static_cast<uint64_t>(parsed.particles.size());
        particle.position_m = Vec3(static_cast<float>(x), static_cast<float>(y), static_cast<float>(z));
        particle.velocity_mPerS = Vec3(static_cast<float>(vx), static_cast<float>(vy), static_cast<float>(vz));
        particle.speciesName = TrimAscii(fields[speciesNameColumn]);
        particle.species = ParseReplaySpeciesName(particle.speciesName);
        particle.weight = hasWeight ? weight : 0.0;
        particle.speed_mPerS = hasSpeed ? speed : static_cast<double>(particle.velocity_mPerS.Magnitude());
        particle.kineticEnergy_keV = hasKineticEnergy ? kineticEnergy : 0.0;
        particle.pitchAngle_deg = hasPitchAngle ? pitchAngle : std::numeric_limits<double>::quiet_NaN();
        parsed.particles.push_back(std::move(particle));
    }

    if (!sawAnyRow) {
        if (errorOut != nullptr) {
            *errorOut = "Snapshot CSV contains no particle rows: " + snapshotCsvPath.string();
        }
        return false;
    }

    *outFrame = std::move(parsed);
    return true;
}

bool ParseSummaryCsv(
    const std::filesystem::path& summaryCsvPath,
    std::vector<ReplaySummaryPoint>* outSummary,
    std::string* errorOut) {
    if (outSummary == nullptr) {
        if (errorOut != nullptr) {
            *errorOut = "ParseSummaryCsv: outSummary is null";
        }
        return false;
    }

    std::ifstream input(summaryCsvPath);
    if (!input.is_open()) {
        if (errorOut != nullptr) {
            *errorOut = "Failed to open summary CSV: " + summaryCsvPath.string();
        }
        return false;
    }

    std::string headerLine;
    if (!std::getline(input, headerLine)) {
        if (errorOut != nullptr) {
            *errorOut = "Summary CSV missing header row: " + summaryCsvPath.string();
        }
        return false;
    }

    const std::vector<std::string> headerFields = SplitCsvLine(headerLine);
    std::unordered_map<std::string, std::size_t> headerToIndex;
    for (std::size_t i = 0; i < headerFields.size(); ++i) {
        headerToIndex[TrimAscii(headerFields[i])] = i;
    }

    std::size_t stepColumn = 0;
    std::size_t timeColumn = 0;
    std::size_t totalIonsColumn = 0;
    std::size_t avgEnergyColumn = 0;
    std::size_t fusionColumn = 0;

    if (!LookupColumn(headerToIndex, "step", &stepColumn, errorOut) ||
        !LookupColumn(headerToIndex, "time_s", &timeColumn, errorOut) ||
        !LookupColumn(headerToIndex, "total_ions", &totalIonsColumn, errorOut) ||
        !LookupColumn(headerToIndex, "avg_energy_kev", &avgEnergyColumn, errorOut) ||
        !LookupColumn(headerToIndex, "fusion_events_total", &fusionColumn, errorOut)) {
        return false;
    }

    const std::size_t maxRequiredColumn =
        std::max({stepColumn, timeColumn, totalIonsColumn, avgEnergyColumn, fusionColumn});

    std::vector<ReplaySummaryPoint> parsed;

    std::string line;
    int lineNumber = 1;
    while (std::getline(input, line)) {
        ++lineNumber;
        if (TrimAscii(line).empty()) {
            continue;
        }

        const std::vector<std::string> fields = SplitCsvLine(line);
        if (fields.size() <= maxRequiredColumn) {
            if (errorOut != nullptr) {
                *errorOut = "Summary CSV row has too few columns at line " + std::to_string(lineNumber);
            }
            return false;
        }

        ReplaySummaryPoint point;
        if (!ParseInt(TrimAscii(fields[stepColumn]), &point.step) ||
            !ParseDouble(TrimAscii(fields[timeColumn]), &point.time_s) ||
            !ParseUInt64(TrimAscii(fields[totalIonsColumn]), &point.totalIons) ||
            !ParseDouble(TrimAscii(fields[avgEnergyColumn]), &point.avgEnergy_keV) ||
            !ParseUInt64(TrimAscii(fields[fusionColumn]), &point.fusionEventsTotal)) {
            if (errorOut != nullptr) {
                *errorOut = "Summary CSV parse error at line " + std::to_string(lineNumber);
            }
            return false;
        }

        parsed.push_back(point);
    }

    std::sort(parsed.begin(), parsed.end(), [](const ReplaySummaryPoint& a, const ReplaySummaryPoint& b) {
        return a.step < b.step;
    });

    *outSummary = std::move(parsed);
    return true;
}

}  // namespace tokamak::viewer
