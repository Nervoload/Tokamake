#include "tokamak/viewer/replay_analytics.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

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

bool ParseUInt32(const std::string& text, uint32_t* outValue) {
    try {
        *outValue = static_cast<uint32_t>(std::stoul(text));
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

bool ParseBool(const std::string& text, bool* outValue) {
    const std::string trimmed = TrimAscii(text);
    if (trimmed == "true" || trimmed == "1" || trimmed == "TRUE") {
        *outValue = true;
        return true;
    }
    if (trimmed == "false" || trimmed == "0" || trimmed == "FALSE") {
        *outValue = false;
        return true;
    }
    return false;
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

template <typename T>
bool OpenCsv(
    const std::filesystem::path& path,
    std::ifstream* input,
    std::unordered_map<std::string, std::size_t>* headerToIndex,
    std::string* errorOut,
    const char* missingHeaderMessage) {
    input->open(path);
    if (!input->is_open()) {
        if (errorOut != nullptr) {
            *errorOut = "Failed to open CSV: " + path.string();
        }
        return false;
    }

    std::string headerLine;
    if (!std::getline(*input, headerLine)) {
        if (errorOut != nullptr) {
            *errorOut = std::string(missingHeaderMessage) + path.string();
        }
        return false;
    }

    const std::vector<std::string> headerFields = SplitCsvLine(headerLine);
    headerToIndex->clear();
    for (std::size_t i = 0; i < headerFields.size(); ++i) {
        (*headerToIndex)[TrimAscii(headerFields[i])] = i;
    }
    return true;
}

template <typename T>
void SortBinsByStepAndIndex(std::vector<T>* rows) {
    std::sort(rows->begin(), rows->end(), [](const T& a, const T& b) {
        if (a.step != b.step) {
            return a.step < b.step;
        }
        return a.binIndex < b.binIndex;
    });
}

template <typename T>
void BuildBinsByStep(const std::vector<T>& rows, std::unordered_map<int, std::vector<T>>* outByStep) {
    outByStep->clear();
    for (const T& row : rows) {
        (*outByStep)[row.step].push_back(row);
    }
}

template <typename T>
void SortSingleSeries(std::vector<T>* rows) {
    std::sort(rows->begin(), rows->end(), [](const T& a, const T& b) {
        if (a.step != b.step) {
            return a.step < b.step;
        }
        return a.time_s < b.time_s;
    });
}

template <typename T>
void BuildSingleIndex(const std::vector<T>& rows, std::unordered_map<int, std::size_t>* outByStep) {
    outByStep->clear();
    for (std::size_t i = 0; i < rows.size(); ++i) {
        (*outByStep)[rows[i].step] = i;
    }
}

bool EnsureFileExists(
    const std::filesystem::path& path,
    bool required,
    std::string* errorOut) {
    if (path.empty()) {
        if (required && errorOut != nullptr) {
            *errorOut = "Required analytics CSV path is empty";
        }
        return !required;
    }
    if (std::filesystem::exists(path)) {
        return true;
    }
    if (required && errorOut != nullptr) {
        *errorOut = "Analytics CSV listed in manifest does not exist: " + path.string();
    }
    return !required;
}

bool ParseRadialProfiles(
    const std::filesystem::path& path,
    std::vector<ReplayRadialProfileBin>* outRows,
    std::string* errorOut) {
    std::ifstream input;
    std::unordered_map<std::string, std::size_t> headerToIndex;
    if (!OpenCsv<ReplayRadialProfileBin>(
            path,
            &input,
            &headerToIndex,
            errorOut,
            "Radial profile CSV missing header row: ")) {
        return false;
    }

    std::size_t stepCol = 0;
    std::size_t timeCol = 0;
    std::size_t binCol = 0;
    std::size_t rInnerCol = 0;
    std::size_t rOuterCol = 0;
    std::size_t rCenterCol = 0;
    std::size_t ionCountCol = 0;
    std::size_t macroWeightCol = 0;
    std::size_t shellVolumeCol = 0;
    std::size_t densityCol = 0;
    std::size_t avgEnergyCol = 0;
    std::size_t fusionEventsCol = 0;
    std::size_t fusionRateCol = 0;
    std::size_t fusionPlaceholderCol = 0;

    if (!LookupColumn(headerToIndex, "step", &stepCol, errorOut) ||
        !LookupColumn(headerToIndex, "time_s", &timeCol, errorOut) ||
        !LookupColumn(headerToIndex, "bin_index", &binCol, errorOut) ||
        !LookupColumn(headerToIndex, "r_inner_m", &rInnerCol, errorOut) ||
        !LookupColumn(headerToIndex, "r_outer_m", &rOuterCol, errorOut) ||
        !LookupColumn(headerToIndex, "r_center_m", &rCenterCol, errorOut) ||
        !LookupColumn(headerToIndex, "ion_count", &ionCountCol, errorOut) ||
        !LookupColumn(headerToIndex, "macro_weight", &macroWeightCol, errorOut) ||
        !LookupColumn(headerToIndex, "shell_volume_m3", &shellVolumeCol, errorOut) ||
        !LookupColumn(headerToIndex, "density_m3", &densityCol, errorOut) ||
        !LookupColumn(headerToIndex, "avg_ion_energy_kev", &avgEnergyCol, errorOut) ||
        !LookupColumn(headerToIndex, "fusion_events_cumulative", &fusionEventsCol, errorOut) ||
        !LookupColumn(headerToIndex, "fusion_rate_m3_s", &fusionRateCol, errorOut) ||
        !LookupColumn(headerToIndex, "fusion_rate_placeholder", &fusionPlaceholderCol, errorOut)) {
        return false;
    }

    const std::size_t maxRequired = std::max(
        {stepCol, timeCol, binCol, rInnerCol, rOuterCol, rCenterCol, ionCountCol, macroWeightCol, shellVolumeCol,
         densityCol, avgEnergyCol, fusionEventsCol, fusionRateCol, fusionPlaceholderCol});

    outRows->clear();
    std::string line;
    int lineNumber = 1;
    while (std::getline(input, line)) {
        ++lineNumber;
        if (TrimAscii(line).empty()) {
            continue;
        }
        const std::vector<std::string> fields = SplitCsvLine(line);
        if (fields.size() <= maxRequired) {
            if (errorOut != nullptr) {
                *errorOut = "Radial profile CSV row has too few columns at line " + std::to_string(lineNumber);
            }
            return false;
        }

        ReplayRadialProfileBin row;
        bool fusionPlaceholder = false;
        if (!ParseInt(TrimAscii(fields[stepCol]), &row.step) ||
            !ParseDouble(TrimAscii(fields[timeCol]), &row.time_s) ||
            !ParseInt(TrimAscii(fields[binCol]), &row.binIndex) ||
            !ParseDouble(TrimAscii(fields[rInnerCol]), &row.rInner_m) ||
            !ParseDouble(TrimAscii(fields[rOuterCol]), &row.rOuter_m) ||
            !ParseDouble(TrimAscii(fields[rCenterCol]), &row.rCenter_m) ||
            !ParseUInt64(TrimAscii(fields[ionCountCol]), &row.ionCount) ||
            !ParseDouble(TrimAscii(fields[macroWeightCol]), &row.macroWeight) ||
            !ParseDouble(TrimAscii(fields[shellVolumeCol]), &row.shellVolume_m3) ||
            !ParseDouble(TrimAscii(fields[densityCol]), &row.density_m3) ||
            !ParseDouble(TrimAscii(fields[avgEnergyCol]), &row.avgIonEnergy_keV) ||
            !ParseUInt64(TrimAscii(fields[fusionEventsCol]), &row.fusionEventsCumulative) ||
            !ParseDouble(TrimAscii(fields[fusionRateCol]), &row.fusionRate_m3_s) ||
            !ParseBool(TrimAscii(fields[fusionPlaceholderCol]), &fusionPlaceholder)) {
            if (errorOut != nullptr) {
                *errorOut = "Radial profile CSV parse error at line " + std::to_string(lineNumber);
            }
            return false;
        }
        row.fusionRatePlaceholder = fusionPlaceholder;
        outRows->push_back(row);
    }

    SortBinsByStepAndIndex(outRows);
    return true;
}

bool ParseMagneticFieldDiagnostics(
    const std::filesystem::path& path,
    std::vector<ReplayMagneticFieldBin>* outRows,
    std::vector<ReplayMagneticFieldStepSummary>* outSeries,
    std::string* errorOut) {
    std::ifstream input;
    std::unordered_map<std::string, std::size_t> headerToIndex;
    if (!OpenCsv<ReplayMagneticFieldBin>(
            path,
            &input,
            &headerToIndex,
            errorOut,
            "Magnetic field diagnostics CSV missing header row: ")) {
        return false;
    }

    std::size_t stepCol = 0;
    std::size_t timeCol = 0;
    std::size_t binCol = 0;
    std::size_t rInnerCol = 0;
    std::size_t rOuterCol = 0;
    std::size_t rCenterCol = 0;
    std::size_t meanBCol = 0;
    std::size_t sampleCountCol = 0;
    std::size_t maxBCol = 0;
    std::size_t recommendedDtCol = 0;
    std::size_t profileKindCol = 0;

    if (!LookupColumn(headerToIndex, "step", &stepCol, errorOut) ||
        !LookupColumn(headerToIndex, "time_s", &timeCol, errorOut) ||
        !LookupColumn(headerToIndex, "bin_index", &binCol, errorOut) ||
        !LookupColumn(headerToIndex, "r_inner_m", &rInnerCol, errorOut) ||
        !LookupColumn(headerToIndex, "r_outer_m", &rOuterCol, errorOut) ||
        !LookupColumn(headerToIndex, "r_center_m", &rCenterCol, errorOut) ||
        !LookupColumn(headerToIndex, "mean_b_t", &meanBCol, errorOut) ||
        !LookupColumn(headerToIndex, "sample_count", &sampleCountCol, errorOut) ||
        !LookupColumn(headerToIndex, "step_max_b_t", &maxBCol, errorOut) ||
        !LookupColumn(headerToIndex, "recommended_dt_s", &recommendedDtCol, errorOut) ||
        !LookupColumn(headerToIndex, "profile_kind", &profileKindCol, errorOut)) {
        return false;
    }

    const std::size_t maxRequired = std::max(
        {stepCol, timeCol, binCol, rInnerCol, rOuterCol, rCenterCol, meanBCol, sampleCountCol, maxBCol,
         recommendedDtCol, profileKindCol});

    outRows->clear();
    outSeries->clear();
    std::unordered_map<int, std::size_t> summaryByStep;

    std::string line;
    int lineNumber = 1;
    while (std::getline(input, line)) {
        ++lineNumber;
        if (TrimAscii(line).empty()) {
            continue;
        }
        const std::vector<std::string> fields = SplitCsvLine(line);
        if (fields.size() <= maxRequired) {
            if (errorOut != nullptr) {
                *errorOut = "Magnetic field diagnostics row has too few columns at line " + std::to_string(lineNumber);
            }
            return false;
        }

        ReplayMagneticFieldBin row;
        if (!ParseInt(TrimAscii(fields[stepCol]), &row.step) ||
            !ParseDouble(TrimAscii(fields[timeCol]), &row.time_s) ||
            !ParseInt(TrimAscii(fields[binCol]), &row.binIndex) ||
            !ParseDouble(TrimAscii(fields[rInnerCol]), &row.rInner_m) ||
            !ParseDouble(TrimAscii(fields[rOuterCol]), &row.rOuter_m) ||
            !ParseDouble(TrimAscii(fields[rCenterCol]), &row.rCenter_m) ||
            !ParseDouble(TrimAscii(fields[meanBCol]), &row.meanField_T) ||
            !ParseUInt64(TrimAscii(fields[sampleCountCol]), &row.sampleCount) ||
            !ParseDouble(TrimAscii(fields[maxBCol]), &row.stepMaxField_T) ||
            !ParseDouble(TrimAscii(fields[recommendedDtCol]), &row.recommendedDt_s)) {
            if (errorOut != nullptr) {
                *errorOut = "Magnetic field diagnostics parse error at line " + std::to_string(lineNumber);
            }
            return false;
        }
        row.profileKind = TrimAscii(fields[profileKindCol]);
        outRows->push_back(row);

        if (summaryByStep.find(row.step) == summaryByStep.end()) {
            summaryByStep[row.step] = outSeries->size();
            outSeries->push_back(
                ReplayMagneticFieldStepSummary{
                    row.step,
                    row.time_s,
                    row.stepMaxField_T,
                    row.recommendedDt_s,
                    row.profileKind});
        }
    }

    SortBinsByStepAndIndex(outRows);
    SortSingleSeries(outSeries);
    return true;
}

bool ParseElectrostaticDiagnostics(
    const std::filesystem::path& path,
    std::vector<ReplayElectrostaticPoint>* outRows,
    std::string* errorOut) {
    std::ifstream input;
    std::unordered_map<std::string, std::size_t> headerToIndex;
    if (!OpenCsv<ReplayElectrostaticPoint>(
            path,
            &input,
            &headerToIndex,
            errorOut,
            "Electrostatic diagnostics CSV missing header row: ")) {
        return false;
    }

    std::size_t stepCol = 0;
    std::size_t timeCol = 0;
    std::size_t modeCol = 0;
    std::size_t boundaryCol = 0;
    std::size_t assignmentCol = 0;
    std::size_t maxECol = 0;
    std::size_t meanECol = 0;
    std::size_t iterationsCol = 0;
    std::size_t convergedCol = 0;
    std::size_t residualCol = 0;

    if (!LookupColumn(headerToIndex, "step", &stepCol, errorOut) ||
        !LookupColumn(headerToIndex, "time_s", &timeCol, errorOut) ||
        !LookupColumn(headerToIndex, "electric_field_mode", &modeCol, errorOut) ||
        !LookupColumn(headerToIndex, "boundary_condition", &boundaryCol, errorOut) ||
        !LookupColumn(headerToIndex, "charge_assignment", &assignmentCol, errorOut) ||
        !LookupColumn(headerToIndex, "max_electric_field_v_per_m", &maxECol, errorOut) ||
        !LookupColumn(headerToIndex, "mean_electric_field_v_per_m", &meanECol, errorOut) ||
        !LookupColumn(headerToIndex, "solver_iterations", &iterationsCol, errorOut) ||
        !LookupColumn(headerToIndex, "solver_converged", &convergedCol, errorOut) ||
        !LookupColumn(headerToIndex, "residual_l2", &residualCol, errorOut)) {
        return false;
    }

    const std::size_t maxRequired = std::max(
        {stepCol, timeCol, modeCol, boundaryCol, assignmentCol, maxECol, meanECol, iterationsCol, convergedCol,
         residualCol});

    outRows->clear();
    std::string line;
    int lineNumber = 1;
    while (std::getline(input, line)) {
        ++lineNumber;
        if (TrimAscii(line).empty()) {
            continue;
        }
        const std::vector<std::string> fields = SplitCsvLine(line);
        if (fields.size() <= maxRequired) {
            if (errorOut != nullptr) {
                *errorOut = "Electrostatic diagnostics row has too few columns at line " + std::to_string(lineNumber);
            }
            return false;
        }

        ReplayElectrostaticPoint row;
        bool solverConverged = false;
        if (!ParseInt(TrimAscii(fields[stepCol]), &row.step) ||
            !ParseDouble(TrimAscii(fields[timeCol]), &row.time_s) ||
            !ParseDouble(TrimAscii(fields[maxECol]), &row.maxElectricField_VPerM) ||
            !ParseDouble(TrimAscii(fields[meanECol]), &row.meanElectricField_VPerM) ||
            !ParseUInt32(TrimAscii(fields[iterationsCol]), &row.solverIterations) ||
            !ParseBool(TrimAscii(fields[convergedCol]), &solverConverged) ||
            !ParseDouble(TrimAscii(fields[residualCol]), &row.residualL2)) {
            if (errorOut != nullptr) {
                *errorOut = "Electrostatic diagnostics parse error at line " + std::to_string(lineNumber);
            }
            return false;
        }
        row.solverConverged = solverConverged;
        row.electricFieldMode = TrimAscii(fields[modeCol]);
        row.boundaryCondition = TrimAscii(fields[boundaryCol]);
        row.chargeAssignment = TrimAscii(fields[assignmentCol]);
        outRows->push_back(std::move(row));
    }

    SortSingleSeries(outRows);
    return true;
}

bool ParseSpeedHistogram(
    const std::filesystem::path& path,
    std::vector<ReplaySpeedHistogramBin>* outRows,
    std::string* errorOut) {
    std::ifstream input;
    std::unordered_map<std::string, std::size_t> headerToIndex;
    if (!OpenCsv<ReplaySpeedHistogramBin>(
            path,
            &input,
            &headerToIndex,
            errorOut,
            "Speed histogram CSV missing header row: ")) {
        return false;
    }

    std::size_t stepCol = 0;
    std::size_t timeCol = 0;
    std::size_t binCol = 0;
    std::size_t lowerCol = 0;
    std::size_t upperCol = 0;
    std::size_t countCol = 0;
    std::size_t totalSamplesCol = 0;

    if (!LookupColumn(headerToIndex, "step", &stepCol, errorOut) ||
        !LookupColumn(headerToIndex, "time_s", &timeCol, errorOut) ||
        !LookupColumn(headerToIndex, "bin_index", &binCol, errorOut) ||
        !LookupColumn(headerToIndex, "speed_min_m_per_s", &lowerCol, errorOut) ||
        !LookupColumn(headerToIndex, "speed_max_m_per_s", &upperCol, errorOut) ||
        !LookupColumn(headerToIndex, "count", &countCol, errorOut) ||
        !LookupColumn(headerToIndex, "total_samples", &totalSamplesCol, errorOut)) {
        return false;
    }

    const std::size_t maxRequired = std::max({stepCol, timeCol, binCol, lowerCol, upperCol, countCol, totalSamplesCol});

    outRows->clear();
    std::string line;
    int lineNumber = 1;
    while (std::getline(input, line)) {
        ++lineNumber;
        if (TrimAscii(line).empty()) {
            continue;
        }
        const std::vector<std::string> fields = SplitCsvLine(line);
        if (fields.size() <= maxRequired) {
            if (errorOut != nullptr) {
                *errorOut = "Speed histogram row has too few columns at line " + std::to_string(lineNumber);
            }
            return false;
        }

        ReplaySpeedHistogramBin row;
        if (!ParseInt(TrimAscii(fields[stepCol]), &row.step) ||
            !ParseDouble(TrimAscii(fields[timeCol]), &row.time_s) ||
            !ParseInt(TrimAscii(fields[binCol]), &row.binIndex) ||
            !ParseDouble(TrimAscii(fields[lowerCol]), &row.speedMin_mPerS) ||
            !ParseDouble(TrimAscii(fields[upperCol]), &row.speedMax_mPerS) ||
            !ParseUInt64(TrimAscii(fields[countCol]), &row.count) ||
            !ParseUInt64(TrimAscii(fields[totalSamplesCol]), &row.totalSamples)) {
            if (errorOut != nullptr) {
                *errorOut = "Speed histogram parse error at line " + std::to_string(lineNumber);
            }
            return false;
        }
        outRows->push_back(row);
    }

    SortBinsByStepAndIndex(outRows);
    return true;
}

bool ParsePitchHistogram(
    const std::filesystem::path& path,
    std::vector<ReplayPitchHistogramBin>* outRows,
    std::string* errorOut) {
    std::ifstream input;
    std::unordered_map<std::string, std::size_t> headerToIndex;
    if (!OpenCsv<ReplayPitchHistogramBin>(
            path,
            &input,
            &headerToIndex,
            errorOut,
            "Pitch histogram CSV missing header row: ")) {
        return false;
    }

    std::size_t stepCol = 0;
    std::size_t timeCol = 0;
    std::size_t binCol = 0;
    std::size_t lowerCol = 0;
    std::size_t upperCol = 0;
    std::size_t countCol = 0;
    std::size_t totalSamplesCol = 0;
    std::size_t invalidSamplesCol = 0;

    if (!LookupColumn(headerToIndex, "step", &stepCol, errorOut) ||
        !LookupColumn(headerToIndex, "time_s", &timeCol, errorOut) ||
        !LookupColumn(headerToIndex, "bin_index", &binCol, errorOut) ||
        !LookupColumn(headerToIndex, "pitch_min_deg", &lowerCol, errorOut) ||
        !LookupColumn(headerToIndex, "pitch_max_deg", &upperCol, errorOut) ||
        !LookupColumn(headerToIndex, "count", &countCol, errorOut) ||
        !LookupColumn(headerToIndex, "total_samples", &totalSamplesCol, errorOut) ||
        !LookupColumn(headerToIndex, "invalid_samples", &invalidSamplesCol, errorOut)) {
        return false;
    }

    const std::size_t maxRequired =
        std::max({stepCol, timeCol, binCol, lowerCol, upperCol, countCol, totalSamplesCol, invalidSamplesCol});

    outRows->clear();
    std::string line;
    int lineNumber = 1;
    while (std::getline(input, line)) {
        ++lineNumber;
        if (TrimAscii(line).empty()) {
            continue;
        }
        const std::vector<std::string> fields = SplitCsvLine(line);
        if (fields.size() <= maxRequired) {
            if (errorOut != nullptr) {
                *errorOut = "Pitch histogram row has too few columns at line " + std::to_string(lineNumber);
            }
            return false;
        }

        ReplayPitchHistogramBin row;
        if (!ParseInt(TrimAscii(fields[stepCol]), &row.step) ||
            !ParseDouble(TrimAscii(fields[timeCol]), &row.time_s) ||
            !ParseInt(TrimAscii(fields[binCol]), &row.binIndex) ||
            !ParseDouble(TrimAscii(fields[lowerCol]), &row.pitchMin_deg) ||
            !ParseDouble(TrimAscii(fields[upperCol]), &row.pitchMax_deg) ||
            !ParseUInt64(TrimAscii(fields[countCol]), &row.count) ||
            !ParseUInt64(TrimAscii(fields[totalSamplesCol]), &row.totalSamples) ||
            !ParseUInt64(TrimAscii(fields[invalidSamplesCol]), &row.invalidSamples)) {
            if (errorOut != nullptr) {
                *errorOut = "Pitch histogram parse error at line " + std::to_string(lineNumber);
            }
            return false;
        }
        outRows->push_back(row);
    }

    SortBinsByStepAndIndex(outRows);
    return true;
}

bool ParseSolverResiduals(
    const std::filesystem::path& path,
    std::vector<ReplaySolverResidualPoint>* outRows,
    std::string* errorOut) {
    std::ifstream input;
    std::unordered_map<std::string, std::size_t> headerToIndex;
    if (!OpenCsv<ReplaySolverResidualPoint>(
            path,
            &input,
            &headerToIndex,
            errorOut,
            "Solver residual CSV missing header row: ")) {
        return false;
    }

    std::size_t stepCol = 0;
    std::size_t timeCol = 0;
    std::size_t availableCol = 0;
    std::size_t residualCol = 0;
    std::size_t solverNameCol = 0;
    std::size_t statusCol = 0;
    std::size_t iterationsCol = 0;
    std::size_t convergedCol = 0;
    std::size_t toleranceCol = 0;
    std::size_t noteCol = 0;

    if (!LookupColumn(headerToIndex, "step", &stepCol, errorOut) ||
        !LookupColumn(headerToIndex, "time_s", &timeCol, errorOut) ||
        !LookupColumn(headerToIndex, "residual_available", &availableCol, errorOut) ||
        !LookupColumn(headerToIndex, "residual_l2", &residualCol, errorOut) ||
        !LookupColumn(headerToIndex, "solver_name", &solverNameCol, errorOut) ||
        !LookupColumn(headerToIndex, "status", &statusCol, errorOut) ||
        !LookupColumn(headerToIndex, "iterations", &iterationsCol, errorOut) ||
        !LookupColumn(headerToIndex, "converged", &convergedCol, errorOut) ||
        !LookupColumn(headerToIndex, "tolerance", &toleranceCol, errorOut) ||
        !LookupColumn(headerToIndex, "note", &noteCol, errorOut)) {
        return false;
    }

    const std::size_t maxRequired = std::max(
        {stepCol, timeCol, availableCol, residualCol, solverNameCol, statusCol, iterationsCol, convergedCol,
         toleranceCol, noteCol});

    outRows->clear();
    std::string line;
    int lineNumber = 1;
    while (std::getline(input, line)) {
        ++lineNumber;
        if (TrimAscii(line).empty()) {
            continue;
        }
        const std::vector<std::string> fields = SplitCsvLine(line);
        if (fields.size() <= maxRequired) {
            if (errorOut != nullptr) {
                *errorOut = "Solver residual row has too few columns at line " + std::to_string(lineNumber);
            }
            return false;
        }

        ReplaySolverResidualPoint row;
        bool residualAvailable = false;
        bool converged = false;
        if (!ParseInt(TrimAscii(fields[stepCol]), &row.step) ||
            !ParseDouble(TrimAscii(fields[timeCol]), &row.time_s) ||
            !ParseBool(TrimAscii(fields[availableCol]), &residualAvailable) ||
            !ParseDouble(TrimAscii(fields[residualCol]), &row.residualL2) ||
            !ParseUInt32(TrimAscii(fields[iterationsCol]), &row.iterations) ||
            !ParseBool(TrimAscii(fields[convergedCol]), &converged) ||
            !ParseDouble(TrimAscii(fields[toleranceCol]), &row.tolerance)) {
            if (errorOut != nullptr) {
                *errorOut = "Solver residual parse error at line " + std::to_string(lineNumber);
            }
            return false;
        }
        row.residualAvailable = residualAvailable;
        row.converged = converged;
        row.solverName = TrimAscii(fields[solverNameCol]);
        row.status = TrimAscii(fields[statusCol]);
        row.note = fields[noteCol];
        outRows->push_back(std::move(row));
    }

    SortSingleSeries(outRows);
    return true;
}

bool ParseFusionReactivity(
    const std::filesystem::path& path,
    std::vector<ReplayFusionReactivityPoint>* outRows,
    std::string* errorOut) {
    std::ifstream input;
    std::unordered_map<std::string, std::size_t> headerToIndex;
    if (!OpenCsv<ReplayFusionReactivityPoint>(
            path,
            &input,
            &headerToIndex,
            errorOut,
            "Fusion reactivity CSV missing header row: ")) {
        return false;
    }

    std::size_t stepCol = 0;
    std::size_t timeCol = 0;
    std::size_t modelCol = 0;
    std::size_t sigmaScaleCol = 0;
    std::size_t probabilityClampCol = 0;
    std::size_t minEnergyCol = 0;
    std::size_t attemptsStepCol = 0;
    std::size_t acceptedStepCol = 0;
    std::size_t weightAttemptedCol = 0;
    std::size_t weightAcceptedCol = 0;
    std::size_t fuelDCol = 0;
    std::size_t fuelTCol = 0;
    std::size_t ashHeCol = 0;
    std::size_t avgSigmaCol = 0;
    std::size_t avgProbabilityCol = 0;
    std::size_t avgRelativeSpeedCol = 0;
    std::size_t attemptsTotalCol = 0;
    std::size_t acceptedTotalCol = 0;

    if (!LookupColumn(headerToIndex, "step", &stepCol, errorOut) ||
        !LookupColumn(headerToIndex, "time_s", &timeCol, errorOut) ||
        !LookupColumn(headerToIndex, "reactivity_model", &modelCol, errorOut) ||
        !LookupColumn(headerToIndex, "cross_section_scale", &sigmaScaleCol, errorOut) ||
        !LookupColumn(headerToIndex, "probability_clamp", &probabilityClampCol, errorOut) ||
        !LookupColumn(headerToIndex, "min_energy_kev", &minEnergyCol, errorOut) ||
        !LookupColumn(headerToIndex, "fusion_attempts_step", &attemptsStepCol, errorOut) ||
        !LookupColumn(headerToIndex, "fusion_accepted_step", &acceptedStepCol, errorOut) ||
        !LookupColumn(headerToIndex, "fusion_weight_attempted_step", &weightAttemptedCol, errorOut) ||
        !LookupColumn(headerToIndex, "fusion_weight_accepted_step", &weightAcceptedCol, errorOut) ||
        !LookupColumn(headerToIndex, "fuel_weight_consumed_d_step", &fuelDCol, errorOut) ||
        !LookupColumn(headerToIndex, "fuel_weight_consumed_t_step", &fuelTCol, errorOut) ||
        !LookupColumn(headerToIndex, "ash_weight_produced_he_step", &ashHeCol, errorOut) ||
        !LookupColumn(headerToIndex, "avg_sigma_m2_step", &avgSigmaCol, errorOut) ||
        !LookupColumn(headerToIndex, "avg_probability_step", &avgProbabilityCol, errorOut) ||
        !LookupColumn(headerToIndex, "avg_relative_speed_m_per_s_step", &avgRelativeSpeedCol, errorOut) ||
        !LookupColumn(headerToIndex, "fusion_attempts_total", &attemptsTotalCol, errorOut) ||
        !LookupColumn(headerToIndex, "fusion_accepted_total", &acceptedTotalCol, errorOut)) {
        return false;
    }

    const std::size_t maxRequired = std::max(
        {stepCol, timeCol, modelCol, sigmaScaleCol, probabilityClampCol, minEnergyCol, attemptsStepCol,
         acceptedStepCol, weightAttemptedCol, weightAcceptedCol, fuelDCol, fuelTCol, ashHeCol, avgSigmaCol,
         avgProbabilityCol, avgRelativeSpeedCol, attemptsTotalCol, acceptedTotalCol});

    outRows->clear();
    std::string line;
    int lineNumber = 1;
    while (std::getline(input, line)) {
        ++lineNumber;
        if (TrimAscii(line).empty()) {
            continue;
        }
        const std::vector<std::string> fields = SplitCsvLine(line);
        if (fields.size() <= maxRequired) {
            if (errorOut != nullptr) {
                *errorOut = "Fusion reactivity row has too few columns at line " + std::to_string(lineNumber);
            }
            return false;
        }

        ReplayFusionReactivityPoint row;
        if (!ParseInt(TrimAscii(fields[stepCol]), &row.step) ||
            !ParseDouble(TrimAscii(fields[timeCol]), &row.time_s) ||
            !ParseDouble(TrimAscii(fields[sigmaScaleCol]), &row.crossSectionScale) ||
            !ParseDouble(TrimAscii(fields[probabilityClampCol]), &row.probabilityClamp) ||
            !ParseDouble(TrimAscii(fields[minEnergyCol]), &row.minEnergy_keV) ||
            !ParseUInt64(TrimAscii(fields[attemptsStepCol]), &row.fusionAttemptsStep) ||
            !ParseUInt64(TrimAscii(fields[acceptedStepCol]), &row.fusionAcceptedStep) ||
            !ParseDouble(TrimAscii(fields[weightAttemptedCol]), &row.fusionWeightAttemptedStep) ||
            !ParseDouble(TrimAscii(fields[weightAcceptedCol]), &row.fusionWeightAcceptedStep) ||
            !ParseDouble(TrimAscii(fields[fuelDCol]), &row.fuelWeightConsumedDStep) ||
            !ParseDouble(TrimAscii(fields[fuelTCol]), &row.fuelWeightConsumedTStep) ||
            !ParseDouble(TrimAscii(fields[ashHeCol]), &row.ashWeightProducedHeStep) ||
            !ParseDouble(TrimAscii(fields[avgSigmaCol]), &row.avgSigma_m2Step) ||
            !ParseDouble(TrimAscii(fields[avgProbabilityCol]), &row.avgProbabilityStep) ||
            !ParseDouble(TrimAscii(fields[avgRelativeSpeedCol]), &row.avgRelativeSpeed_mPerSStep) ||
            !ParseUInt64(TrimAscii(fields[attemptsTotalCol]), &row.fusionAttemptsTotal) ||
            !ParseUInt64(TrimAscii(fields[acceptedTotalCol]), &row.fusionAcceptedTotal)) {
            if (errorOut != nullptr) {
                *errorOut = "Fusion reactivity parse error at line " + std::to_string(lineNumber);
            }
            return false;
        }
        row.reactivityModel = TrimAscii(fields[modelCol]);
        outRows->push_back(std::move(row));
    }

    SortSingleSeries(outRows);
    return true;
}

bool ParseWallInteraction(
    const std::filesystem::path& path,
    std::vector<ReplayWallInteractionPoint>* outRows,
    std::string* errorOut) {
    std::ifstream input;
    std::unordered_map<std::string, std::size_t> headerToIndex;
    if (!OpenCsv<ReplayWallInteractionPoint>(
            path,
            &input,
            &headerToIndex,
            errorOut,
            "Wall interaction CSV missing header row: ")) {
        return false;
    }

    std::size_t stepCol = 0;
    std::size_t timeCol = 0;
    std::size_t wallModeCol = 0;
    std::size_t recycleCol = 0;
    std::size_t hitStepCol = 0;
    std::size_t impactStepCol = 0;
    std::size_t lossStepCol = 0;
    std::size_t hitTotalCol = 0;
    std::size_t impactTotalCol = 0;
    std::size_t lossTotalCol = 0;

    if (!LookupColumn(headerToIndex, "step", &stepCol, errorOut) ||
        !LookupColumn(headerToIndex, "time_s", &timeCol, errorOut) ||
        !LookupColumn(headerToIndex, "wall_mode", &wallModeCol, errorOut) ||
        !LookupColumn(headerToIndex, "recycle_fraction", &recycleCol, errorOut) ||
        !LookupColumn(headerToIndex, "wall_hit_count_step", &hitStepCol, errorOut) ||
        !LookupColumn(headerToIndex, "wall_impact_energy_j_step", &impactStepCol, errorOut) ||
        !LookupColumn(headerToIndex, "wall_loss_weight_step", &lossStepCol, errorOut) ||
        !LookupColumn(headerToIndex, "wall_hit_count_total", &hitTotalCol, errorOut) ||
        !LookupColumn(headerToIndex, "wall_impact_energy_j_total", &impactTotalCol, errorOut) ||
        !LookupColumn(headerToIndex, "wall_loss_weight_total", &lossTotalCol, errorOut)) {
        return false;
    }

    const std::size_t maxRequired = std::max(
        {stepCol, timeCol, wallModeCol, recycleCol, hitStepCol, impactStepCol, lossStepCol, hitTotalCol,
         impactTotalCol, lossTotalCol});

    outRows->clear();
    std::string line;
    int lineNumber = 1;
    while (std::getline(input, line)) {
        ++lineNumber;
        if (TrimAscii(line).empty()) {
            continue;
        }
        const std::vector<std::string> fields = SplitCsvLine(line);
        if (fields.size() <= maxRequired) {
            if (errorOut != nullptr) {
                *errorOut = "Wall interaction row has too few columns at line " + std::to_string(lineNumber);
            }
            return false;
        }

        ReplayWallInteractionPoint row;
        if (!ParseInt(TrimAscii(fields[stepCol]), &row.step) ||
            !ParseDouble(TrimAscii(fields[timeCol]), &row.time_s) ||
            !ParseDouble(TrimAscii(fields[recycleCol]), &row.recycleFraction) ||
            !ParseUInt64(TrimAscii(fields[hitStepCol]), &row.wallHitCountStep) ||
            !ParseDouble(TrimAscii(fields[impactStepCol]), &row.wallImpactEnergy_JStep) ||
            !ParseDouble(TrimAscii(fields[lossStepCol]), &row.wallLossWeightStep) ||
            !ParseUInt64(TrimAscii(fields[hitTotalCol]), &row.wallHitCountTotal) ||
            !ParseDouble(TrimAscii(fields[impactTotalCol]), &row.wallImpactEnergy_JTotal) ||
            !ParseDouble(TrimAscii(fields[lossTotalCol]), &row.wallLossWeightTotal)) {
            if (errorOut != nullptr) {
                *errorOut = "Wall interaction parse error at line " + std::to_string(lineNumber);
            }
            return false;
        }
        row.wallMode = TrimAscii(fields[wallModeCol]);
        outRows->push_back(std::move(row));
    }

    SortSingleSeries(outRows);
    return true;
}

}  // namespace

bool ReplayAnalytics::LoadFromManifest(const ReplayManifest& manifest, std::string* errorOut) {
    Clear();

    const std::filesystem::path radialPath = manifest.runDirectory / manifest.files.radialProfilesCsv;
    const std::filesystem::path magneticPath = manifest.runDirectory / manifest.files.magneticFieldDiagnosticsCsv;
    const std::filesystem::path electrostaticPath = manifest.runDirectory / manifest.files.electrostaticDiagnosticsCsv;
    const std::filesystem::path speedHistogramPath = manifest.runDirectory / manifest.files.speedHistogramCsv;
    const std::filesystem::path pitchHistogramPath = manifest.runDirectory / manifest.files.pitchHistogramCsv;
    const std::filesystem::path solverResidualPath = manifest.runDirectory / manifest.files.solverResidualCsv;
    const std::filesystem::path fusionReactivityPath = manifest.files.fusionReactivityDiagnosticsCsv.empty()
        ? std::filesystem::path()
        : (manifest.runDirectory / manifest.files.fusionReactivityDiagnosticsCsv);
    const std::filesystem::path wallInteractionPath = manifest.files.wallInteractionBridgeCsv.empty()
        ? std::filesystem::path()
        : (manifest.runDirectory / manifest.files.wallInteractionBridgeCsv);

    if (!EnsureFileExists(radialPath, true, errorOut) ||
        !EnsureFileExists(magneticPath, true, errorOut) ||
        !EnsureFileExists(electrostaticPath, true, errorOut) ||
        !EnsureFileExists(speedHistogramPath, true, errorOut) ||
        !EnsureFileExists(pitchHistogramPath, true, errorOut) ||
        !EnsureFileExists(solverResidualPath, true, errorOut) ||
        !EnsureFileExists(fusionReactivityPath, false, errorOut) ||
        !EnsureFileExists(wallInteractionPath, false, errorOut)) {
        return false;
    }

    if (!ParseRadialProfiles(radialPath, &radialProfiles_, errorOut) ||
        !ParseMagneticFieldDiagnostics(magneticPath, &magneticFieldBins_, &magneticFieldSeries_, errorOut) ||
        !ParseElectrostaticDiagnostics(electrostaticPath, &electrostaticSeries_, errorOut) ||
        !ParseSpeedHistogram(speedHistogramPath, &speedHistogramBins_, errorOut) ||
        !ParsePitchHistogram(pitchHistogramPath, &pitchHistogramBins_, errorOut) ||
        !ParseSolverResiduals(solverResidualPath, &solverResidualSeries_, errorOut)) {
        Clear();
        return false;
    }

    if (!fusionReactivityPath.empty() &&
        !ParseFusionReactivity(fusionReactivityPath, &fusionReactivitySeries_, errorOut)) {
        Clear();
        return false;
    }
    if (!wallInteractionPath.empty() &&
        !ParseWallInteraction(wallInteractionPath, &wallInteractionSeries_, errorOut)) {
        Clear();
        return false;
    }

    BuildBinsByStep(radialProfiles_, &radialProfilesByStep_);
    BuildBinsByStep(magneticFieldBins_, &magneticFieldBinsByStep_);
    BuildSingleIndex(magneticFieldSeries_, &magneticFieldSeriesByStep_);
    BuildSingleIndex(electrostaticSeries_, &electrostaticSeriesByStep_);
    BuildBinsByStep(speedHistogramBins_, &speedHistogramBinsByStep_);
    BuildBinsByStep(pitchHistogramBins_, &pitchHistogramBinsByStep_);
    BuildSingleIndex(solverResidualSeries_, &solverResidualSeriesByStep_);
    BuildSingleIndex(fusionReactivitySeries_, &fusionReactivitySeriesByStep_);
    BuildSingleIndex(wallInteractionSeries_, &wallInteractionSeriesByStep_);

    return true;
}

void ReplayAnalytics::Clear() {
    radialProfiles_.clear();
    radialProfilesByStep_.clear();
    magneticFieldBins_.clear();
    magneticFieldBinsByStep_.clear();
    magneticFieldSeries_.clear();
    magneticFieldSeriesByStep_.clear();
    electrostaticSeries_.clear();
    electrostaticSeriesByStep_.clear();
    speedHistogramBins_.clear();
    speedHistogramBinsByStep_.clear();
    pitchHistogramBins_.clear();
    pitchHistogramBinsByStep_.clear();
    solverResidualSeries_.clear();
    solverResidualSeriesByStep_.clear();
    fusionReactivitySeries_.clear();
    fusionReactivitySeriesByStep_.clear();
    wallInteractionSeries_.clear();
    wallInteractionSeriesByStep_.clear();
}

bool ReplayAnalytics::HasAnyData() const {
    return !radialProfiles_.empty() ||
           !magneticFieldSeries_.empty() ||
           !electrostaticSeries_.empty() ||
           !speedHistogramBins_.empty() ||
           !pitchHistogramBins_.empty() ||
           !solverResidualSeries_.empty() ||
           !fusionReactivitySeries_.empty() ||
           !wallInteractionSeries_.empty();
}

const std::vector<ReplayRadialProfileBin>* ReplayAnalytics::RadialProfileForStep(int step) const {
    return LookupBins(radialProfilesByStep_, step);
}

const std::vector<ReplayMagneticFieldBin>* ReplayAnalytics::MagneticFieldBinsForStep(int step) const {
    return LookupBins(magneticFieldBinsByStep_, step);
}

const ReplayElectrostaticPoint* ReplayAnalytics::ElectrostaticForStep(int step) const {
    return LookupPoint(electrostaticSeriesByStep_, electrostaticSeries_, step);
}

const std::vector<ReplaySpeedHistogramBin>* ReplayAnalytics::SpeedHistogramForStep(int step) const {
    return LookupBins(speedHistogramBinsByStep_, step);
}

const std::vector<ReplayPitchHistogramBin>* ReplayAnalytics::PitchHistogramForStep(int step) const {
    return LookupBins(pitchHistogramBinsByStep_, step);
}

const ReplaySolverResidualPoint* ReplayAnalytics::SolverResidualForStep(int step) const {
    return LookupPoint(solverResidualSeriesByStep_, solverResidualSeries_, step);
}

const ReplayFusionReactivityPoint* ReplayAnalytics::FusionReactivityForStep(int step) const {
    return LookupPoint(fusionReactivitySeriesByStep_, fusionReactivitySeries_, step);
}

const ReplayWallInteractionPoint* ReplayAnalytics::WallInteractionForStep(int step) const {
    return LookupPoint(wallInteractionSeriesByStep_, wallInteractionSeries_, step);
}

}  // namespace tokamak::viewer
