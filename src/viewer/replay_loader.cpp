#include "tokamak/viewer/replay_loader.hpp"

#include <algorithm>
#include <filesystem>
#include <regex>
#include <utility>

namespace tokamak::viewer {
namespace {

bool ParseStepFromSnapshotFilename(const std::filesystem::path& path, int* outStep) {
    const std::string filename = path.filename().string();
    const std::regex re(R"(particles_step_(\d+)\.csv)");
    std::smatch match;
    if (!std::regex_match(filename, match, re)) {
        return false;
    }

    try {
        *outStep = std::stoi(match[1].str());
        return true;
    } catch (...) {
        return false;
    }
}

}  // namespace

bool ReplayLoader::OpenFromManifest(const std::filesystem::path& manifestPath) {
    ReplayManifest parsedManifest;
    std::string parseError;
    if (!ParseManifestV2File(manifestPath, &parsedManifest, &parseError)) {
        SetError(parseError);
        return false;
    }

    return InitializeFromManifest(parsedManifest);
}

bool ReplayLoader::OpenFromRunDirectory(const std::filesystem::path& runDirectory) {
    if (!std::filesystem::exists(runDirectory)) {
        SetError("Run directory does not exist: " + runDirectory.string());
        return false;
    }
    if (!std::filesystem::is_directory(runDirectory)) {
        SetError("Provided run path is not a directory: " + runDirectory.string());
        return false;
    }

    return OpenFromManifest(runDirectory / "manifest_v2.json");
}

void ReplayLoader::Clear() {
    manifest_ = ReplayManifest();
    runConfig_ = ReplayRunConfig();
    frameEntries_.clear();
    orderedSteps_.clear();
    stepToOrderedIndex_.clear();
    summaryRows_.clear();
    summaryRowByStep_.clear();
    cache_.clear();
    cacheLruOrder_.clear();
    lastError_.clear();
}

bool ReplayLoader::InitializeFromManifest(const ReplayManifest& manifest) {
    Clear();
    manifest_ = manifest;

    const std::filesystem::path runConfigPath = manifest_.runDirectory / manifest_.files.runConfigJson;
    std::string runConfigError;
    if (!ParseRunConfigV2File(runConfigPath, &runConfig_, &runConfigError)) {
        SetError(runConfigError);
        return false;
    }

    frameEntries_.reserve(manifest_.files.particleSnapshotCsvFiles.size());
    for (const std::string& relativePath : manifest_.files.particleSnapshotCsvFiles) {
        const std::filesystem::path snapshotPath = manifest_.runDirectory / relativePath;
        if (!std::filesystem::exists(snapshotPath)) {
            SetError("Snapshot file listed in manifest does not exist: " + snapshotPath.string());
            return false;
        }

        int step = -1;
        if (!ParseStepFromSnapshotFilename(snapshotPath, &step)) {
            SetError("Unable to parse step index from snapshot filename: " + snapshotPath.filename().string());
            return false;
        }

        frameEntries_.push_back(FrameEntry{step, snapshotPath});
    }

    std::sort(frameEntries_.begin(), frameEntries_.end(), [](const FrameEntry& a, const FrameEntry& b) {
        return a.step < b.step;
    });

    orderedSteps_.reserve(frameEntries_.size());
    for (std::size_t i = 0; i < frameEntries_.size(); ++i) {
        orderedSteps_.push_back(frameEntries_[i].step);
        stepToOrderedIndex_[frameEntries_[i].step] = i;
    }

    if (!LoadSummary()) {
        return false;
    }

    return true;
}

bool ReplayLoader::LoadSummary() {
    const std::filesystem::path summaryPath = manifest_.runDirectory / manifest_.files.summaryCsv;
    std::string error;
    if (!ParseSummaryCsv(summaryPath, &summaryRows_, &error)) {
        SetError(error);
        return false;
    }

    summaryRowByStep_.clear();
    for (std::size_t i = 0; i < summaryRows_.size(); ++i) {
        summaryRowByStep_[summaryRows_[i].step] = i;
    }

    return true;
}

bool ReplayLoader::LoadFrameUncached(std::size_t orderedIndex, ReplayFrame* outFrame) {
    if (orderedIndex >= frameEntries_.size()) {
        SetError("Frame index out of range: " + std::to_string(orderedIndex));
        return false;
    }

    const std::filesystem::path& path = frameEntries_[orderedIndex].snapshotPath;
    std::string error;
    if (!ParseParticleSnapshotCsv(path, outFrame, &error)) {
        SetError(error);
        return false;
    }

    return true;
}

void ReplayLoader::TouchCacheEntry(std::size_t orderedIndex) {
    auto it = std::find(cacheLruOrder_.begin(), cacheLruOrder_.end(), orderedIndex);
    if (it != cacheLruOrder_.end()) {
        cacheLruOrder_.erase(it);
    }
    cacheLruOrder_.insert(cacheLruOrder_.begin(), orderedIndex);
}

void ReplayLoader::InsertCacheEntry(std::size_t orderedIndex, const ReplayFrame& frame) {
    for (CachedFrame& cached : cache_) {
        if (cached.orderedIndex == orderedIndex) {
            cached.frame = frame;
            TouchCacheEntry(orderedIndex);
            return;
        }
    }

    if (cache_.size() >= cacheCapacity_ && !cacheLruOrder_.empty()) {
        const std::size_t evictIndex = cacheLruOrder_.back();
        cacheLruOrder_.pop_back();

        cache_.erase(
            std::remove_if(cache_.begin(), cache_.end(), [&](const CachedFrame& cached) {
                return cached.orderedIndex == evictIndex;
            }),
            cache_.end());
    }

    cache_.push_back(CachedFrame{orderedIndex, frame});
    TouchCacheEntry(orderedIndex);
}

bool ReplayLoader::LoadFrameByOrderedIndex(std::size_t orderedIndex, ReplayFrame* outFrame) {
    if (outFrame == nullptr) {
        SetError("LoadFrameByOrderedIndex: outFrame is null");
        return false;
    }

    for (const CachedFrame& cached : cache_) {
        if (cached.orderedIndex == orderedIndex) {
            *outFrame = cached.frame;
            TouchCacheEntry(orderedIndex);
            return true;
        }
    }

    ReplayFrame parsed;
    if (!LoadFrameUncached(orderedIndex, &parsed)) {
        return false;
    }

    InsertCacheEntry(orderedIndex, parsed);
    *outFrame = std::move(parsed);
    return true;
}

bool ReplayLoader::LoadFrameByStep(int step, ReplayFrame* outFrame) {
    const auto it = stepToOrderedIndex_.find(step);
    if (it == stepToOrderedIndex_.end()) {
        SetError("No frame exists for step: " + std::to_string(step));
        return false;
    }
    return LoadFrameByOrderedIndex(it->second, outFrame);
}

const ReplaySummaryPoint* ReplayLoader::SummaryForStep(int step) const {
    const auto it = summaryRowByStep_.find(step);
    if (it == summaryRowByStep_.end()) {
        return nullptr;
    }
    return &summaryRows_[it->second];
}

void ReplayLoader::SetError(const std::string& message) {
    lastError_ = message;
}

}  // namespace tokamak::viewer
