#pragma once

#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "tokamak/viewer/replay_manifest.hpp"
#include "tokamak/viewer/replay_snapshot.hpp"

namespace tokamak::viewer {

class ReplayLoader {
public:
    bool OpenFromManifest(const std::filesystem::path& manifestPath);
    bool OpenFromRunDirectory(const std::filesystem::path& runDirectory);
    void Clear();

    std::size_t FrameCount() const { return frameEntries_.size(); }
    bool HasData() const { return !frameEntries_.empty(); }

    const ReplayManifest& Manifest() const { return manifest_; }
    const ReplayRunConfig& RunConfig() const { return runConfig_; }

    const std::vector<int>& OrderedSteps() const { return orderedSteps_; }

    bool LoadFrameByOrderedIndex(std::size_t orderedIndex, ReplayFrame* outFrame);
    bool LoadFrameByStep(int step, ReplayFrame* outFrame);

    const ReplaySummaryPoint* SummaryForStep(int step) const;

    const std::string& LastError() const { return lastError_; }

private:
    struct FrameEntry {
        int step = -1;
        std::filesystem::path snapshotPath;
    };

    bool InitializeFromManifest(const ReplayManifest& manifest);
    bool LoadSummary();
    bool LoadFrameUncached(std::size_t orderedIndex, ReplayFrame* outFrame);

    struct CachedFrame {
        std::size_t orderedIndex = 0;
        ReplayFrame frame;
    };

    void TouchCacheEntry(std::size_t orderedIndex);
    void InsertCacheEntry(std::size_t orderedIndex, const ReplayFrame& frame);

    void SetError(const std::string& message);

    ReplayManifest manifest_;
    ReplayRunConfig runConfig_;
    std::vector<FrameEntry> frameEntries_;
    std::vector<int> orderedSteps_;
    std::unordered_map<int, std::size_t> stepToOrderedIndex_;
    std::vector<ReplaySummaryPoint> summaryRows_;
    std::unordered_map<int, std::size_t> summaryRowByStep_;

    std::vector<CachedFrame> cache_;
    std::vector<std::size_t> cacheLruOrder_;
    std::size_t cacheCapacity_ = 3;

    std::string lastError_;
};

}  // namespace tokamak::viewer
