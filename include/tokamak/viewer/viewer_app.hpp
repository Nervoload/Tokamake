#pragma once

#include <filesystem>
#include <string>

#include "tokamak/viewer/replay_loader.hpp"

namespace tokamak::viewer {

struct ViewerCliOptions {
    std::filesystem::path manifestPath;
    std::filesystem::path runDirectory;
    float pointSizePixels = 3.0f;
    float playbackRate = 1.0f;
    int startStep = -1;
};

class ViewerApp {
public:
    explicit ViewerApp(ViewerCliOptions options);

    int Run();

private:
    ViewerCliOptions options_;
    ReplayLoader loader_;
};

bool ParseViewerCliArgs(int argc, char** argv, ViewerCliOptions* outOptions, std::string* errorOut);
void PrintViewerUsage(const char* argv0);

}  // namespace tokamak::viewer
