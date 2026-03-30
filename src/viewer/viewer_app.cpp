#include "tokamak/viewer/viewer_app.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <utility>

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

#include "tokamak/viewer/camera.hpp"
#include "tokamak/viewer/gl_renderer.hpp"

namespace tokamak::viewer {
namespace {

struct WindowInputState {
    double scrollDeltaY = 0.0;
};

void ScrollCallback(GLFWwindow* window, double /*xoffset*/, double yoffset) {
    WindowInputState* input = static_cast<WindowInputState*>(glfwGetWindowUserPointer(window));
    if (input != nullptr) {
        input->scrollDeltaY += yoffset;
    }
}

bool ParseFloatValue(const char* text, float* outValue) {
    char* end = nullptr;
    const float value = std::strtof(text, &end);
    if (end == text || *end != '\0' || !std::isfinite(value)) {
        return false;
    }
    *outValue = value;
    return true;
}

bool ParseIntValue(const char* text, int* outValue) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value < std::numeric_limits<int>::min() ||
        value > std::numeric_limits<int>::max()) {
        return false;
    }
    *outValue = static_cast<int>(value);
    return true;
}

std::size_t ClampToSizeT(uint64_t value) {
    const uint64_t maxValue = static_cast<uint64_t>(std::numeric_limits<std::size_t>::max());
    return static_cast<std::size_t>(std::min(value, maxValue));
}

constexpr double kMinDisplayFrameDuration_s = 1.0 / 24.0;

float Clamp01(float value) {
    return std::max(0.0f, std::min(value, 1.0f));
}

double Clamp01Double(double value) {
    return std::max(0.0, std::min(value, 1.0));
}

double SmoothStep01(double value) {
    const double clamped = Clamp01Double(value);
    return clamped * clamped * (3.0 - (2.0 * clamped));
}

struct PowerOnVisualState {
    float startupProgress = 1.0f;
    float fusionGateProgress = 1.0f;
    float ignitionIntensity = 0.0f;
    uint64_t recentFusionEvents = 0;
    std::string label = "steady-state";
};

PowerOnVisualState ComputePowerOnVisualState(
    double displayedTime_s,
    const ReplayRunConfig& runConfig,
    const ReplaySummaryPoint* summary,
    const ReplaySummaryPoint* windowStartSummary) {
    PowerOnVisualState state;

    if (runConfig.startupRampDuration_s > 0.0) {
        state.startupProgress = static_cast<float>(
            SmoothStep01(displayedTime_s / runConfig.startupRampDuration_s));
    }

    const double fusionGateStart_s =
        std::max(0.0, runConfig.startupRampDuration_s) + std::max(0.0, runConfig.fusionStartDelay_s);
    if (displayedTime_s < fusionGateStart_s) {
        state.fusionGateProgress = 0.0f;
    } else if (runConfig.fusionRampDuration_s > 0.0) {
        state.fusionGateProgress = static_cast<float>(
            SmoothStep01((displayedTime_s - fusionGateStart_s) / runConfig.fusionRampDuration_s));
    } else {
        state.fusionGateProgress = 1.0f;
    }

    if (summary != nullptr) {
        uint64_t windowStartTotal = 0;
        if (windowStartSummary != nullptr &&
            windowStartSummary->fusionEventsTotal <= summary->fusionEventsTotal) {
            windowStartTotal = windowStartSummary->fusionEventsTotal;
        }
        state.recentFusionEvents = summary->fusionEventsTotal - windowStartTotal;

        const float recentFusionSignal = Clamp01(
            static_cast<float>(std::log10(1.0 + static_cast<double>(state.recentFusionEvents)) / 2.0));
        const float cumulativeFusionSignal = Clamp01(
            static_cast<float>(std::log10(1.0 + static_cast<double>(summary->fusionEventsTotal)) / 2.4));
        const float energySignal = Clamp01(static_cast<float>((summary->avgEnergy_keV - 90.0) / 140.0));
        const float fusionEvidence =
            std::max((0.65f * recentFusionSignal) + (0.35f * cumulativeFusionSignal), 0.20f * energySignal);
        state.ignitionIntensity = state.fusionGateProgress * fusionEvidence;
    }

    if (state.startupProgress < 0.995f) {
        state.label = "powering magnetics + beam";
    } else if (state.fusionGateProgress < 0.05f) {
        state.label = "thermalizing core";
    } else if (state.recentFusionEvents == 0) {
        state.label = (state.fusionGateProgress < 0.995f) ? "fusion conditions forming" : "threshold not yet sustained";
    } else if (state.ignitionIntensity < 0.35f) {
        state.label = "fusion onset";
    } else if (state.ignitionIntensity < 0.75f) {
        state.label = "burn building";
    } else {
        state.label = "ignition established";
    }

    return state;
}

}  // namespace

ViewerApp::ViewerApp(ViewerCliOptions options)
    : options_(std::move(options)) {}

bool ParseViewerCliArgs(int argc, char** argv, ViewerCliOptions* outOptions, std::string* errorOut) {
    if (outOptions == nullptr) {
        if (errorOut != nullptr) {
            *errorOut = "Internal error: outOptions is null";
        }
        return false;
    }

    ViewerCliOptions parsed;

    for (int i = 1; i < argc; ++i) {
        const std::string_view arg(argv[i]);
        auto needValue = [&](const char* optionName) -> const char* {
            if (i + 1 >= argc) {
                if (errorOut != nullptr) {
                    *errorOut = std::string("Missing value for ") + optionName;
                }
                return nullptr;
            }
            return argv[++i];
        };

        if (arg == "--manifest") {
            const char* value = needValue("--manifest");
            if (value == nullptr) {
                return false;
            }
            parsed.manifestPath = value;
            continue;
        }

        if (arg == "--run-dir") {
            const char* value = needValue("--run-dir");
            if (value == nullptr) {
                return false;
            }
            parsed.runDirectory = value;
            continue;
        }

        if (arg == "--point-size") {
            const char* value = needValue("--point-size");
            if (value == nullptr) {
                return false;
            }
            float parsedValue = 0.0f;
            if (!ParseFloatValue(value, &parsedValue) || parsedValue <= 0.0f) {
                if (errorOut != nullptr) {
                    *errorOut = std::string("Invalid --point-size: ") + value;
                }
                return false;
            }
            parsed.pointSizePixels = parsedValue;
            continue;
        }

        if (arg == "--playback-rate") {
            const char* value = needValue("--playback-rate");
            if (value == nullptr) {
                return false;
            }
            float parsedValue = 0.0f;
            if (!ParseFloatValue(value, &parsedValue) || parsedValue <= 0.0f) {
                if (errorOut != nullptr) {
                    *errorOut = std::string("Invalid --playback-rate: ") + value;
                }
                return false;
            }
            parsed.playbackRate = parsedValue;
            continue;
        }

        if (arg == "--start-step") {
            const char* value = needValue("--start-step");
            if (value == nullptr) {
                return false;
            }
            int parsedValue = 0;
            if (!ParseIntValue(value, &parsedValue) || parsedValue < 0) {
                if (errorOut != nullptr) {
                    *errorOut = std::string("Invalid --start-step: ") + value;
                }
                return false;
            }
            parsed.startStep = parsedValue;
            continue;
        }

        if (arg == "--help" || arg == "-h") {
            if (errorOut != nullptr) {
                *errorOut = "help";
            }
            return false;
        }

        if (errorOut != nullptr) {
            *errorOut = std::string("Unknown option: ") + std::string(arg);
        }
        return false;
    }

    const bool hasManifest = !parsed.manifestPath.empty();
    const bool hasRunDir = !parsed.runDirectory.empty();
    if (hasManifest == hasRunDir) {
        if (errorOut != nullptr) {
            *errorOut = "Specify exactly one of --manifest <path> or --run-dir <path>";
        }
        return false;
    }

    *outOptions = parsed;
    return true;
}

void PrintViewerUsage(const char* argv0) {
    std::cout << "Usage: " << argv0 << " [options]\n"
              << "  --manifest <path/to/manifest_v2.json>\n"
              << "  --run-dir <path/to/run_directory>\n"
              << "  --point-size <float>\n"
              << "  --playback-rate <float>\n"
              << "  --start-step <int>\n"
              << "  --help\n";
}

int ViewerApp::Run() {
    const bool opened = !options_.manifestPath.empty()
        ? loader_.OpenFromManifest(options_.manifestPath)
        : loader_.OpenFromRunDirectory(options_.runDirectory);

    if (!opened) {
        std::cerr << "Failed to open replay: " << loader_.LastError() << "\n";
        return 1;
    }

    if (!loader_.HasData()) {
        std::cerr << "No replay frames are available in manifest.\n";
        return 1;
    }

    const auto& orderedSteps = loader_.OrderedSteps();
    std::size_t currentFrameIndex = 0;
    if (options_.startStep >= 0) {
        const auto it = std::lower_bound(orderedSteps.begin(), orderedSteps.end(), options_.startStep);
        if (it == orderedSteps.end() || *it != options_.startStep) {
            std::cerr << "Requested --start-step not found: " << options_.startStep << "\n";
            return 1;
        }
        currentFrameIndex = static_cast<std::size_t>(std::distance(orderedSteps.begin(), it));
    }

    ReplayFrame currentFrame;
    if (!loader_.LoadFrameByOrderedIndex(currentFrameIndex, &currentFrame)) {
        std::cerr << "Failed to load initial frame: " << loader_.LastError() << "\n";
        return 1;
    }

    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        return 1;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE);
#endif

    GLFWwindow* window = glfwCreateWindow(1600, 900, "Tokamak Viewer (Replay)", nullptr, nullptr);
    if (window == nullptr) {
        glfwTerminate();
        std::cerr << "Failed to create GLFW window\n";
        return 1;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    if (!gladLoadGLLoader(reinterpret_cast<GLADloadproc>(glfwGetProcAddress))) {
        std::cerr << "Failed to initialize GLAD\n";
        glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");

    WindowInputState inputState;
    glfwSetWindowUserPointer(window, &inputState);
    glfwSetScrollCallback(window, ScrollCallback);

    const ReplayRunConfig& runConfig = loader_.RunConfig();
    std::size_t initialParticleCapacity = std::max<std::size_t>(
        std::max<std::size_t>(currentFrame.particles.size(), ClampToSizeT(currentFrame.sampledParticles)),
        1);
    if (runConfig.hasMaxParticlesPerSnapshot) {
        initialParticleCapacity = std::max(initialParticleCapacity, ClampToSizeT(runConfig.maxParticlesPerSnapshot));
    }
    if (orderedSteps.size() <= 256 && initialParticleCapacity <= 12000) {
        loader_.SetCacheCapacity(orderedSteps.size());
        for (std::size_t i = 0; i < orderedSteps.size(); ++i) {
            if (!loader_.PrefetchFrameByOrderedIndex(i)) {
                std::cerr << "Frame prefetch error: " << loader_.LastError() << "\n";
                break;
            }
        }
    }

    OrbitCamera camera;
    GlRenderer renderer;
    std::string rendererError;
    if (!renderer.Initialize(initialParticleCapacity, &rendererError)) {
        std::cerr << "Failed to initialize renderer: " << rendererError << "\n";
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();
        glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }

    if (runConfig.hasTokamakGeometry) {
        renderer.SetTorusGeometry(runConfig.majorRadius_m, runConfig.minorRadius_m);
    } else {
        renderer.SetTorusGeometry(2.0f, 0.5f);
        std::cerr << "Warning: tokamak geometry not found in run_config_v2.json, using default torus geometry.\n";
    }

    bool paused = false;
    float playbackRate = options_.playbackRate;
    float pointSize = options_.pointSizePixels;
    double displayFrameProgress_s = 0.0;
    ReplayFrame nextFrame;
    bool hasNextFrame = false;

    auto loadAdjacentFrame = [&](std::size_t baseIndex) -> bool {
        if (baseIndex + 1 >= orderedSteps.size()) {
            hasNextFrame = false;
            nextFrame = ReplayFrame();
            return true;
        }
        if (!loader_.LoadFrameByOrderedIndex(baseIndex + 1, &nextFrame)) {
            std::cerr << "Frame load error: " << loader_.LastError() << "\n";
            hasNextFrame = false;
            return false;
        }
        hasNextFrame = true;
        return true;
    };

    auto prefetchUpcomingFrames = [&](std::size_t baseIndex) {
        constexpr std::size_t kPrefetchDepth = 12;
        for (std::size_t offset = 2; offset < (2 + kPrefetchDepth); ++offset) {
            const std::size_t prefetchIndex = baseIndex + offset;
            if (prefetchIndex >= orderedSteps.size()) {
                break;
            }
            if (!loader_.PrefetchFrameByOrderedIndex(prefetchIndex)) {
                std::cerr << "Frame prefetch error: " << loader_.LastError() << "\n";
                break;
            }
        }
    };

    auto jumpToFrame = [&](std::size_t orderedIndex, bool pausePlayback) -> bool {
        ReplayFrame loadedFrame;
        if (!loader_.LoadFrameByOrderedIndex(orderedIndex, &loadedFrame)) {
            std::cerr << "Frame load error: " << loader_.LastError() << "\n";
            return false;
        }
        currentFrameIndex = orderedIndex;
        currentFrame = std::move(loadedFrame);
        displayFrameProgress_s = 0.0;
        if (!loadAdjacentFrame(currentFrameIndex)) {
            return false;
        }
        renderer.UploadFrame(currentFrame, hasNextFrame ? &nextFrame : nullptr);
        prefetchUpcomingFrames(currentFrameIndex);
        paused = pausePlayback || !hasNextFrame;
        return true;
    };

    if (!jumpToFrame(currentFrameIndex, false)) {
        renderer.Shutdown();
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();
        glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }

    auto lastFrameTime = std::chrono::steady_clock::now();

    bool leftDownPrev = false;

    while (!glfwWindowShouldClose(window)) {
        const auto frameStart = std::chrono::steady_clock::now();
        const double deltaSeconds = std::chrono::duration<double>(frameStart - lastFrameTime).count();
        lastFrameTime = frameStart;

        glfwPollEvents();

        int frameBufferW = 0;
        int frameBufferH = 0;
        glfwGetFramebufferSize(window, &frameBufferW, &frameBufferH);
        glViewport(0, 0, frameBufferW, frameBufferH);
        camera.SetViewport(frameBufferW, frameBufferH);

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        const ImGuiIO& io = ImGui::GetIO();
        double mouseX = 0.0;
        double mouseY = 0.0;
        glfwGetCursorPos(window, &mouseX, &mouseY);
        const bool leftDown = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;

        if (!io.WantCaptureMouse) {
            if (leftDown && !leftDownPrev) {
                camera.BeginRotate(mouseX, mouseY);
            }
            if (!leftDown && leftDownPrev) {
                camera.EndRotate();
            }
            camera.OnCursorMove(mouseX, mouseY);
            if (inputState.scrollDeltaY != 0.0) {
                camera.OnScroll(inputState.scrollDeltaY);
                inputState.scrollDeltaY = 0.0;
            }
        }
        leftDownPrev = leftDown;

        if (!paused && hasNextFrame) {
            const double physicalFrameSpan_s = std::max(0.0, nextFrame.time_s - currentFrame.time_s);
            const double displayFrameSpan_s = std::max(
                kMinDisplayFrameDuration_s,
                physicalFrameSpan_s) / std::max(0.1, static_cast<double>(playbackRate));
            displayFrameProgress_s = std::min(displayFrameProgress_s + deltaSeconds, displayFrameSpan_s);

            if (displayFrameProgress_s >= displayFrameSpan_s) {
                displayFrameProgress_s = 0.0;
                ++currentFrameIndex;
                currentFrame = std::move(nextFrame);
                if (!loadAdjacentFrame(currentFrameIndex)) {
                    paused = true;
                } else {
                    prefetchUpcomingFrames(currentFrameIndex);
                    renderer.UploadFrame(currentFrame, hasNextFrame ? &nextFrame : nullptr);
                    if (!hasNextFrame) {
                        paused = true;
                    }
                }
            }
        }

        float interpolationAlpha = 0.0f;
        double displayedTime_s = currentFrame.time_s;
        double blendFraction = 0.0;
        if (!paused && hasNextFrame) {
            const double physicalFrameSpan_s = std::max(0.0, nextFrame.time_s - currentFrame.time_s);
            const double displayFrameSpan_s = std::max(
                kMinDisplayFrameDuration_s,
                physicalFrameSpan_s) / std::max(0.1, static_cast<double>(playbackRate));
            if (displayFrameSpan_s > 0.0) {
                blendFraction = std::max(0.0, std::min(displayFrameProgress_s / displayFrameSpan_s, 1.0));
            }
            interpolationAlpha = static_cast<float>(blendFraction);
            displayedTime_s = currentFrame.time_s + (blendFraction * physicalFrameSpan_s);
        } else {
            if (!hasNextFrame) {
                paused = true;
            }
            displayFrameProgress_s = 0.0;
        }

        const ReplaySummaryPoint* summary = loader_.SummaryForStep(currentFrame.step);
        const std::size_t fusionWindowStartIndex =
            (currentFrameIndex > 4) ? (currentFrameIndex - 4) : 0;
        const ReplaySummaryPoint* fusionWindowStartSummary =
            loader_.SummaryForStep(orderedSteps[fusionWindowStartIndex]);
        const PowerOnVisualState powerOnState = ComputePowerOnVisualState(
            displayedTime_s,
            runConfig,
            summary,
            fusionWindowStartSummary);
        const float ignitionIntensity = powerOnState.ignitionIntensity;

        if (ImGui::Begin("Replay Controls")) {
            double sampledWeight = 0.0;
            double sampledWeightedEnergyKeV = 0.0;
            double sampledEnergeticWeight = 0.0;
            double sampledEnergyKeVUnweighted = 0.0;
            double maxSpeed = 0.0;
            std::size_t energeticSamples = 0;
            for (const ReplayParticle& particle : currentFrame.particles) {
                sampledWeight += particle.weight;
                maxSpeed = std::max(maxSpeed, particle.speed_mPerS);
                if (particle.kineticEnergy_keV > 0.0) {
                    sampledEnergyKeVUnweighted += particle.kineticEnergy_keV;
                    ++energeticSamples;
                    if (particle.weight > 0.0) {
                        sampledWeightedEnergyKeV += particle.kineticEnergy_keV * particle.weight;
                        sampledEnergeticWeight += particle.weight;
                    }
                }
            }

            ImGui::TextWrapped(
                "Replaying sampled macro-particles from artifact snapshots. Colors reflect species and energy; "
                "motion is smoothed between snapshots for readability.");
            ImGui::TextWrapped("Drag to orbit, scroll to zoom, and use Playback rate to speed up or slow down the replay.");
            ImGui::Separator();
            ImGui::Text("Run: %s", loader_.Manifest().runId.c_str());
            ImGui::Text("Snapshots: %zu total, showing %zu / %zu", orderedSteps.size(), currentFrameIndex + 1, orderedSteps.size());
            ImGui::Text("Solver step: %d", currentFrame.step);
            ImGui::Text("Replay time: %.6f s", displayedTime_s);
            ImGui::Text("Snapshot time: %.6f s", currentFrame.time_s);
            ImGui::Text("Visible macro-particles: %zu", currentFrame.particles.size());
            ImGui::Text("Particle slots in snapshot: %llu", static_cast<unsigned long long>(currentFrame.totalParticles));
            ImGui::Text("Represented ions in sample: %.3e", sampledWeight);
            ImGui::Text(
                "Average sampled energy: %.3f keV",
                sampledEnergeticWeight > 0.0
                    ? (sampledWeightedEnergyKeV / sampledEnergeticWeight)
                    : (energeticSamples > 0 ? (sampledEnergyKeVUnweighted / static_cast<double>(energeticSamples)) : 0.0));
            ImGui::Text("Max sampled speed: %.3e m/s", maxSpeed);
            ImGui::Text("Blend to next snapshot: %.0f%%", blendFraction * 100.0);
            ImGui::Text("Core state: %s", powerOnState.label.c_str());
            ImGui::Text("Power ramp complete: %.0f%%", powerOnState.startupProgress * 100.0f);
            ImGui::Text("Fusion onset gate: %.0f%%", powerOnState.fusionGateProgress * 100.0f);
            ImGui::Text("Ignition glow: %.0f%%", ignitionIntensity * 100.0f);

            if (summary != nullptr) {
                ImGui::Separator();
                ImGui::Text("Whole-plasma ions: %llu", static_cast<unsigned long long>(summary->totalIons));
                ImGui::Text("Whole-plasma avg energy: %.6f keV", summary->avgEnergy_keV);
                ImGui::Text("Whole-plasma fusion events: %llu", static_cast<unsigned long long>(summary->fusionEventsTotal));
                ImGui::Text(
                    "Recent fusion events (%zu snapshots): %llu",
                    std::min<std::size_t>(5, currentFrameIndex + 1),
                    static_cast<unsigned long long>(powerOnState.recentFusionEvents));
            } else {
                ImGui::TextUnformatted("Summary row not found for this step.");
            }

            if (ImGui::Button(paused ? "Play" : "Pause")) {
                paused = !paused;
            }
            ImGui::SameLine();
            if (ImGui::Button("Restart")) {
                jumpToFrame(0, false);
            }

            ImGui::SliderFloat("Playback rate", &playbackRate, 0.1f, 8.0f, "%.2fx");
            ImGui::SliderFloat("Point size", &pointSize, 1.0f, 8.0f, "%.1f px");

            int frameSlider = static_cast<int>(currentFrameIndex);
            if (ImGui::SliderInt("Frame", &frameSlider, 0, static_cast<int>(orderedSteps.size() - 1))) {
                jumpToFrame(static_cast<std::size_t>(frameSlider), true);
            }

            if (!runConfig.hasTokamakGeometry) {
                ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.2f, 1.0f), "Geometry fallback active (default torus).");
            }
        }
        ImGui::End();

        glClearColor(
            0.03f + (0.06f * ignitionIntensity),
            0.05f + (0.01f * ignitionIntensity),
            0.08f + (0.10f * ignitionIntensity),
            1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        renderer.Draw(camera.ViewProjectionMatrix(), pointSize, interpolationAlpha, ignitionIntensity);

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    renderer.Shutdown();
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}

}  // namespace tokamak::viewer
