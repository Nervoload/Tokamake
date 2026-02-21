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

    std::size_t currentFrameIndex = 0;
    if (options_.startStep >= 0) {
        const auto& steps = loader_.OrderedSteps();
        const auto it = std::lower_bound(steps.begin(), steps.end(), options_.startStep);
        if (it == steps.end() || *it != options_.startStep) {
            std::cerr << "Requested --start-step not found: " << options_.startStep << "\n";
            return 1;
        }
        currentFrameIndex = static_cast<std::size_t>(std::distance(steps.begin(), it));
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

    OrbitCamera camera;
    GlRenderer renderer;
    std::string rendererError;
    if (!renderer.Initialize(std::max<std::size_t>(currentFrame.particles.size(), 1), &rendererError)) {
        std::cerr << "Failed to initialize renderer: " << rendererError << "\n";
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();
        glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }

    const ReplayRunConfig& runConfig = loader_.RunConfig();
    if (runConfig.hasTokamakGeometry) {
        renderer.SetTorusGeometry(runConfig.majorRadius_m, runConfig.minorRadius_m);
    } else {
        renderer.SetTorusGeometry(2.0f, 0.5f);
        std::cerr << "Warning: tokamak geometry not found in run_config_v2.json, using default torus geometry.\n";
    }
    renderer.UploadFrame(currentFrame);

    bool paused = false;
    float playbackRate = options_.playbackRate;
    float pointSize = options_.pointSizePixels;
    double playbackAccumulator = 0.0;

    const auto appStart = std::chrono::steady_clock::now();
    auto lastFrameTime = appStart;

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

        const auto& orderedSteps = loader_.OrderedSteps();
        if (!paused && currentFrameIndex + 1 < orderedSteps.size()) {
            playbackAccumulator += deltaSeconds * static_cast<double>(playbackRate);
            while (playbackAccumulator >= (1.0 / 30.0) && currentFrameIndex + 1 < orderedSteps.size()) {
                playbackAccumulator -= (1.0 / 30.0);
                ++currentFrameIndex;
                if (!loader_.LoadFrameByOrderedIndex(currentFrameIndex, &currentFrame)) {
                    paused = true;
                    std::cerr << "Frame load error: " << loader_.LastError() << "\n";
                    break;
                }
                renderer.UploadFrame(currentFrame);
            }
        }

        if (ImGui::Begin("Replay Controls")) {
            ImGui::Text("Run: %s", loader_.Manifest().runId.c_str());
            ImGui::Text("Frame %zu / %zu", currentFrameIndex + 1, orderedSteps.size());
            ImGui::Text("Step: %d", currentFrame.step);
            ImGui::Text("Time: %.6f s", currentFrame.time_s);
            ImGui::Text("Sampled particles: %zu", currentFrame.particles.size());
            ImGui::Text("Total particles (reported): %llu", static_cast<unsigned long long>(currentFrame.totalParticles));

            const ReplaySummaryPoint* summary = loader_.SummaryForStep(currentFrame.step);
            if (summary != nullptr) {
                ImGui::Text("Summary total ions: %llu", static_cast<unsigned long long>(summary->totalIons));
                ImGui::Text("Summary avg energy: %.6f keV", summary->avgEnergy_keV);
                ImGui::Text("Summary fusion events: %llu", static_cast<unsigned long long>(summary->fusionEventsTotal));
            } else {
                ImGui::TextUnformatted("Summary row not found for this step.");
            }

            if (ImGui::Button(paused ? "Play" : "Pause")) {
                paused = !paused;
            }
            ImGui::SameLine();
            if (ImGui::Button("Restart")) {
                paused = false;
                currentFrameIndex = 0;
                playbackAccumulator = 0.0;
                if (loader_.LoadFrameByOrderedIndex(currentFrameIndex, &currentFrame)) {
                    renderer.UploadFrame(currentFrame);
                }
            }

            ImGui::SliderFloat("Playback rate", &playbackRate, 0.1f, 8.0f, "%.2fx");
            ImGui::SliderFloat("Point size", &pointSize, 1.0f, 8.0f, "%.1f px");

            int frameSlider = static_cast<int>(currentFrameIndex);
            if (ImGui::SliderInt("Frame", &frameSlider, 0, static_cast<int>(orderedSteps.size() - 1))) {
                currentFrameIndex = static_cast<std::size_t>(frameSlider);
                playbackAccumulator = 0.0;
                paused = true;
                if (loader_.LoadFrameByOrderedIndex(currentFrameIndex, &currentFrame)) {
                    renderer.UploadFrame(currentFrame);
                }
            }

            if (!runConfig.hasTokamakGeometry) {
                ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.2f, 1.0f), "Geometry fallback active (default torus).");
            }
        }
        ImGui::End();

        glClearColor(0.03f, 0.05f, 0.08f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        renderer.Draw(camera.ViewProjectionMatrix(), pointSize);

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
