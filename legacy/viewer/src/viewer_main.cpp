#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

#include "render/camera.h"
#include "render/particle_renderer.h"
#include "render/scene_renderer.h"
#include "sim_snapshot.h"
#include "snapshot_ring_buffer.h"
#include "snapshot_stream.h"
#include "tokamak_engine.h"
#include "ui/telemetry_panels.h"

namespace {

struct Options {
    bool replayMode = false;
    std::string replayPath;
    std::string recordPath;
    std::string exportCsvPrefix;

    Scenario scenario = Scenario::NBI_IGNITION;
    float timeStep = 1.0e-7f;
    uint32_t maxSampledParticles = 100000;
    uint32_t publishStride = 5;
    size_t ringCapacity = 16;
};

struct WindowInputState {
    double scrollDeltaY = 0.0;
};

struct LiveRuntime {
    explicit LiveRuntime(size_t ringCapacity, size_t reserveParticles)
        : ring(ringCapacity, reserveParticles) {}

    SnapshotRingBuffer ring;
    SnapshotRequest request;
    std::atomic<bool> running{true};

    float timeStep = 1.0e-7f;
    uint32_t publishStride = 1;
    Scenario scenario = Scenario::NBI_IGNITION;
    std::string recordPath;

    std::mutex errorMutex;
    std::string writerError;
};

void PrintUsage() {
    std::cout
        << "Usage: tokamak_viewer [options]\n"
        << "  --replay <file.tksnap>         Run replay mode from snapshot file\n"
        << "  --record <file.tksnap>         Record live snapshots while sim runs\n"
        << "  --export-csv <prefix>          Export replay file to CSV prefix\n"
        << "  --max-sampled <N>              Max sampled particles to visualize (default: 100000)\n"
        << "  --publish-stride <N>           Publish snapshot every N sim steps (default: 5)\n"
        << "  --dt <seconds>                 Simulation step dt (default: 1e-7)\n"
        << "  --scenario <cold|ignition|failure>\n";
}

bool ParseUnsigned(const char* text, uint32_t& out) {
    char* end = nullptr;
    const unsigned long value = std::strtoul(text, &end, 10);
    if (end == text || *end != '\0' || value > std::numeric_limits<uint32_t>::max()) {
        return false;
    }
    out = static_cast<uint32_t>(value);
    return true;
}

bool ParseDouble(const char* text, double& out) {
    char* end = nullptr;
    const double value = std::strtod(text, &end);
    if (end == text || *end != '\0') {
        return false;
    }
    out = value;
    return true;
}

bool ParseOptions(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        auto requireValue = [&](const char* name) -> const char* {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << name << "\n";
                return nullptr;
            }
            return argv[++i];
        };

        if (arg == "--replay") {
            const char* value = requireValue("--replay");
            if (value == nullptr) {
                return false;
            }
            options.replayMode = true;
            options.replayPath = value;
        } else if (arg == "--record") {
            const char* value = requireValue("--record");
            if (value == nullptr) {
                return false;
            }
            options.recordPath = value;
        } else if (arg == "--export-csv") {
            const char* value = requireValue("--export-csv");
            if (value == nullptr) {
                return false;
            }
            options.exportCsvPrefix = value;
        } else if (arg == "--max-sampled") {
            const char* value = requireValue("--max-sampled");
            if (value == nullptr) {
                return false;
            }
            uint32_t parsed = 0;
            if (!ParseUnsigned(value, parsed)) {
                std::cerr << "Invalid --max-sampled value\n";
                return false;
            }
            options.maxSampledParticles = parsed;
        } else if (arg == "--publish-stride") {
            const char* value = requireValue("--publish-stride");
            if (value == nullptr) {
                return false;
            }
            uint32_t parsed = 0;
            if (!ParseUnsigned(value, parsed) || parsed == 0) {
                std::cerr << "Invalid --publish-stride value\n";
                return false;
            }
            options.publishStride = parsed;
        } else if (arg == "--dt") {
            const char* value = requireValue("--dt");
            if (value == nullptr) {
                return false;
            }
            double parsed = 0.0;
            if (!ParseDouble(value, parsed) || parsed <= 0.0) {
                std::cerr << "Invalid --dt value\n";
                return false;
            }
            options.timeStep = static_cast<float>(parsed);
        } else if (arg == "--scenario") {
            const char* value = requireValue("--scenario");
            if (value == nullptr) {
                return false;
            }
            if (std::strcmp(value, "cold") == 0) {
                options.scenario = Scenario::COLD_VACUUM;
            } else if (std::strcmp(value, "ignition") == 0) {
                options.scenario = Scenario::NBI_IGNITION;
            } else if (std::strcmp(value, "failure") == 0) {
                options.scenario = Scenario::MAGNETIC_FAILURE;
            } else {
                std::cerr << "Invalid --scenario value\n";
                return false;
            }
        } else if (arg == "--help" || arg == "-h") {
            PrintUsage();
            return false;
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            return false;
        }
    }

    if (options.replayMode && options.replayPath.empty()) {
        std::cerr << "Replay mode requires --replay <file>\n";
        return false;
    }

    return true;
}

void ScrollCallback(GLFWwindow* window, double /*xoffset*/, double yoffset) {
    WindowInputState* input = static_cast<WindowInputState*>(glfwGetWindowUserPointer(window));
    if (input != nullptr) {
        input->scrollDeltaY += yoffset;
    }
}

bool LoadReplayFrames(const std::string& path, std::vector<SimulationSnapshot>& frames, std::string& error) {
    SnapshotStreamReader reader;
    if (!reader.Open(path)) {
        error = reader.LastError();
        return false;
    }

    SimulationSnapshot frame;
    while (reader.ReadNext(frame)) {
        frames.push_back(frame);
    }

    if (!reader.ReachedEof()) {
        error = reader.LastError().empty() ? "Replay read failed" : reader.LastError();
        return false;
    }

    if (frames.empty()) {
        error = "Replay file has no frames";
        return false;
    }

    return true;
}

void RunLiveSimulation(LiveRuntime* runtime) {
    TokamakEngine engine(runtime->scenario);

    SnapshotStreamWriter writer;
    bool recording = false;
    if (!runtime->recordPath.empty()) {
        recording = writer.Open(runtime->recordPath);
        if (!recording) {
            std::lock_guard<std::mutex> lock(runtime->errorMutex);
            runtime->writerError = writer.LastError();
        }
    }

    SimulationSnapshot snapshot;
    snapshot.Reserve(runtime->request.maxSampledParticles);

    while (runtime->running.load(std::memory_order_relaxed)) {
        engine.Step(runtime->timeStep);

        if ((engine.StepIndex() % runtime->publishStride) != 0) {
            continue;
        }

        engine.ExportSnapshot(snapshot, runtime->request);
        runtime->ring.Publish(snapshot);

        if (recording && !writer.Write(snapshot)) {
            std::lock_guard<std::mutex> lock(runtime->errorMutex);
            runtime->writerError = writer.LastError();
            recording = false;
            writer.Close();
        }
    }

    writer.Close();
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!ParseOptions(argc, argv, options)) {
        PrintUsage();
        return 1;
    }

    if (options.replayMode && !options.exportCsvPrefix.empty()) {
        std::string error;
        if (!ExportSnapshotFileToCsv(options.replayPath, options.exportCsvPrefix, 64, &error)) {
            std::cerr << "CSV export failed: " << error << "\n";
            return 1;
        }
        std::cout << "CSV export complete: " << options.exportCsvPrefix << "_aggregates.csv and _particles.csv\n";
    }

    std::vector<SimulationSnapshot> replayFrames;
    if (options.replayMode) {
        std::string loadError;
        if (!LoadReplayFrames(options.replayPath, replayFrames, loadError)) {
            std::cerr << "Failed to load replay: " << loadError << "\n";
            return 1;
        }
        std::cout << "Loaded " << replayFrames.size() << " replay frames from " << options.replayPath << "\n";
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

    GLFWwindow* window = glfwCreateWindow(1600, 900, "Tokamak Viewer", nullptr, nullptr);
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

    ImGui_ImplGlfw_InitForOpenGL(window, false);
    ImGui_ImplOpenGL3_Init("#version 330");

    WindowInputState inputState;
    glfwSetWindowUserPointer(window, &inputState);
    glfwSetScrollCallback(window, ScrollCallback);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_PROGRAM_POINT_SIZE);

    OrbitCamera camera;
    ParticleRenderer particleRenderer;
    SceneRenderer sceneRenderer;
    TelemetryUiState telemetryUi;

    if (!particleRenderer.Initialize(options.maxSampledParticles)) {
        std::cerr << "Failed to initialize particle renderer\n";
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();
        glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }

    if (!sceneRenderer.Initialize()) {
        std::cerr << "Failed to initialize scene renderer\n";
        particleRenderer.Shutdown();
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();
        glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }

    bool sceneInitialized = false;
    SimulationSnapshot currentSnapshot;
    currentSnapshot.Reserve(options.maxSampledParticles);

    std::thread simulationThread;
    std::unique_ptr<LiveRuntime> liveRuntime;

    if (!options.replayMode) {
        liveRuntime = std::make_unique<LiveRuntime>(options.ringCapacity, options.maxSampledParticles);
        liveRuntime->timeStep = options.timeStep;
        liveRuntime->publishStride = options.publishStride;
        liveRuntime->scenario = options.scenario;
        liveRuntime->recordPath = options.recordPath;
        liveRuntime->request.maxSampledParticles = options.maxSampledParticles;
        liveRuntime->request.includeVelocities = true;
        liveRuntime->request.includeMassAndCharge = true;

        simulationThread = std::thread([&]() {
            RunLiveSimulation(liveRuntime.get());
        });
    } else {
        currentSnapshot = replayFrames.front();
        sceneRenderer.UpdateGeometry(currentSnapshot.config, currentSnapshot.nbi);
        sceneInitialized = true;
        particleRenderer.Upload(currentSnapshot.sampledParticles);
    }

    size_t replayIndex = 0;
    double replayAccumulator = 0.0;
    bool replayPaused = false;
    float replayRate = 1.0f;

    bool leftWasDown = false;
    bool rightWasDown = false;

    const auto appStart = std::chrono::steady_clock::now();
    auto lastFrameTime = appStart;
    double fpsSmoothed = 60.0;
    double lastRenderMs = 0.0;
    bool haveSnapshot = options.replayMode;

    while (!glfwWindowShouldClose(window)) {
        const auto frameStart = std::chrono::steady_clock::now();
        const double deltaSeconds = std::chrono::duration<double>(frameStart - lastFrameTime).count();
        lastFrameTime = frameStart;

        glfwPollEvents();

        int framebufferW = 0;
        int framebufferH = 0;
        glfwGetFramebufferSize(window, &framebufferW, &framebufferH);
        glViewport(0, 0, framebufferW, framebufferH);
        camera.SetViewport(framebufferW, framebufferH);

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGuiIO& io = ImGui::GetIO();

        double mouseX = 0.0;
        double mouseY = 0.0;
        glfwGetCursorPos(window, &mouseX, &mouseY);

        const bool leftDown = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        const bool rightDown = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;

        if (!io.WantCaptureMouse) {
            if (leftDown && !leftWasDown) {
                camera.BeginDrag(0, mouseX, mouseY);
            }
            if (!leftDown && leftWasDown) {
                camera.EndDrag(0);
            }

            if (rightDown && !rightWasDown) {
                camera.BeginDrag(1, mouseX, mouseY);
            }
            if (!rightDown && rightWasDown) {
                camera.EndDrag(1);
            }

            camera.OnCursorMove(mouseX, mouseY);
            if (inputState.scrollDeltaY != 0.0) {
                camera.OnScroll(inputState.scrollDeltaY);
                inputState.scrollDeltaY = 0.0;
            }
        } else {
            inputState.scrollDeltaY = 0.0;
        }

        leftWasDown = leftDown;
        rightWasDown = rightDown;

        if (options.replayMode) {
            if (ImGui::Begin("Replay Controls")) {
                if (ImGui::Button(replayPaused ? "Play" : "Pause")) {
                    replayPaused = !replayPaused;
                }
                ImGui::SameLine();
                if (ImGui::Button("Restart")) {
                    replayIndex = 0;
                    replayAccumulator = 0.0;
                    replayPaused = false;
                }

                ImGui::SliderFloat("Speed", &replayRate, 0.25f, 4.0f, "%.2fx");

                int frameIndexInt = static_cast<int>(replayIndex);
                if (ImGui::SliderInt(
                        "Frame",
                        &frameIndexInt,
                        0,
                        static_cast<int>(replayFrames.size() - 1))) {
                    replayIndex = static_cast<size_t>(frameIndexInt);
                    replayPaused = true;
                    replayAccumulator = 0.0;
                }
            }
            ImGui::End();

            if (!replayPaused && replayIndex + 1 < replayFrames.size()) {
                replayAccumulator += deltaSeconds * static_cast<double>(replayRate);

                while (replayIndex + 1 < replayFrames.size()) {
                    const double stepDt = replayFrames[replayIndex + 1].simTimeSeconds -
                        replayFrames[replayIndex].simTimeSeconds;
                    const double frameDt = (stepDt > 0.0) ? stepDt : (1.0 / 60.0);
                    if (replayAccumulator < frameDt) {
                        break;
                    }
                    replayAccumulator -= frameDt;
                    replayIndex++;
                }

                currentSnapshot = replayFrames[replayIndex];
                particleRenderer.Upload(currentSnapshot.sampledParticles);
                haveSnapshot = true;
            }
        } else if (liveRuntime != nullptr) {
            SimulationSnapshot latest;
            latest.Reserve(options.maxSampledParticles);
            if (liveRuntime->ring.ConsumeLatest(latest)) {
                currentSnapshot = latest;
                particleRenderer.Upload(currentSnapshot.sampledParticles);
                haveSnapshot = true;
                sceneInitialized = false;
            }
        }

        if (haveSnapshot && !sceneInitialized) {
            sceneRenderer.UpdateGeometry(currentSnapshot.config, currentSnapshot.nbi);
            sceneInitialized = true;
        }

        glClearColor(0.03f, 0.05f, 0.08f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        const Mat4 viewProjection = camera.ViewProjectionMatrix();
        if (sceneInitialized) {
            sceneRenderer.Draw(viewProjection);
        }

        particleRenderer.Draw(viewProjection, 3.0f);

        ViewerMetrics metrics;
        metrics.renderMs = lastRenderMs;
        if (deltaSeconds > 0.0) {
            const double fpsNow = 1.0 / deltaSeconds;
            fpsSmoothed = (fpsSmoothed * 0.9) + (fpsNow * 0.1);
            metrics.fps = fpsSmoothed;
        }

        if (options.replayMode) {
            metrics.replayMode = true;
            metrics.paused = replayPaused;
            metrics.playbackRate = replayRate;
            metrics.producedSnapshots = replayFrames.size();
            metrics.consumedSnapshots = replayIndex;
        } else if (liveRuntime != nullptr) {
            metrics.replayMode = false;
            metrics.producedSnapshots = liveRuntime->ring.ProducedCount();
            metrics.consumedSnapshots = liveRuntime->ring.ConsumedCount();
            metrics.droppedSnapshots = liveRuntime->ring.DroppedCount();
            metrics.pendingSnapshots = liveRuntime->ring.PendingCount();

            std::lock_guard<std::mutex> lock(liveRuntime->errorMutex);
            if (!liveRuntime->writerError.empty()) {
                ImGui::Begin("Recorder Status");
                ImGui::TextWrapped("Recording disabled: %s", liveRuntime->writerError.c_str());
                ImGui::End();
            }
        }

        if (haveSnapshot) {
            DrawTelemetryPanels(currentSnapshot, metrics, telemetryUi);
        }

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);

        const auto frameEnd = std::chrono::steady_clock::now();
        lastRenderMs = std::chrono::duration<double, std::milli>(frameEnd - frameStart).count();
    }

    if (liveRuntime != nullptr) {
        liveRuntime->running.store(false, std::memory_order_relaxed);
    }

    if (simulationThread.joinable()) {
        simulationThread.join();
    }

    sceneRenderer.Shutdown();
    particleRenderer.Shutdown();

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
