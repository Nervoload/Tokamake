#pragma once

#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "sim_snapshot.h"

class TokamakEngine {
public:
    explicit TokamakEngine(Scenario scenario, uint64_t seed = 0);
    ~TokamakEngine();

    void Step(float dt);
    bool ExportSnapshot(SimulationSnapshot& out, const SnapshotRequest& request) const;
    void PrintTelemetry(uint64_t step) const;

    uint64_t StepIndex() const { return stepIndex_; }
    double TimeSeconds() const { return timeSeconds_; }

    const TokamakConfig& Config() const { return config_; }
    const NBIConfig& NeutralBeamConfig() const { return nbi_; }

private:
    static constexpr float PI = 3.14159265359f;
    static constexpr float MU_0 = 1.25663706e-6f;
    static constexpr float E_CHARGE = 1.60217663e-19f;
    static constexpr float MASS_D = 3.3435e-27f;
    static constexpr float MASS_T = 5.0082e-27f;
    static constexpr float MASS_HE = 6.6464e-27f;
    static constexpr size_t MAX_PARTICLES = 2500000;

    struct ParticleSystem {
        std::vector<Vec3> positions;
        std::vector<Vec3> velocities;
        std::vector<float> masses;
        std::vector<float> charges;
        std::vector<float> qOverM;
        std::vector<ParticleType> species;

        float macroWeight = 1.0e12f;
        uint64_t fusionCountTotal = 0;

        void AddParticle(Vec3 pos, Vec3 vel, float mass, float charge, ParticleType type);
        void MarkDead(size_t index);
        void Compact();
    };

    struct SpatialGrid {
        float cellSize = 0.2f;
        int gridWidth = 0;
        int totalCells = 0;
        std::vector<uint32_t> cellCounts;
        std::vector<uint32_t> cellOffsets;
        std::vector<uint32_t> sortedParticleIDs;

        SpatialGrid(float reactorSize, float cSize);
        int GetCellIndex(const Vec3& pos, float offset) const;
    };

    struct TelemetryAccumulator {
        uint32_t dCount = 0;
        uint32_t tCount = 0;
        uint32_t heCount = 0;
        double totalKineticEnergy = 0.0;
    };

    Vec3 CalculateBField(const Vec3& pos) const;
    Vec3 CalculateEField(const Vec3& pos) const;

    void PushParticles(float dt);
    void SortParticlesIntoGrid();
    void RunCollisions(float dt);
    void RunNBI(float dt);

    TelemetryAccumulator ComputeTelemetry() const;
    void FillTelemetrySummary(TelemetrySummary& summary, const TelemetryAccumulator& accum) const;

    TokamakConfig config_;
    NBIConfig nbi_;
    ParticleSystem particles_;
    SpatialGrid* grid_ = nullptr;

    mutable double lastStepMs_ = 0.0;
    uint64_t stepIndex_ = 0;
    double timeSeconds_ = 0.0;
    std::mt19937 rng_;
};
