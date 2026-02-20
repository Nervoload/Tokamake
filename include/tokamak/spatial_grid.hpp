#pragma once

#include <cstdint>
#include <vector>

#include "tokamak/particle_system.hpp"

namespace tokamak {

class SpatialGrid {
public:
    SpatialGrid(float reactorSize_m, float requestedCellSize_m);

    int GetCellIndex(const Vec3& position, float offset_m, bool* wasClamped = nullptr) const;
    bool IsValidCellIndex(int cellIndex) const;

    void EnsureParticleCapacity(std::size_t particleCount);
    void ResetCounts();
    void BuildOffsets();
    void ResetWriteHeads();

    float cellSize_m;
    int gridWidth;
    int totalCells;

    std::vector<uint32_t> cellCounts;
    std::vector<uint32_t> cellOffsets;
    std::vector<uint32_t> sortedParticleIDs;
    std::vector<uint32_t> writeHeads;
};

uint64_t SortParticlesIntoGrid(const ParticleSystem& particles, SpatialGrid& grid, float gridOffset_m);

}  // namespace tokamak
