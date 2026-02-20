#include "tokamak/spatial_grid.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

namespace tokamak {

SpatialGrid::SpatialGrid(float reactorSize_m, float requestedCellSize_m)
    : cellSize_m(std::max(requestedCellSize_m, reactorSize_m / 50.0f)),
      gridWidth(static_cast<int>(std::ceil((reactorSize_m * 2.0f) / cellSize_m))),
      totalCells(gridWidth * gridWidth * gridWidth),
      cellCounts(static_cast<std::size_t>(totalCells), 0),
      cellOffsets(static_cast<std::size_t>(totalCells), 0),
      writeHeads(static_cast<std::size_t>(totalCells), 0) {}

int SpatialGrid::GetCellIndex(const Vec3& position, float offset_m, bool* wasClamped) const {
    int x = static_cast<int>(std::floor((position.x + offset_m) / cellSize_m));
    int y = static_cast<int>(std::floor((position.y + offset_m) / cellSize_m));
    int z = static_cast<int>(std::floor((position.z + offset_m) / cellSize_m));

    const bool clamped = (x < 0 || x >= gridWidth || y < 0 || y >= gridWidth || z < 0 || z >= gridWidth);
    if (wasClamped != nullptr) {
        *wasClamped = clamped;
    }

    x = std::max(0, std::min(x, gridWidth - 1));
    y = std::max(0, std::min(y, gridWidth - 1));
    z = std::max(0, std::min(z, gridWidth - 1));

    return x + (gridWidth * (y + (gridWidth * z)));
}

bool SpatialGrid::IsValidCellIndex(int cellIndex) const {
    return cellIndex >= 0 && cellIndex < totalCells;
}

void SpatialGrid::EnsureParticleCapacity(std::size_t particleCount) {
    if (sortedParticleIDs.size() < particleCount) {
        sortedParticleIDs.resize(particleCount);
    }
}

void SpatialGrid::ResetCounts() {
    std::fill(cellCounts.begin(), cellCounts.end(), 0);
}

void SpatialGrid::BuildOffsets() {
    if (cellOffsets.empty()) {
        return;
    }

    cellOffsets[0] = 0;
    for (int i = 1; i < totalCells; ++i) {
        cellOffsets[static_cast<std::size_t>(i)] =
            cellOffsets[static_cast<std::size_t>(i - 1)] + cellCounts[static_cast<std::size_t>(i - 1)];
    }
}

void SpatialGrid::ResetWriteHeads() {
    std::copy(cellOffsets.begin(), cellOffsets.end(), writeHeads.begin());
}

uint64_t SortParticlesIntoGrid(const ParticleSystem& particles, SpatialGrid& grid, float gridOffset_m) {
    const auto& positions = particles.Positions();
    const std::size_t particleCount = positions.size();
    uint64_t outOfDomainClampCount = 0;

    grid.EnsureParticleCapacity(particleCount);
    grid.ResetCounts();

    for (std::size_t i = 0; i < particleCount; ++i) {
        bool wasClamped = false;
        const int cellIndex = grid.GetCellIndex(positions[i], gridOffset_m, &wasClamped);
        assert(grid.IsValidCellIndex(cellIndex));
        if (wasClamped) {
            ++outOfDomainClampCount;
        }
        ++grid.cellCounts[static_cast<std::size_t>(cellIndex)];
    }

    grid.BuildOffsets();
    grid.ResetWriteHeads();

    for (std::size_t i = 0; i < particleCount; ++i) {
        const int cellIndex = grid.GetCellIndex(positions[i], gridOffset_m);
        assert(grid.IsValidCellIndex(cellIndex));

        const auto destIndex = grid.writeHeads[static_cast<std::size_t>(cellIndex)]++;
        assert(destIndex < grid.sortedParticleIDs.size());
        grid.sortedParticleIDs[destIndex] = static_cast<uint32_t>(i);
    }

    return outOfDomainClampCount;
}

}  // namespace tokamak
