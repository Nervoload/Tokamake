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
    bool clamped = false;
    const double maxIndex = static_cast<double>(gridWidth - 1);
    auto clampCellIndex = [&](double coordinate) -> int {
        if (!std::isfinite(coordinate) || coordinate < 0.0) {
            clamped = true;
            return 0;
        }
        if (coordinate > maxIndex) {
            clamped = true;
            return gridWidth - 1;
        }
        return static_cast<int>(coordinate);
    };

    const double xCoord = std::floor((static_cast<double>(position.x) + static_cast<double>(offset_m)) / static_cast<double>(cellSize_m));
    const double yCoord = std::floor((static_cast<double>(position.y) + static_cast<double>(offset_m)) / static_cast<double>(cellSize_m));
    const double zCoord = std::floor((static_cast<double>(position.z) + static_cast<double>(offset_m)) / static_cast<double>(cellSize_m));

    const int x = clampCellIndex(xCoord);
    const int y = clampCellIndex(yCoord);
    const int z = clampCellIndex(zCoord);

    if (wasClamped != nullptr) {
        *wasClamped = clamped;
    }

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
