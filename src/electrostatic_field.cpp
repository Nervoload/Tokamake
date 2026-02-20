#include "tokamak/electrostatic_field.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace tokamak {
namespace {

constexpr double kEpsilon0_C2PerN_m2 = 8.8541878128e-12;

int ClampIndex(int value, int minValue, int maxValue) {
    return std::max(minValue, std::min(value, maxValue));
}

std::size_t CellIndex(int x, int y, int z, int width) {
    return static_cast<std::size_t>(x + (width * (y + (width * z))));
}

bool IsBoundaryCell(int x, int y, int z, int width) {
    return (x == 0) || (y == 0) || (z == 0) || (x == (width - 1)) || (y == (width - 1)) || (z == (width - 1));
}

void ClearAndResizeToCellCount(
    const ElectrostaticMeshGeometry& geometry,
    std::vector<double>* values) {
    if (values == nullptr) {
        return;
    }
    values->assign(ElectrostaticCellCount(geometry), 0.0);
}

void EnforceDirichletBoundary(
    int width,
    double boundaryPotential_V,
    std::vector<double>* potential_V) {
    if (potential_V == nullptr || potential_V->empty()) {
        return;
    }
    for (int z = 0; z < width; ++z) {
        for (int y = 0; y < width; ++y) {
            for (int x = 0; x < width; ++x) {
                if (!IsBoundaryCell(x, y, z, width)) {
                    continue;
                }
                (*potential_V)[CellIndex(x, y, z, width)] = boundaryPotential_V;
            }
        }
    }
}

void EnforceNeumannBoundary(int width, std::vector<double>* potential_V) {
    if (potential_V == nullptr || potential_V->empty() || width < 2) {
        return;
    }

    for (int z = 0; z < width; ++z) {
        for (int y = 0; y < width; ++y) {
            for (int x = 0; x < width; ++x) {
                if (!IsBoundaryCell(x, y, z, width)) {
                    continue;
                }
                const int sourceX = ClampIndex(x, 1, width - 2);
                const int sourceY = ClampIndex(y, 1, width - 2);
                const int sourceZ = ClampIndex(z, 1, width - 2);
                (*potential_V)[CellIndex(x, y, z, width)] = (*potential_V)[CellIndex(sourceX, sourceY, sourceZ, width)];
            }
        }
    }
}

double NeighborPotential(
    int x,
    int y,
    int z,
    int width,
    const std::vector<double>& potential_V,
    ElectrostaticBoundaryCondition boundaryCondition,
    double dirichletBoundaryPotential_V,
    double fallbackCurrentValue_V) {
    if (x >= 0 && x < width && y >= 0 && y < width && z >= 0 && z < width) {
        return potential_V[CellIndex(x, y, z, width)];
    }

    if (boundaryCondition == ElectrostaticBoundaryCondition::DirichletZero) {
        return dirichletBoundaryPotential_V;
    }

    return fallbackCurrentValue_V;
}

double ComputeResidualL2(
    const ElectrostaticMeshGeometry& geometry,
    ElectrostaticBoundaryCondition boundaryCondition,
    double dirichletBoundaryPotential_V,
    const std::vector<double>& chargeDensity_CPerM3,
    const std::vector<double>& potential_V) {
    const int width = geometry.gridWidth;
    if (width <= 0 || geometry.cellSize_m <= 0.0 ||
        chargeDensity_CPerM3.size() != potential_V.size()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double inverseCellSq = 1.0 / (geometry.cellSize_m * geometry.cellSize_m);
    double residualSqSum = 0.0;
    uint64_t sampleCount = 0;

    for (int z = 0; z < width; ++z) {
        for (int y = 0; y < width; ++y) {
            for (int x = 0; x < width; ++x) {
                if (boundaryCondition == ElectrostaticBoundaryCondition::DirichletZero &&
                    IsBoundaryCell(x, y, z, width)) {
                    continue;
                }

                const std::size_t idx = CellIndex(x, y, z, width);
                const double center = potential_V[idx];
                const double sumNeighbors =
                    NeighborPotential(x + 1, y, z, width, potential_V, boundaryCondition, dirichletBoundaryPotential_V, center) +
                    NeighborPotential(x - 1, y, z, width, potential_V, boundaryCondition, dirichletBoundaryPotential_V, center) +
                    NeighborPotential(x, y + 1, z, width, potential_V, boundaryCondition, dirichletBoundaryPotential_V, center) +
                    NeighborPotential(x, y - 1, z, width, potential_V, boundaryCondition, dirichletBoundaryPotential_V, center) +
                    NeighborPotential(x, y, z + 1, width, potential_V, boundaryCondition, dirichletBoundaryPotential_V, center) +
                    NeighborPotential(x, y, z - 1, width, potential_V, boundaryCondition, dirichletBoundaryPotential_V, center);

                const double laplacian = (sumNeighbors - (6.0 * center)) * inverseCellSq;
                const double residual = laplacian + (chargeDensity_CPerM3[idx] / kEpsilon0_C2PerN_m2);
                residualSqSum += residual * residual;
                ++sampleCount;
            }
        }
    }

    if (sampleCount == 0) {
        return 0.0;
    }
    return std::sqrt(residualSqSum / static_cast<double>(sampleCount));
}

}  // namespace

std::size_t ElectrostaticCellCount(const ElectrostaticMeshGeometry& geometry) {
    if (geometry.gridWidth <= 0) {
        return 0;
    }
    const int width = geometry.gridWidth;
    return static_cast<std::size_t>(width) * static_cast<std::size_t>(width) * static_cast<std::size_t>(width);
}

void DepositChargeDensity(
    const std::vector<Vec3>& positions,
    const std::vector<float>& charges_C,
    float macroWeight,
    const ElectrostaticMeshGeometry& geometry,
    ChargeAssignmentScheme assignmentScheme,
    std::vector<double>* chargeDensity_CPerM3) {
    ClearAndResizeToCellCount(geometry, chargeDensity_CPerM3);
    if (chargeDensity_CPerM3 == nullptr || geometry.gridWidth <= 0 || geometry.cellSize_m <= 0.0) {
        return;
    }

    const int width = geometry.gridWidth;
    const double inverseCellVolume = 1.0 / (geometry.cellSize_m * geometry.cellSize_m * geometry.cellSize_m);
    const std::size_t count = std::min(positions.size(), charges_C.size());

    for (std::size_t i = 0; i < count; ++i) {
        const double charge_C = static_cast<double>(charges_C[i]) * static_cast<double>(macroWeight);
        if (!std::isfinite(charge_C) || charge_C == 0.0) {
            continue;
        }

        const double gridX = (static_cast<double>(positions[i].x) + geometry.offset_m) / geometry.cellSize_m;
        const double gridY = (static_cast<double>(positions[i].y) + geometry.offset_m) / geometry.cellSize_m;
        const double gridZ = (static_cast<double>(positions[i].z) + geometry.offset_m) / geometry.cellSize_m;

        if (assignmentScheme == ChargeAssignmentScheme::NGP) {
            const int x = ClampIndex(static_cast<int>(std::llround(gridX)), 0, width - 1);
            const int y = ClampIndex(static_cast<int>(std::llround(gridY)), 0, width - 1);
            const int z = ClampIndex(static_cast<int>(std::llround(gridZ)), 0, width - 1);
            (*chargeDensity_CPerM3)[CellIndex(x, y, z, width)] += charge_C * inverseCellVolume;
            continue;
        }

        const int x0 = static_cast<int>(std::floor(gridX));
        const int y0 = static_cast<int>(std::floor(gridY));
        const int z0 = static_cast<int>(std::floor(gridZ));
        const double tx = gridX - static_cast<double>(x0);
        const double ty = gridY - static_cast<double>(y0);
        const double tz = gridZ - static_cast<double>(z0);

        struct Contribution {
            int x;
            int y;
            int z;
            double weight;
        };

        Contribution contributions[8];
        int validCount = 0;
        double validWeightSum = 0.0;

        for (int dz = 0; dz <= 1; ++dz) {
            for (int dy = 0; dy <= 1; ++dy) {
                for (int dx = 0; dx <= 1; ++dx) {
                    const int ix = x0 + dx;
                    const int iy = y0 + dy;
                    const int iz = z0 + dz;
                    const double wx = (dx == 0) ? (1.0 - tx) : tx;
                    const double wy = (dy == 0) ? (1.0 - ty) : ty;
                    const double wz = (dz == 0) ? (1.0 - tz) : tz;
                    const double weight = wx * wy * wz;
                    if (ix < 0 || ix >= width || iy < 0 || iy >= width || iz < 0 || iz >= width || weight <= 0.0) {
                        continue;
                    }
                    contributions[validCount++] = Contribution{ix, iy, iz, weight};
                    validWeightSum += weight;
                }
            }
        }

        if (validCount == 0 || validWeightSum <= 0.0) {
            const int x = ClampIndex(x0, 0, width - 1);
            const int y = ClampIndex(y0, 0, width - 1);
            const int z = ClampIndex(z0, 0, width - 1);
            (*chargeDensity_CPerM3)[CellIndex(x, y, z, width)] += charge_C * inverseCellVolume;
            continue;
        }

        for (int c = 0; c < validCount; ++c) {
            const Contribution& contribution = contributions[c];
            const double normalizedWeight = contribution.weight / validWeightSum;
            (*chargeDensity_CPerM3)[CellIndex(contribution.x, contribution.y, contribution.z, width)] +=
                charge_C * normalizedWeight * inverseCellVolume;
        }
    }
}

void ApplyNeutralizingBackground(
    double neutralizingBackgroundFraction,
    std::vector<double>* chargeDensity_CPerM3) {
    if (chargeDensity_CPerM3 == nullptr || chargeDensity_CPerM3->empty()) {
        return;
    }
    if (!std::isfinite(neutralizingBackgroundFraction) || neutralizingBackgroundFraction <= 0.0) {
        return;
    }

    const double clampedFraction = std::max(0.0, std::min(neutralizingBackgroundFraction, 1.0));
    double sum = 0.0;
    for (const double value : *chargeDensity_CPerM3) {
        sum += value;
    }
    const double mean = sum / static_cast<double>(chargeDensity_CPerM3->size());
    const double subtract = clampedFraction * mean;
    for (double& value : *chargeDensity_CPerM3) {
        value -= subtract;
    }
}

ElectrostaticSolveResult SolvePoissonSor(
    const ElectrostaticMeshGeometry& geometry,
    const ElectrostaticSolveConfig& solveConfig,
    double dirichletBoundaryPotential_V,
    const std::vector<double>& chargeDensity_CPerM3,
    std::vector<double>* potential_V) {
    ElectrostaticSolveResult result;
    result.residual.solverKind = ResidualSolverKind::Sor;
    result.residual.tolerance = solveConfig.tolerance;
    result.residual.status = ResidualStatus::Unavailable;

    const std::size_t expectedCellCount = ElectrostaticCellCount(geometry);
    if (potential_V == nullptr || expectedCellCount == 0 || geometry.cellSize_m <= 0.0 ||
        chargeDensity_CPerM3.size() != expectedCellCount) {
        result.residual.status = ResidualStatus::Failed;
        return result;
    }

    potential_V->resize(expectedCellCount, 0.0);
    if (solveConfig.boundaryCondition == ElectrostaticBoundaryCondition::DirichletZero) {
        EnforceDirichletBoundary(geometry.gridWidth, dirichletBoundaryPotential_V, potential_V);
    }

    const int width = geometry.gridWidth;
    const double cellSq = geometry.cellSize_m * geometry.cellSize_m;
    const double omega = std::max(0.1, std::min(1.99, solveConfig.sorOmega));
    const uint32_t maxIterations = std::max<uint32_t>(1, solveConfig.maxIterations);
    const double tolerance = std::max(1.0e-14, solveConfig.tolerance);

    double latestResidual = std::numeric_limits<double>::quiet_NaN();

    for (uint32_t iter = 0; iter < maxIterations; ++iter) {
        for (int z = 0; z < width; ++z) {
            for (int y = 0; y < width; ++y) {
                for (int x = 0; x < width; ++x) {
                    if (solveConfig.boundaryCondition == ElectrostaticBoundaryCondition::DirichletZero &&
                        IsBoundaryCell(x, y, z, width)) {
                        continue;
                    }

                    const std::size_t idx = CellIndex(x, y, z, width);
                    const double center = (*potential_V)[idx];
                    const double neighborSum =
                        NeighborPotential(x + 1, y, z, width, *potential_V, solveConfig.boundaryCondition, dirichletBoundaryPotential_V, center) +
                        NeighborPotential(x - 1, y, z, width, *potential_V, solveConfig.boundaryCondition, dirichletBoundaryPotential_V, center) +
                        NeighborPotential(x, y + 1, z, width, *potential_V, solveConfig.boundaryCondition, dirichletBoundaryPotential_V, center) +
                        NeighborPotential(x, y - 1, z, width, *potential_V, solveConfig.boundaryCondition, dirichletBoundaryPotential_V, center) +
                        NeighborPotential(x, y, z + 1, width, *potential_V, solveConfig.boundaryCondition, dirichletBoundaryPotential_V, center) +
                        NeighborPotential(x, y, z - 1, width, *potential_V, solveConfig.boundaryCondition, dirichletBoundaryPotential_V, center);

                    const double poissonJacobi = (neighborSum + (cellSq * chargeDensity_CPerM3[idx] / kEpsilon0_C2PerN_m2)) / 6.0;
                    (*potential_V)[idx] = ((1.0 - omega) * center) + (omega * poissonJacobi);
                }
            }
        }

        if (solveConfig.boundaryCondition == ElectrostaticBoundaryCondition::DirichletZero) {
            EnforceDirichletBoundary(width, dirichletBoundaryPotential_V, potential_V);
        } else {
            EnforceNeumannBoundary(width, potential_V);
        }

        latestResidual = ComputeResidualL2(
            geometry,
            solveConfig.boundaryCondition,
            dirichletBoundaryPotential_V,
            chargeDensity_CPerM3,
            *potential_V);
        result.residualHistoryL2.push_back(latestResidual);

        if (std::isfinite(latestResidual) && latestResidual <= tolerance) {
            result.residual.status = ResidualStatus::Measured;
            result.residual.residualAvailable = true;
            result.residual.residualL2 = latestResidual;
            result.residual.iterations = iter + 1;
            result.residual.converged = true;
            return result;
        }
    }

    result.residual.iterations = maxIterations;
    result.residual.converged = false;
    if (std::isfinite(latestResidual)) {
        result.residual.status = ResidualStatus::Measured;
        result.residual.residualAvailable = true;
        result.residual.residualL2 = latestResidual;
    } else {
        result.residual.status = ResidualStatus::Failed;
        result.residual.residualAvailable = false;
    }
    return result;
}

void ReconstructElectricField(
    const ElectrostaticMeshGeometry& geometry,
    const std::vector<double>& potential_V,
    std::vector<Vec3>* electricField_VPerM) {
    if (electricField_VPerM == nullptr) {
        return;
    }
    electricField_VPerM->assign(ElectrostaticCellCount(geometry), Vec3(0.0f, 0.0f, 0.0f));
    if (geometry.gridWidth <= 0 || geometry.cellSize_m <= 0.0 || potential_V.size() != electricField_VPerM->size()) {
        return;
    }

    const int width = geometry.gridWidth;
    const double inv2h = 1.0 / (2.0 * geometry.cellSize_m);
    const double invh = 1.0 / geometry.cellSize_m;

    for (int z = 0; z < width; ++z) {
        for (int y = 0; y < width; ++y) {
            for (int x = 0; x < width; ++x) {
                const auto idx = CellIndex(x, y, z, width);

                double dPhiDx = 0.0;
                double dPhiDy = 0.0;
                double dPhiDz = 0.0;

                if (x == 0) {
                    dPhiDx = (potential_V[CellIndex(x + 1, y, z, width)] - potential_V[idx]) * invh;
                } else if (x == (width - 1)) {
                    dPhiDx = (potential_V[idx] - potential_V[CellIndex(x - 1, y, z, width)]) * invh;
                } else {
                    dPhiDx = (potential_V[CellIndex(x + 1, y, z, width)] - potential_V[CellIndex(x - 1, y, z, width)]) * inv2h;
                }

                if (y == 0) {
                    dPhiDy = (potential_V[CellIndex(x, y + 1, z, width)] - potential_V[idx]) * invh;
                } else if (y == (width - 1)) {
                    dPhiDy = (potential_V[idx] - potential_V[CellIndex(x, y - 1, z, width)]) * invh;
                } else {
                    dPhiDy = (potential_V[CellIndex(x, y + 1, z, width)] - potential_V[CellIndex(x, y - 1, z, width)]) * inv2h;
                }

                if (z == 0) {
                    dPhiDz = (potential_V[CellIndex(x, y, z + 1, width)] - potential_V[idx]) * invh;
                } else if (z == (width - 1)) {
                    dPhiDz = (potential_V[idx] - potential_V[CellIndex(x, y, z - 1, width)]) * invh;
                } else {
                    dPhiDz = (potential_V[CellIndex(x, y, z + 1, width)] - potential_V[CellIndex(x, y, z - 1, width)]) * inv2h;
                }

                (*electricField_VPerM)[idx] = Vec3(
                    static_cast<float>(-dPhiDx),
                    static_cast<float>(-dPhiDy),
                    static_cast<float>(-dPhiDz));
            }
        }
    }
}

Vec3 SampleElectricField(
    const ElectrostaticMeshGeometry& geometry,
    const std::vector<Vec3>& electricField_VPerM,
    const Vec3& position,
    ChargeAssignmentScheme assignmentScheme) {
    if (geometry.gridWidth <= 0 || geometry.cellSize_m <= 0.0 || electricField_VPerM.empty()) {
        return Vec3(0.0f, 0.0f, 0.0f);
    }

    const int width = geometry.gridWidth;
    const double gridX = (static_cast<double>(position.x) + geometry.offset_m) / geometry.cellSize_m;
    const double gridY = (static_cast<double>(position.y) + geometry.offset_m) / geometry.cellSize_m;
    const double gridZ = (static_cast<double>(position.z) + geometry.offset_m) / geometry.cellSize_m;

    if (assignmentScheme == ChargeAssignmentScheme::NGP) {
        const int x = ClampIndex(static_cast<int>(std::llround(gridX)), 0, width - 1);
        const int y = ClampIndex(static_cast<int>(std::llround(gridY)), 0, width - 1);
        const int z = ClampIndex(static_cast<int>(std::llround(gridZ)), 0, width - 1);
        return electricField_VPerM[CellIndex(x, y, z, width)];
    }

    const int x0 = static_cast<int>(std::floor(gridX));
    const int y0 = static_cast<int>(std::floor(gridY));
    const int z0 = static_cast<int>(std::floor(gridZ));
    const double tx = gridX - static_cast<double>(x0);
    const double ty = gridY - static_cast<double>(y0);
    const double tz = gridZ - static_cast<double>(z0);

    Vec3 weightedField(0.0f, 0.0f, 0.0f);
    double weightSum = 0.0;

    for (int dz = 0; dz <= 1; ++dz) {
        for (int dy = 0; dy <= 1; ++dy) {
            for (int dx = 0; dx <= 1; ++dx) {
                const int ix = x0 + dx;
                const int iy = y0 + dy;
                const int iz = z0 + dz;
                if (ix < 0 || ix >= width || iy < 0 || iy >= width || iz < 0 || iz >= width) {
                    continue;
                }

                const double wx = (dx == 0) ? (1.0 - tx) : tx;
                const double wy = (dy == 0) ? (1.0 - ty) : ty;
                const double wz = (dz == 0) ? (1.0 - tz) : tz;
                const double weight = wx * wy * wz;
                if (weight <= 0.0) {
                    continue;
                }

                const Vec3 cellField = electricField_VPerM[CellIndex(ix, iy, iz, width)];
                weightedField += cellField * static_cast<float>(weight);
                weightSum += weight;
            }
        }
    }

    if (weightSum <= 0.0) {
        const int x = ClampIndex(x0, 0, width - 1);
        const int y = ClampIndex(y0, 0, width - 1);
        const int z = ClampIndex(z0, 0, width - 1);
        return electricField_VPerM[CellIndex(x, y, z, width)];
    }

    return weightedField / static_cast<float>(weightSum);
}

}  // namespace tokamak
