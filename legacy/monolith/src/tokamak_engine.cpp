#include "tokamak_engine.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>

void TokamakEngine::ParticleSystem::AddParticle(
    Vec3 pos,
    Vec3 vel,
    float mass,
    float charge,
    ParticleType type) {
    if (positions.size() >= MAX_PARTICLES) {
        return;
    }

    positions.push_back(pos);
    velocities.push_back(vel);
    masses.push_back(mass);
    charges.push_back(charge);
    qOverM.push_back(charge / mass);
    species.push_back(type);
}

void TokamakEngine::ParticleSystem::MarkDead(size_t index) {
    species[index] = ParticleType::Dead;
}

void TokamakEngine::ParticleSystem::Compact() {
    for (size_t i = 0; i < positions.size();) {
        if (species[i] == ParticleType::Dead) {
            positions[i] = positions.back();
            positions.pop_back();

            velocities[i] = velocities.back();
            velocities.pop_back();

            masses[i] = masses.back();
            masses.pop_back();

            charges[i] = charges.back();
            charges.pop_back();

            qOverM[i] = qOverM.back();
            qOverM.pop_back();

            species[i] = species.back();
            species.pop_back();
        } else {
            ++i;
        }
    }
}

TokamakEngine::SpatialGrid::SpatialGrid(float reactorSize, float cSize) {
    cellSize = std::max(cSize, reactorSize / 50.0f);

    gridWidth = static_cast<int>(std::ceil((reactorSize * 2.0f) / cellSize));
    totalCells = gridWidth * gridWidth * gridWidth;

    cellCounts.resize(totalCells, 0);
    cellOffsets.resize(totalCells, 0);
}

int TokamakEngine::SpatialGrid::GetCellIndex(const Vec3& pos, float offset) const {
    int x = static_cast<int>(std::floor((pos.x + offset) / cellSize));
    int y = static_cast<int>(std::floor((pos.y + offset) / cellSize));
    int z = static_cast<int>(std::floor((pos.z + offset) / cellSize));

    x = std::max(0, std::min(x, gridWidth - 1));
    y = std::max(0, std::min(y, gridWidth - 1));
    z = std::max(0, std::min(z, gridWidth - 1));

    return x + gridWidth * (y + gridWidth * z);
}

TokamakEngine::TokamakEngine(Scenario scenario, uint64_t seed) {
    if (seed == 0) {
        rng_.seed(std::random_device{}());
    } else {
        rng_.seed(static_cast<uint32_t>(seed));
    }

    grid_ = new SpatialGrid(config_.majorRadius + config_.minorRadius, 0.2f);

    std::cout << "--- INITIALIZING TOKAMAK SIMULATION ---\n";

    switch (scenario) {
        case Scenario::COLD_VACUUM:
            std::cout << "SCENARIO: Cold Vacuum. No NBI Heating.\n";
            nbi_.isActive = false;
            break;
        case Scenario::NBI_IGNITION:
            std::cout << "SCENARIO: NBI Ignition. High power heating to target fusion.\n";
            nbi_.isActive = true;
            nbi_.beamEnergyKeV = 120.0f;
            break;
        case Scenario::MAGNETIC_FAILURE:
            std::cout << "SCENARIO: Magnetic Failure. Poloidal field collapsing.\n";
            nbi_.isActive = true;
            config_.plasmaCurrent = 0.0f;
            break;
    }

    std::uniform_real_distribution<float> angleDist(0.0f, 2.0f * PI);
    std::uniform_real_distribution<float> rDist(0.0f, config_.minorRadius * 0.8f);

    for (int i = 0; i < 1000; ++i) {
        const float a1 = angleDist(rng_);
        const float a2 = angleDist(rng_);
        const float rad = rDist(rng_);

        Vec3 pos(
            (config_.majorRadius + rad * std::cos(a2)) * std::cos(a1),
            (config_.majorRadius + rad * std::cos(a2)) * std::sin(a1),
            rad * std::sin(a2));
        Vec3 vel(0.0f, 0.0f, 0.0f);

        particles_.AddParticle(
            pos,
            vel,
            (i % 2 == 0) ? MASS_D : MASS_T,
            E_CHARGE,
            (i % 2 == 0) ? ParticleType::Deuterium : ParticleType::Tritium);
    }
}

TokamakEngine::~TokamakEngine() {
    delete grid_;
    grid_ = nullptr;
}

Vec3 TokamakEngine::CalculateBField(const Vec3& pos) const {
    const float rMajor = std::sqrt(pos.x * pos.x + pos.y * pos.y);
    if (rMajor < 0.001f) {
        return Vec3(0.0f, 0.0f, 0.0f);
    }

    const float bToroidalMag = (MU_0 * config_.toroidalCoilTurns * config_.toroidalCurrent) /
        (2.0f * PI * rMajor);
    const Vec3 bToroidal(-bToroidalMag * (pos.y / rMajor), bToroidalMag * (pos.x / rMajor), 0.0f);

    const float rMinor = std::sqrt(
        std::pow(rMajor - config_.majorRadius, 2.0f) + (pos.z * pos.z));
    if (rMinor < 0.001f) {
        return bToroidal;
    }

    const float bPoloidalMag = (MU_0 * config_.plasmaCurrent) / (2.0f * PI * rMinor);
    const float bR = -bPoloidalMag * (pos.z / rMinor);
    const float bZ = bPoloidalMag * ((rMajor - config_.majorRadius) / rMinor);
    const Vec3 bPoloidal(bR * (pos.x / rMajor), bR * (pos.y / rMajor), bZ);

    return bToroidal + bPoloidal;
}

Vec3 TokamakEngine::CalculateEField(const Vec3& pos) const {
    const float rMajor = std::sqrt(pos.x * pos.x + pos.y * pos.y);
    if (rMajor < 0.001f) {
        return Vec3(0.0f, 0.0f, 0.0f);
    }

    const Vec3 coreCenter(
        config_.majorRadius * (pos.x / rMajor),
        config_.majorRadius * (pos.y / rMajor),
        0.0f);
    const Vec3 dirToCenter = coreCenter - pos;
    return dirToCenter * 50.0f;
}

void TokamakEngine::PushParticles(float dt) {
    const float halfDt = dt * 0.5f;

    for (size_t i = 0; i < particles_.positions.size(); ++i) {
        const Vec3 pos = particles_.positions[i];
        const Vec3 vel = particles_.velocities[i];
        const float qOverM = particles_.qOverM[i];

        const Vec3 eField = CalculateEField(pos);
        const Vec3 bField = CalculateBField(pos);

        const Vec3 vMinus = vel + eField * (qOverM * halfDt);
        const Vec3 t = bField * (qOverM * halfDt);
        const Vec3 vPrime = vMinus + Vec3::cross(vMinus, t);
        const Vec3 s = t * (2.0f / (1.0f + Vec3::dot(t, t)));
        Vec3 vNew = vMinus + Vec3::cross(vPrime, s);
        vNew = vNew + eField * (qOverM * halfDt);

        Vec3 posNew = pos + (vNew * dt);

        const float rPos = std::sqrt((posNew.x * posNew.x) + (posNew.y * posNew.y));
        const float rTube = std::sqrt(
            std::pow(rPos - config_.majorRadius, 2.0f) + (posNew.z * posNew.z));

        if (rTube >= config_.minorRadius) {
            const float safeRPos = std::max(rPos, 0.001f);
            const Vec3 coreCenter(
                config_.majorRadius * (posNew.x / safeRPos),
                config_.majorRadius * (posNew.y / safeRPos),
                0.0f);
            const Vec3 normal = (posNew - coreCenter).normalize();
            vNew = vNew - normal * (2.0f * Vec3::dot(vNew, normal));
            posNew = coreCenter + normal * (config_.minorRadius * 0.99f);
        }

        particles_.velocities[i] = vNew;
        particles_.positions[i] = posNew;
    }
}

void TokamakEngine::SortParticlesIntoGrid() {
    const size_t numParticles = particles_.positions.size();
    const float offset = config_.majorRadius + config_.minorRadius;

    if (grid_->sortedParticleIDs.size() < numParticles) {
        grid_->sortedParticleIDs.resize(numParticles);
    }

    std::fill(grid_->cellCounts.begin(), grid_->cellCounts.end(), 0);

    for (size_t i = 0; i < numParticles; ++i) {
        const int cellIdx = grid_->GetCellIndex(particles_.positions[i], offset);
        grid_->cellCounts[cellIdx]++;
    }

    grid_->cellOffsets[0] = 0;
    for (int i = 1; i < grid_->totalCells; ++i) {
        grid_->cellOffsets[i] = grid_->cellOffsets[i - 1] + grid_->cellCounts[i - 1];
    }

    std::vector<uint32_t> currentOffsets = grid_->cellOffsets;
    for (size_t i = 0; i < numParticles; ++i) {
        const int cellIdx = grid_->GetCellIndex(particles_.positions[i], offset);
        const uint32_t dst = currentOffsets[cellIdx]++;
        grid_->sortedParticleIDs[dst] = static_cast<uint32_t>(i);
    }
}

void TokamakEngine::RunCollisions(float dt) {
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    for (int cell = 0; cell < grid_->totalCells; ++cell) {
        const uint32_t count = grid_->cellCounts[cell];
        if (count < 2) {
            continue;
        }

        const uint32_t startIdx = grid_->cellOffsets[cell];

        for (uint32_t i = 0; i + 1 < count; i += 2) {
            const size_t p1 = grid_->sortedParticleIDs[startIdx + i];
            const size_t p2 = grid_->sortedParticleIDs[startIdx + i + 1];

            const bool isDT =
                ((particles_.species[p1] == ParticleType::Deuterium) &&
                 (particles_.species[p2] == ParticleType::Tritium)) ||
                ((particles_.species[p1] == ParticleType::Tritium) &&
                 (particles_.species[p2] == ParticleType::Deuterium));
            if (!isDT) {
                continue;
            }

            const Vec3 vRel = particles_.velocities[p1] - particles_.velocities[p2];
            const float kineticJoules =
                0.5f * ((MASS_D * MASS_T) / (MASS_D + MASS_T)) * Vec3::dot(vRel, vRel);
            const float kineticKeV = kineticJoules / (1000.0f * E_CHARGE);

            if (kineticKeV <= 15.0f) {
                continue;
            }

            const float probability = 0.05f * dt * particles_.macroWeight;
            if (dist(rng_) >= probability) {
                continue;
            }

            particles_.fusionCountTotal++;

            const Vec3 centerPos = (particles_.positions[p1] + particles_.positions[p2]) * 0.5f;
            const Vec3 centerVel =
                (particles_.velocities[p1] * MASS_D + particles_.velocities[p2] * MASS_T) /
                (MASS_D + MASS_T);

            const float alphaEnergyJ = 3.5e6f * E_CHARGE;
            const float alphaSpeed = std::sqrt((2.0f * alphaEnergyJ) / MASS_HE);

            Vec3 randomDir(dist(rng_) - 0.5f, dist(rng_) - 0.5f, dist(rng_) - 0.5f);
            if (randomDir.mag() < 1.0e-5f) {
                randomDir = Vec3(1.0f, 0.0f, 0.0f);
            }

            const Vec3 heVel = centerVel + (randomDir.normalize() * alphaSpeed);

            particles_.MarkDead(p1);
            particles_.MarkDead(p2);
            particles_.AddParticle(centerPos, heVel, MASS_HE, E_CHARGE * 2.0f, ParticleType::Helium);

            // Prevent chained modifications in the same cell pass.
            break;
        }
    }
}

void TokamakEngine::RunNBI(float /*dt*/) {
    if (!nbi_.isActive) {
        return;
    }

    const float energyJoules = nbi_.beamEnergyKeV * 1000.0f * E_CHARGE;
    const float speed = std::sqrt((2.0f * energyJoules) / MASS_D);
    const Vec3 beamVel = nbi_.injectionNormal * speed;

    std::uniform_real_distribution<float> dist(-0.1f, 0.1f);

    for (int i = 0; i < nbi_.particlesPerStep; ++i) {
        Vec3 spawnPos = nbi_.injectorPos +
            (nbi_.injectionNormal * config_.minorRadius) +
            Vec3(dist(rng_), dist(rng_), dist(rng_));

        particles_.AddParticle(spawnPos, beamVel, MASS_D, E_CHARGE, ParticleType::Deuterium);

        spawnPos = spawnPos + Vec3(dist(rng_), dist(rng_), dist(rng_));
        particles_.AddParticle(spawnPos, beamVel, MASS_T, E_CHARGE, ParticleType::Tritium);
    }
}

TokamakEngine::TelemetryAccumulator TokamakEngine::ComputeTelemetry() const {
    TelemetryAccumulator accum;

    for (size_t i = 0; i < particles_.positions.size(); ++i) {
        const ParticleType type = particles_.species[i];
        if (type == ParticleType::Dead) {
            continue;
        }

        if (type == ParticleType::Deuterium) {
            accum.dCount++;
        } else if (type == ParticleType::Tritium) {
            accum.tCount++;
        } else if (type == ParticleType::Helium) {
            accum.heCount++;
        }

        const float speedSq = Vec3::dot(particles_.velocities[i], particles_.velocities[i]);
        accum.totalKineticEnergy += 0.5 * static_cast<double>(particles_.masses[i]) * speedSq;
    }

    return accum;
}

void TokamakEngine::FillTelemetrySummary(
    TelemetrySummary& summary,
    const TelemetryAccumulator& accum) const {
    summary.deuteriumCount = accum.dCount;
    summary.tritiumCount = accum.tCount;
    summary.heliumCount = accum.heCount;
    summary.totalIons = accum.dCount + accum.tCount + accum.heCount;
    summary.fusionEvents = particles_.fusionCountTotal;

    if (summary.totalIons == 0) {
        summary.avgTempKeV = 0.0;
    } else {
        summary.avgTempKeV =
            (accum.totalKineticEnergy / static_cast<double>(summary.totalIons)) /
            (1000.0 * E_CHARGE);
    }
}

void TokamakEngine::Step(float dt) {
    const auto start = std::chrono::high_resolution_clock::now();

    RunNBI(dt);
    PushParticles(dt);
    SortParticlesIntoGrid();
    RunCollisions(dt);
    particles_.Compact();

    timeSeconds_ += dt;
    stepIndex_++;

    const auto end = std::chrono::high_resolution_clock::now();
    lastStepMs_ = std::chrono::duration<double, std::milli>(end - start).count();
}

bool TokamakEngine::ExportSnapshot(SimulationSnapshot& out, const SnapshotRequest& request) const {
    const auto exportStart = std::chrono::high_resolution_clock::now();

    const size_t total = particles_.positions.size();
    const uint32_t maxSample = std::max<uint32_t>(1, request.maxSampledParticles);
    const uint32_t sampleCount = static_cast<uint32_t>(std::min<size_t>(total, maxSample));
    const uint32_t stride = (sampleCount == 0)
        ? 1
        : static_cast<uint32_t>((total + sampleCount - 1) / sampleCount);

    out.stepIndex = stepIndex_;
    out.simTimeSeconds = timeSeconds_;
    out.config = config_;
    out.nbi = nbi_;
    out.totalParticleCount = static_cast<uint32_t>(total);
    out.sampleStride = std::max<uint32_t>(stride, 1);

    const TelemetryAccumulator telemetry = ComputeTelemetry();
    FillTelemetrySummary(out.telemetry, telemetry);

    out.sampledParticles.clear();
    if (out.sampledParticles.capacity() < sampleCount) {
        out.sampledParticles.reserve(sampleCount);
    }

    if (sampleCount > 0) {
        for (size_t i = 0; i < total && out.sampledParticles.size() < sampleCount; i += out.sampleStride) {
            ParticleSample sample;
            sample.x = particles_.positions[i].x;
            sample.y = particles_.positions[i].y;
            sample.z = particles_.positions[i].z;

            if (request.includeVelocities) {
                sample.vx = particles_.velocities[i].x;
                sample.vy = particles_.velocities[i].y;
                sample.vz = particles_.velocities[i].z;
            }

            if (request.includeMassAndCharge) {
                sample.mass = particles_.masses[i];
                sample.charge = particles_.charges[i];
                sample.qOverM = particles_.qOverM[i];
            }

            sample.species = static_cast<uint8_t>(particles_.species[i]);
            out.sampledParticles.push_back(sample);
        }
    }

    const auto exportEnd = std::chrono::high_resolution_clock::now();
    out.performance.simStepMs = lastStepMs_;
    out.performance.exportMs = std::chrono::duration<double, std::milli>(exportEnd - exportStart).count();

    return true;
}

void TokamakEngine::PrintTelemetry(uint64_t step) const {
    TelemetrySummary summary;
    FillTelemetrySummary(summary, ComputeTelemetry());

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "[Step " << std::setw(6) << step
              << " | Time " << (timeSeconds_ * 1000.0) << " ms] "
              << "Total Ions: " << std::setw(7) << summary.totalIons << " | "
              << "D: " << summary.deuteriumCount << " "
              << "T: " << summary.tritiumCount << " "
              << "He(Ash): " << summary.heliumCount << " | "
              << "Avg Temp: " << std::setw(8) << summary.avgTempKeV << " keV | "
              << "Fusion Events: " << summary.fusionEvents
              << std::endl;
}
