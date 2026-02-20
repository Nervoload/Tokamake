#include "snapshot_stream.h"

#include <chrono>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <limits>

namespace {

constexpr uint32_t kSnapshotVersion = 1;
constexpr uint32_t kEndianMarker = 0x01020304u;
constexpr char kMagic[8] = {'T', 'K', 'S', 'N', 'A', 'P', '1', '\0'};

struct TokamakConfigWire {
    float majorRadius;
    float minorRadius;
    float toroidalCurrent;
    int32_t toroidalCoilTurns;
    float plasmaCurrent;
};

struct NBIConfigWire {
    uint32_t isActive;
    float beamEnergyKeV;
    int32_t particlesPerStep;
    float injectorPos[3];
    float injectionNormal[3];
};

struct SnapshotFileHeader {
    char magic[8];
    uint32_t version;
    uint32_t endianMarker;
    TokamakConfigWire config;
    NBIConfigWire nbi;
    uint64_t createdUnixNs;
    uint32_t reserved[8];
};

struct SnapshotFrameHeader {
    uint64_t stepIndex;
    double simTimeSeconds;

    uint32_t totalParticleCount;
    uint32_t sampledParticleCount;
    uint32_t sampleStride;
    uint32_t reserved0;

    uint32_t deuteriumCount;
    uint32_t tritiumCount;
    uint32_t heliumCount;
    uint32_t totalIons;

    double avgTempKeV;
    uint64_t fusionEvents;

    double simStepMs;
    double exportMs;
    double frameSubmitMs;
    uint64_t reserved1;
};

struct ParticleSampleWire {
    float x;
    float y;
    float z;
    float vx;
    float vy;
    float vz;
    float mass;
    float charge;
    float qOverM;
    uint32_t species;
};

TokamakConfigWire ToWire(const TokamakConfig& config) {
    TokamakConfigWire wire{};
    wire.majorRadius = config.majorRadius;
    wire.minorRadius = config.minorRadius;
    wire.toroidalCurrent = config.toroidalCurrent;
    wire.toroidalCoilTurns = config.toroidalCoilTurns;
    wire.plasmaCurrent = config.plasmaCurrent;
    return wire;
}

NBIConfigWire ToWire(const NBIConfig& nbi) {
    NBIConfigWire wire{};
    wire.isActive = nbi.isActive ? 1u : 0u;
    wire.beamEnergyKeV = nbi.beamEnergyKeV;
    wire.particlesPerStep = nbi.particlesPerStep;
    wire.injectorPos[0] = nbi.injectorPos.x;
    wire.injectorPos[1] = nbi.injectorPos.y;
    wire.injectorPos[2] = nbi.injectorPos.z;
    wire.injectionNormal[0] = nbi.injectionNormal.x;
    wire.injectionNormal[1] = nbi.injectionNormal.y;
    wire.injectionNormal[2] = nbi.injectionNormal.z;
    return wire;
}

TokamakConfig FromWire(const TokamakConfigWire& wire) {
    TokamakConfig out;
    out.majorRadius = wire.majorRadius;
    out.minorRadius = wire.minorRadius;
    out.toroidalCurrent = wire.toroidalCurrent;
    out.toroidalCoilTurns = wire.toroidalCoilTurns;
    out.plasmaCurrent = wire.plasmaCurrent;
    return out;
}

NBIConfig FromWire(const NBIConfigWire& wire) {
    NBIConfig out;
    out.isActive = wire.isActive != 0;
    out.beamEnergyKeV = wire.beamEnergyKeV;
    out.particlesPerStep = wire.particlesPerStep;
    out.injectorPos = Vec3(wire.injectorPos[0], wire.injectorPos[1], wire.injectorPos[2]);
    out.injectionNormal = Vec3(
        wire.injectionNormal[0],
        wire.injectionNormal[1],
        wire.injectionNormal[2]);
    return out;
}

ParticleSampleWire ToWire(const ParticleSample& sample) {
    ParticleSampleWire wire{};
    wire.x = sample.x;
    wire.y = sample.y;
    wire.z = sample.z;
    wire.vx = sample.vx;
    wire.vy = sample.vy;
    wire.vz = sample.vz;
    wire.mass = sample.mass;
    wire.charge = sample.charge;
    wire.qOverM = sample.qOverM;
    wire.species = static_cast<uint32_t>(sample.species);
    return wire;
}

ParticleSample FromWire(const ParticleSampleWire& wire) {
    ParticleSample sample;
    sample.x = wire.x;
    sample.y = wire.y;
    sample.z = wire.z;
    sample.vx = wire.vx;
    sample.vy = wire.vy;
    sample.vz = wire.vz;
    sample.mass = wire.mass;
    sample.charge = wire.charge;
    sample.qOverM = wire.qOverM;
    sample.species = static_cast<uint8_t>(wire.species);
    return sample;
}

template <typename T>
bool WriteBinary(std::ofstream& file, const T& value) {
    file.write(reinterpret_cast<const char*>(&value), sizeof(T));
    return file.good();
}

template <typename T>
bool ReadBinary(std::ifstream& file, T& value) {
    file.read(reinterpret_cast<char*>(&value), sizeof(T));
    return file.good();
}

uint64_t UnixNowNs() {
    const auto now = std::chrono::system_clock::now().time_since_epoch();
    return static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(now).count());
}

}  // namespace

SnapshotStreamWriter::~SnapshotStreamWriter() {
    Close();
}

bool SnapshotStreamWriter::Open(const std::string& path) {
    Close();
    lastError_.clear();
    file_.open(path, std::ios::binary | std::ios::trunc);
    headerWritten_ = false;
    if (!file_.is_open()) {
        SetError("Failed to open snapshot file for writing: " + path);
        return false;
    }
    return true;
}

bool SnapshotStreamWriter::WriteHeader(const SimulationSnapshot& snapshot) {
    SnapshotFileHeader header{};
    std::memcpy(header.magic, kMagic, sizeof(kMagic));
    header.version = kSnapshotVersion;
    header.endianMarker = kEndianMarker;
    header.config = ToWire(snapshot.config);
    header.nbi = ToWire(snapshot.nbi);
    header.createdUnixNs = UnixNowNs();

    if (!WriteBinary(file_, header)) {
        SetError("Failed to write snapshot file header");
        return false;
    }

    headerWritten_ = true;
    return true;
}

bool SnapshotStreamWriter::Write(const SimulationSnapshot& snapshot) {
    if (!file_.is_open()) {
        SetError("Snapshot writer is not open");
        return false;
    }

    if (!headerWritten_ && !WriteHeader(snapshot)) {
        return false;
    }

    SnapshotFrameHeader frame{};
    frame.stepIndex = snapshot.stepIndex;
    frame.simTimeSeconds = snapshot.simTimeSeconds;
    frame.totalParticleCount = snapshot.totalParticleCount;
    frame.sampledParticleCount = static_cast<uint32_t>(snapshot.sampledParticles.size());
    frame.sampleStride = snapshot.sampleStride;
    frame.deuteriumCount = snapshot.telemetry.deuteriumCount;
    frame.tritiumCount = snapshot.telemetry.tritiumCount;
    frame.heliumCount = snapshot.telemetry.heliumCount;
    frame.totalIons = snapshot.telemetry.totalIons;
    frame.avgTempKeV = snapshot.telemetry.avgTempKeV;
    frame.fusionEvents = snapshot.telemetry.fusionEvents;
    frame.simStepMs = snapshot.performance.simStepMs;
    frame.exportMs = snapshot.performance.exportMs;
    frame.frameSubmitMs = snapshot.performance.frameSubmitMs;

    if (!WriteBinary(file_, frame)) {
        SetError("Failed to write snapshot frame header");
        return false;
    }

    for (const ParticleSample& sample : snapshot.sampledParticles) {
        const ParticleSampleWire wire = ToWire(sample);
        if (!WriteBinary(file_, wire)) {
            SetError("Failed to write particle sample");
            return false;
        }
    }

    file_.flush();
    if (!file_.good()) {
        SetError("Failed to flush snapshot frame");
        return false;
    }

    return true;
}

void SnapshotStreamWriter::Close() {
    if (file_.is_open()) {
        file_.close();
    }
    headerWritten_ = false;
}

void SnapshotStreamWriter::SetError(const std::string& message) {
    lastError_ = message;
}

SnapshotStreamReader::~SnapshotStreamReader() {
    Close();
}

bool SnapshotStreamReader::Open(const std::string& path) {
    Close();
    lastError_.clear();

    path_ = path;
    file_.open(path, std::ios::binary);
    headerRead_ = false;
    reachedEof_ = false;

    if (!file_.is_open()) {
        SetError("Failed to open snapshot file for reading: " + path);
        return false;
    }

    return ReadHeader();
}

bool SnapshotStreamReader::ReadHeader() {
    SnapshotFileHeader header{};
    if (!ReadBinary(file_, header)) {
        SetError("Failed to read snapshot file header");
        return false;
    }

    if (std::memcmp(header.magic, kMagic, sizeof(kMagic)) != 0) {
        SetError("Invalid snapshot magic");
        return false;
    }

    if (header.version != kSnapshotVersion) {
        SetError("Unsupported snapshot version");
        return false;
    }

    if (header.endianMarker != kEndianMarker) {
        SetError("Snapshot endianness mismatch");
        return false;
    }

    config_ = FromWire(header.config);
    nbi_ = FromWire(header.nbi);
    headerRead_ = true;
    return true;
}

bool SnapshotStreamReader::ReadNext(SimulationSnapshot& snapshot) {
    if (!file_.is_open()) {
        SetError("Snapshot reader is not open");
        return false;
    }

    if (!headerRead_ && !ReadHeader()) {
        return false;
    }

    lastError_.clear();
    reachedEof_ = false;
    SnapshotFrameHeader frame{};
    if (!ReadBinary(file_, frame)) {
        if (file_.eof()) {
            reachedEof_ = true;
            return false;
        }
        SetError("Failed to read snapshot frame header");
        return false;
    }

    snapshot.stepIndex = frame.stepIndex;
    snapshot.simTimeSeconds = frame.simTimeSeconds;
    snapshot.config = config_;
    snapshot.nbi = nbi_;

    snapshot.totalParticleCount = frame.totalParticleCount;
    snapshot.sampleStride = frame.sampleStride;
    snapshot.telemetry.deuteriumCount = frame.deuteriumCount;
    snapshot.telemetry.tritiumCount = frame.tritiumCount;
    snapshot.telemetry.heliumCount = frame.heliumCount;
    snapshot.telemetry.totalIons = frame.totalIons;
    snapshot.telemetry.avgTempKeV = frame.avgTempKeV;
    snapshot.telemetry.fusionEvents = frame.fusionEvents;

    snapshot.performance.simStepMs = frame.simStepMs;
    snapshot.performance.exportMs = frame.exportMs;
    snapshot.performance.frameSubmitMs = frame.frameSubmitMs;

    snapshot.sampledParticles.resize(frame.sampledParticleCount);
    for (uint32_t i = 0; i < frame.sampledParticleCount; ++i) {
        ParticleSampleWire wire{};
        if (!ReadBinary(file_, wire)) {
            SetError("Failed to read particle sample payload");
            return false;
        }
        snapshot.sampledParticles[i] = FromWire(wire);
    }

    return true;
}

bool SnapshotStreamReader::Rewind() {
    if (path_.empty()) {
        SetError("Cannot rewind before opening a snapshot file");
        return false;
    }
    return Open(path_);
}

void SnapshotStreamReader::Close() {
    if (file_.is_open()) {
        file_.close();
    }
    headerRead_ = false;
    reachedEof_ = false;
    path_.clear();
}

void SnapshotStreamReader::SetError(const std::string& message) {
    lastError_ = message;
}

bool ExportSnapshotFileToCsv(
    const std::string& snapshotPath,
    const std::string& csvPrefix,
    size_t particleRowsPerFrame,
    std::string* errorOut) {
    SnapshotStreamReader reader;
    if (!reader.Open(snapshotPath)) {
        if (errorOut != nullptr) {
            *errorOut = reader.LastError();
        }
        return false;
    }

    const std::string aggregatePath = csvPrefix + "_aggregates.csv";
    const std::string particlesPath = csvPrefix + "_particles.csv";

    std::ofstream aggregateFile(aggregatePath, std::ios::trunc);
    std::ofstream particlesFile(particlesPath, std::ios::trunc);

    if (!aggregateFile.is_open() || !particlesFile.is_open()) {
        if (errorOut != nullptr) {
            *errorOut = "Failed to open CSV output files";
        }
        return false;
    }

    aggregateFile
        << "step,sim_time_s,total_particles,d_count,t_count,he_count,avg_temp_kev,fusion_events,sim_step_ms,export_ms\n";
    particlesFile
        << "step,sample_index,species,x,y,z,vx,vy,vz,mass,charge,q_over_m\n";

    SimulationSnapshot snapshot;
    uint64_t frameCount = 0;
    while (reader.ReadNext(snapshot)) {
        aggregateFile << snapshot.stepIndex << ','
                      << std::setprecision(std::numeric_limits<double>::max_digits10)
                      << snapshot.simTimeSeconds << ','
                      << snapshot.totalParticleCount << ','
                      << snapshot.telemetry.deuteriumCount << ','
                      << snapshot.telemetry.tritiumCount << ','
                      << snapshot.telemetry.heliumCount << ','
                      << snapshot.telemetry.avgTempKeV << ','
                      << snapshot.telemetry.fusionEvents << ','
                      << snapshot.performance.simStepMs << ','
                      << snapshot.performance.exportMs << '\n';

        const size_t rowCount = std::min(particleRowsPerFrame, snapshot.sampledParticles.size());
        for (size_t i = 0; i < rowCount; ++i) {
            const ParticleSample& p = snapshot.sampledParticles[i];
            particlesFile << snapshot.stepIndex << ','
                          << i << ','
                          << static_cast<uint32_t>(p.species) << ','
                          << p.x << ',' << p.y << ',' << p.z << ','
                          << p.vx << ',' << p.vy << ',' << p.vz << ','
                          << p.mass << ',' << p.charge << ',' << p.qOverM << '\n';
        }

        ++frameCount;
    }

    if (frameCount == 0) {
        if (errorOut != nullptr) {
            *errorOut = "No frames were read from snapshot file";
        }
        return false;
    }

    return true;
}
