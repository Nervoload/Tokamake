#pragma once

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <string>

#include "sim_snapshot.h"

class SnapshotStreamWriter {
public:
    SnapshotStreamWriter() = default;
    ~SnapshotStreamWriter();

    bool Open(const std::string& path);
    bool Write(const SimulationSnapshot& snapshot);
    void Close();

    bool IsOpen() const { return file_.is_open(); }
    const std::string& LastError() const { return lastError_; }

private:
    bool WriteHeader(const SimulationSnapshot& snapshot);
    void SetError(const std::string& message);

    std::ofstream file_;
    std::string lastError_;
    bool headerWritten_ = false;
};

class SnapshotStreamReader {
public:
    SnapshotStreamReader() = default;
    ~SnapshotStreamReader();

    bool Open(const std::string& path);
    bool ReadNext(SimulationSnapshot& snapshot);
    bool Rewind();
    void Close();

    bool IsOpen() const { return file_.is_open(); }
    bool ReachedEof() const { return reachedEof_; }
    const std::string& LastError() const { return lastError_; }

private:
    bool ReadHeader();
    void SetError(const std::string& message);

    std::ifstream file_;
    std::string path_;
    std::string lastError_;

    TokamakConfig config_;
    NBIConfig nbi_;
    bool headerRead_ = false;
    bool reachedEof_ = false;
};

bool ExportSnapshotFileToCsv(
    const std::string& snapshotPath,
    const std::string& csvPrefix,
    size_t particleRowsPerFrame,
    std::string* errorOut = nullptr);
