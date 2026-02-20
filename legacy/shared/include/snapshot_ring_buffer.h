#pragma once

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "sim_snapshot.h"

// Single-producer/single-consumer ring buffer.
// Producer never blocks and drops new snapshots when full.
class SnapshotRingBuffer {
public:
    SnapshotRingBuffer(size_t capacity, size_t reserveParticlesPerSlot)
        : capacity_(capacity == 0 ? 1 : capacity), slots_(capacity_) {
        for (SimulationSnapshot& slot : slots_) {
            slot.Reserve(reserveParticlesPerSlot);
        }
    }

    bool Publish(const SimulationSnapshot& snapshot) {
        const size_t write = writeIndex_.load(std::memory_order_relaxed);
        const size_t next = (write + 1) % capacity_;
        const size_t read = readIndex_.load(std::memory_order_acquire);

        if (next == read) {
            dropped_.fetch_add(1, std::memory_order_relaxed);
            return false;
        }

        slots_[write] = snapshot;
        writeIndex_.store(next, std::memory_order_release);
        produced_.fetch_add(1, std::memory_order_relaxed);
        return true;
    }

    bool ConsumeLatest(SimulationSnapshot& out) {
        const size_t read = readIndex_.load(std::memory_order_relaxed);
        const size_t write = writeIndex_.load(std::memory_order_acquire);

        if (read == write) {
            return false;
        }

        const size_t latest = (write + capacity_ - 1) % capacity_;
        out = slots_[latest];
        readIndex_.store(write, std::memory_order_release);

        consumed_.fetch_add(1, std::memory_order_relaxed);
        return true;
    }

    size_t PendingCount() const {
        const size_t read = readIndex_.load(std::memory_order_acquire);
        const size_t write = writeIndex_.load(std::memory_order_acquire);
        return (write + capacity_ - read) % capacity_;
    }

    uint64_t ProducedCount() const { return produced_.load(std::memory_order_relaxed); }
    uint64_t ConsumedCount() const { return consumed_.load(std::memory_order_relaxed); }
    uint64_t DroppedCount() const { return dropped_.load(std::memory_order_relaxed); }

private:
    const size_t capacity_;
    std::vector<SimulationSnapshot> slots_;

    std::atomic<size_t> writeIndex_{0};
    std::atomic<size_t> readIndex_{0};

    std::atomic<uint64_t> produced_{0};
    std::atomic<uint64_t> consumed_{0};
    std::atomic<uint64_t> dropped_{0};
};
