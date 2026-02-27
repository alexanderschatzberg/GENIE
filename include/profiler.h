#ifndef PROFILER_H
#define PROFILER_H

#include <string>
#include <vector>
#include <cstdint>

// ProfileEntry holds aggregated timing information for a named region
struct ProfileEntry {
    std::string name;
    uint64_t calls;
    double total_seconds;
    double children_seconds;  // time spent in direct children
    uint64_t total_bytes;     // Optional: for I/O tracking

    ProfileEntry() : calls(0), total_seconds(0.0), children_seconds(0.0), total_bytes(0) {}
    ProfileEntry(const std::string& n, uint64_t c, double t, double ch = 0.0, uint64_t b = 0)
        : name(n), calls(c), total_seconds(t), children_seconds(ch), total_bytes(b) {}

    double self_seconds() const { return total_seconds - children_seconds; }
};

#ifdef GENIE_PROFILE

#include <unordered_map>
#include <mutex>
#include <chrono>

// Profiler singleton - manages timing data collection
class Profiler {
public:
    static Profiler& instance() {
        static Profiler profiler;
        return profiler;
    }

    // Start timing a named region (leaf name; full path built from stack)
    void start(const char* name);

    // Stop timing the current region (leaf name used for validation)
    void stop(const char* name);

    // Add bytes transferred (for I/O tracking, uses current stack + leaf)
    void add_bytes(const char* name, uint64_t bytes);

    // Get all profile entries
    std::vector<ProfileEntry> entries() const;

    // Clear all profiling data
    void clear();

private:
    Profiler() {}
    ~Profiler() {}

    // Prevent copying
    Profiler(const Profiler&) = delete;
    Profiler& operator=(const Profiler&) = delete;

    // Thread-local storage for start times (keyed by full path)
    struct TimerState {
        std::chrono::time_point<std::chrono::steady_clock> start_time;
        bool active;
        TimerState() : active(false) {}
    };

    // Build colon-separated path from stack
    static std::string build_path(const std::vector<std::string>& stack);

    // Map from region name to aggregated stats
    mutable std::mutex mutex_;
    std::unordered_map<std::string, ProfileEntry> entries_;

    // Thread-local timer stack (leaf names) and start-time map (keyed by full path)
    static thread_local std::vector<std::string> timer_stack_;
    static thread_local std::unordered_map<std::string, TimerState> thread_timers_;
};

// RAII timer - automatically starts on construction, stops on destruction
struct ScopedTimer {
    const char* name;

    explicit ScopedTimer(const char* n) : name(n) {
        Profiler::instance().start(name);
    }

    ~ScopedTimer() {
        Profiler::instance().stop(name);
    }

    // Prevent copying
    ScopedTimer(const ScopedTimer&) = delete;
    ScopedTimer& operator=(const ScopedTimer&) = delete;
};

// Function to dump profile data to a file (JSON or CSV)
void dump_profile(const std::vector<ProfileEntry>& entries, const std::string& filename);

#else  // GENIE_PROFILE not defined - all profiling is no-ops

// Dummy implementations when profiling is disabled
struct ScopedTimer {
    explicit ScopedTimer(const char*) {}
};

class Profiler {
public:
    static Profiler& instance() {
        static Profiler p;
        return p;
    }

    void start(const char*) {}
    void stop(const char*) {}
    void add_bytes(const char*, uint64_t) {}
    std::vector<ProfileEntry> entries() const { return std::vector<ProfileEntry>(); }
    void clear() {}
};

inline void dump_profile(const std::vector<ProfileEntry>&, const std::string&) {}

#endif  // GENIE_PROFILE

#endif  // PROFILER_H
