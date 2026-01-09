#include "profiler.h"

#ifdef GENIE_PROFILE

#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

// Thread-local storage for active timers
thread_local std::unordered_map<std::string, Profiler::TimerState> Profiler::thread_timers_;

void Profiler::start(const char* name) {
    std::string key(name);
    auto& state = thread_timers_[key];
    state.start_time = std::chrono::steady_clock::now();
    state.active = true;
}

void Profiler::stop(const char* name) {
    auto stop_time = std::chrono::steady_clock::now();
    std::string key(name);

    auto it = thread_timers_.find(key);
    if (it == thread_timers_.end() || !it->second.active) {
        // Timer wasn't started or already stopped
        return;
    }

    // Calculate elapsed time
    auto elapsed = std::chrono::duration<double>(stop_time - it->second.start_time).count();
    it->second.active = false;

    // Aggregate into shared entries (thread-safe)
    {
        std::lock_guard<std::mutex> lock(mutex_);
        auto& entry = entries_[key];
        entry.name = key;
        entry.calls++;
        entry.total_seconds += elapsed;
    }
}

void Profiler::add_bytes(const char* name, uint64_t bytes) {
    std::string key(name);
    std::lock_guard<std::mutex> lock(mutex_);
    auto& entry = entries_[key];
    entry.name = key;
    entry.total_bytes += bytes;
}

std::vector<ProfileEntry> Profiler::entries() const {
    std::lock_guard<std::mutex> lock(mutex_);
    std::vector<ProfileEntry> result;
    result.reserve(entries_.size());

    for (const auto& pair : entries_) {
        result.push_back(pair.second);
    }

    // Sort by total time descending
    std::sort(result.begin(), result.end(),
              [](const ProfileEntry& a, const ProfileEntry& b) {
                  return a.total_seconds > b.total_seconds;
              });

    return result;
}

void Profiler::clear() {
    std::lock_guard<std::mutex> lock(mutex_);
    entries_.clear();
}

void dump_profile(const std::vector<ProfileEntry>& entries, const std::string& filename) {
    if (entries.empty()) {
        std::cerr << "Warning: No profiling data collected" << std::endl;
        return;
    }

    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: Could not open profile output file: " << filename << std::endl;
        return;
    }

    // Determine format from extension
    bool is_json = (filename.size() >= 5 && filename.substr(filename.size() - 5) == ".json");
    bool is_csv = (filename.size() >= 4 && filename.substr(filename.size() - 4) == ".csv");

    if (is_json) {
        // JSON format
        out << "[\n";
        for (size_t i = 0; i < entries.size(); ++i) {
            const auto& entry = entries[i];
            out << "  {\n";
            out << "    \"name\": \"" << entry.name << "\",\n";
            out << "    \"calls\": " << entry.calls << ",\n";
            out << "    \"total_seconds\": " << std::fixed << std::setprecision(6) << entry.total_seconds;
            if (entry.total_bytes > 0) {
                out << ",\n    \"total_bytes\": " << entry.total_bytes;
            }
            out << "\n  }";
            if (i + 1 < entries.size()) {
                out << ",";
            }
            out << "\n";
        }
        out << "]\n";
    } else {
        // CSV format (default)
        out << "name,calls,total_seconds,avg_seconds,total_bytes\n";
        for (const auto& entry : entries) {
            double avg = (entry.calls > 0) ? (entry.total_seconds / entry.calls) : 0.0;
            out << entry.name << ","
                << entry.calls << ","
                << std::fixed << std::setprecision(6) << entry.total_seconds << ","
                << std::fixed << std::setprecision(6) << avg << ","
                << entry.total_bytes << "\n";
        }
    }

    out.close();
    std::cout << "Profiling data written to: " << filename << std::endl;
}

#endif  // GENIE_PROFILE
