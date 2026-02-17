#include "profiler.h"

#ifdef GENIE_PROFILE

#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

// Thread-local storage for active timers and timer stack
thread_local std::unordered_map<std::string, Profiler::TimerState> Profiler::thread_timers_;
thread_local std::vector<std::string> Profiler::timer_stack_;

std::string Profiler::build_path(const std::vector<std::string>& stack) {
    std::string path;
    for (size_t i = 0; i < stack.size(); ++i) {
        if (i > 0) path += ':';
        path += stack[i];
    }
    return path;
}

void Profiler::start(const char* name) {
    // Push leaf name onto stack and build full path
    timer_stack_.push_back(std::string(name));
    std::string full_path = build_path(timer_stack_);

    auto& state = thread_timers_[full_path];
    state.start_time = std::chrono::steady_clock::now();
    state.active = true;
}

void Profiler::stop(const char* name) {
    auto stop_time = std::chrono::steady_clock::now();

    // Build the full path from the current stack (should end with name)
    std::string full_path = build_path(timer_stack_);

    auto it = thread_timers_.find(full_path);
    if (it == thread_timers_.end() || !it->second.active) {
        // Timer wasn't started or already stopped; still pop the stack
        if (!timer_stack_.empty()) timer_stack_.pop_back();
        return;
    }

    // Calculate elapsed time
    double elapsed = std::chrono::duration<double>(stop_time - it->second.start_time).count();
    it->second.active = false;

    // Pop the leaf name from the stack
    timer_stack_.pop_back();

    // Build parent's full path (the stack after popping)
    std::string parent_path = build_path(timer_stack_);

    // Aggregate into shared entries (thread-safe)
    {
        std::lock_guard<std::mutex> lock(mutex_);
        auto& entry = entries_[full_path];
        entry.name = full_path;
        entry.calls++;
        entry.total_seconds += elapsed;

        // Add elapsed time to parent's children_seconds
        if (!parent_path.empty()) {
            auto& parent_entry = entries_[parent_path];
            parent_entry.name = parent_path;
            parent_entry.children_seconds += elapsed;
        }
    }
}

void Profiler::add_bytes(const char* name, uint64_t bytes) {
    // Build full path from current stack + leaf name if the leaf matches the
    // top of the stack (i.e., we're inside that timer). Otherwise, build path
    // using the stack as-is plus the provided name.
    // Since add_bytes is called while the timer is active, the leaf should
    // already be on the stack. Use the current full path.
    std::string key;
    if (!timer_stack_.empty() && timer_stack_.back() == std::string(name)) {
        key = build_path(timer_stack_);
    } else {
        // Leaf not on stack — build path with stack + leaf
        std::vector<std::string> tmp = timer_stack_;
        tmp.push_back(std::string(name));
        key = build_path(tmp);
    }
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

    // Sort: "_" first, then by self_seconds descending
    std::sort(result.begin(), result.end(),
              [](const ProfileEntry& a, const ProfileEntry& b) {
                  if (a.name == "_") return true;
                  if (b.name == "_") return false;
                  return a.self_seconds() > b.self_seconds();
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

    if (is_json) {
        // JSON format
        out << "[\n";
        for (size_t i = 0; i < entries.size(); ++i) {
            const auto& entry = entries[i];
            out << "  {\n";
            out << "    \"name\": \"" << entry.name << "\",\n";
            out << "    \"calls\": " << entry.calls << ",\n";
            out << "    \"total_seconds\": " << std::fixed << std::setprecision(6) << entry.total_seconds << ",\n";
            out << "    \"self_seconds\": " << std::fixed << std::setprecision(6) << entry.self_seconds() << ",\n";
            out << "    \"children_seconds\": " << std::fixed << std::setprecision(6) << entry.children_seconds;
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
        out << "name,calls,total_seconds,self_seconds,children_seconds,avg_self_seconds,total_bytes\n";
        for (const auto& entry : entries) {
            double avg_self = (entry.calls > 0) ? (entry.self_seconds() / entry.calls) : 0.0;
            out << entry.name << ","
                << entry.calls << ","
                << std::fixed << std::setprecision(6) << entry.total_seconds << ","
                << std::fixed << std::setprecision(6) << entry.self_seconds() << ","
                << std::fixed << std::setprecision(6) << entry.children_seconds << ","
                << std::fixed << std::setprecision(6) << avg_self << ","
                << entry.total_bytes << "\n";
        }
    }

    out.close();
    std::cout << "Profiling data written to: " << filename << std::endl;
}

#endif  // GENIE_PROFILE
