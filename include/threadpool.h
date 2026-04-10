#ifndef GENIE_THREADPOOL_H
#define GENIE_THREADPOOL_H

#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <atomic>

// Simple fork-join thread pool for GENIE's parallel mailman operations.
// Threads are created once and reused across all multiply calls.
// Uses a generation counter to prevent spurious re-execution.
class ForkJoinPool {
public:
    explicit ForkJoinPool(int num_threads)
        : stop_(false), generation_(0), num_tasks_(0), completed_(0) {
        for (int i = 0; i < num_threads; i++) {
            workers_.emplace_back([this, i] { worker_loop(i); });
        }
    }

    ~ForkJoinPool() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            stop_ = true;
        }
        cv_work_.notify_all();
        for (auto &w : workers_) {
            if (w.joinable()) w.join();
        }
    }

    ForkJoinPool(const ForkJoinPool&) = delete;
    ForkJoinPool& operator=(const ForkJoinPool&) = delete;

    // Dispatch n_tasks parallel tasks and wait for all to complete.
    // Each task receives its index (0..n_tasks-1).
    void parallel_for(int n_tasks, const std::function<void(int)> &task) {
        if (n_tasks <= 0) return;
        if (n_tasks == 1) {
            task(0);
            return;
        }

        int actual = (n_tasks <= static_cast<int>(workers_.size()))
                     ? n_tasks : static_cast<int>(workers_.size());

        {
            std::lock_guard<std::mutex> lock(mutex_);
            current_task_ = &task;
            num_tasks_ = actual;
            completed_.store(0, std::memory_order_release);
            generation_++;
        }
        cv_work_.notify_all();

        // Wait for all tasks to complete
        {
            std::unique_lock<std::mutex> lock(mutex_);
            cv_done_.wait(lock, [this, actual] {
                return completed_.load(std::memory_order_acquire) >= actual;
            });
        }
    }

    int size() const { return static_cast<int>(workers_.size()); }

private:
    void worker_loop(int id) {
        unsigned long long last_gen = 0;
        while (true) {
            unsigned long long gen;
            const std::function<void(int)> *task;
            int ntasks;
            {
                std::unique_lock<std::mutex> lock(mutex_);
                cv_work_.wait(lock, [this, &last_gen] {
                    return stop_ || generation_ > last_gen;
                });
                if (stop_) return;
                gen = generation_;
                task = current_task_;
                ntasks = num_tasks_;
            }
            last_gen = gen;

            // Only execute if this worker's id is within range
            if (id < ntasks) {
                (*task)(id);
                // Signal completion
                if (completed_.fetch_add(1, std::memory_order_acq_rel) + 1 >= ntasks) {
                    cv_done_.notify_one();
                }
            }
        }
    }

    std::vector<std::thread> workers_;
    std::mutex mutex_;
    std::condition_variable cv_work_;
    std::condition_variable cv_done_;
    bool stop_;
    unsigned long long generation_;
    int num_tasks_;
    std::atomic<int> completed_;
    const std::function<void(int)> *current_task_ = nullptr;
};

#endif // GENIE_THREADPOOL_H
