#include <atomic>
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <vector>

#include <omp.h>

// Minimal, algorithm-independent probe for GCC/libgomp ThreadSanitizer
// recognition of the implicit synchronization at an OpenMP region boundary.
int main()
{
    struct SharedState {
        std::atomic<int> visits{0};
        std::mutex mutex;
        std::vector<int> per_thread;
    } shared;

    shared.per_thread.resize(static_cast<std::size_t>(omp_get_max_threads()));

    // Repeated regions and a stack-captured worker lambda mirror the launcher
    // and lifetime shape used by Stage-4B without including its algorithm.
    for (int phase = 0; phase < 8; ++phase) {
        std::atomic<int> phase_visits{0};
        auto worker = [&](int thread) {
            std::vector<int> local(16, thread + phase);
            phase_visits.fetch_add(1, std::memory_order_seq_cst);
            shared.visits.fetch_add(1, std::memory_order_seq_cst);
            std::lock_guard<std::mutex> guard(shared.mutex);
            shared.per_thread[static_cast<std::size_t>(thread)] += local[0];
        };

#pragma omp parallel shared(shared, phase_visits, worker)
        { worker(omp_get_thread_num()); }

        if (phase_visits.load(std::memory_order_seq_cst) !=
            omp_get_max_threads()) {
            std::cerr << "OpenMP phase invariant failed\n";
            return EXIT_FAILURE;
        }
    }

    // The implicit barrier/end of the parallel region must happen-before these
    // main-thread reads and the subsequent destruction of shared.
    int sum = 0;
    {
        std::lock_guard<std::mutex> guard(shared.mutex);
        for (int value : shared.per_thread) sum += value;
    }
    if (shared.visits.load(std::memory_order_seq_cst) !=
            8 * omp_get_max_threads() ||
        sum <= 0) {
        std::cerr << "OpenMP runtime probe invariant failed\n";
        return EXIT_FAILURE;
    }
    std::cout << "OpenMP runtime probe passed\n";
    return EXIT_SUCCESS;
}
