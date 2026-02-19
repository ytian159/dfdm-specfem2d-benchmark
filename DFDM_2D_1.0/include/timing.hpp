#include <chrono>

// Function to start timing
std::chrono::time_point<std::chrono::high_resolution_clock> start_timing() {
    return std::chrono::high_resolution_clock::now();
}

// Function to stop timing and return duration in milliseconds
long long stop_timing_ms(std::chrono::time_point<std::chrono::high_resolution_clock> start) {
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    return duration.count();
}

// Function to stop timing and return duration in microseconds
long long stop_timing_us(std::chrono::time_point<std::chrono::high_resolution_clock> start) {
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    return duration.count();
}

// Function to stop timing and return duration in nanoseconds
long long stop_timing_ns(std::chrono::time_point<std::chrono::high_resolution_clock> start) {
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    return duration.count();
}