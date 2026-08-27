#pragma once

#include <chrono>

using Clock = std::chrono::steady_clock;

inline double elapsed_seconds(Clock::time_point begin, Clock::time_point end)
{
    return std::chrono::duration<double>(end - begin).count();
}
