#pragma once
#include <chrono>

using timepoint_t = std::chrono::high_resolution_clock::time_point;
using duration_t = std::chrono::high_resolution_clock::duration;

timepoint_t get_time();
double duration_to_sec(duration_t const& tp);
double get_secs_since(timepoint_t const& tp);
