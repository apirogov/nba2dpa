#include "bench.hh"

timepoint_t get_time() { return std::chrono::high_resolution_clock::now(); }

double duration_to_sec(duration_t const& tp) {
  return std::chrono::duration_cast<std::chrono::duration<double>>(tp).count();
}

double get_secs_since(timepoint_t const& tp) {
  return duration_to_sec(get_time()-tp);
}

