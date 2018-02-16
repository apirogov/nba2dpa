#pragma once
#include <chrono>
#include <string>

#include <spdlog/spdlog.h>

using timepoint_t = std::chrono::high_resolution_clock::time_point;
using duration_t = std::chrono::high_resolution_clock::duration;

timepoint_t get_time();
double duration_to_sec(duration_t const& tp);
double get_secs_since(timepoint_t const& tp);

//any function can be benchmarked like this:
//  bench(logger, "function name", WRAP(function_call(args)));
#define WRAP(x) ([&](){return std::move(x);})
template<typename F>
auto bench(std::shared_ptr<spdlog::logger> log, std::string name, F f, bool enabled=true) {
  timepoint_t starttime;
  if (log && enabled) {
    starttime = get_time();
    log->info("bench: {} started", name);
  }

  auto r = f();

  if (log && enabled) {
    log->info("bench: {} completed ({:.4f} s)"
              , name, get_secs_since(starttime));
  }

  return std::move(r);
}
