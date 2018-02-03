#include "nba2dpa.hpp"
#include <spdlog/spdlog.h>
#include <iostream>
#include <memory>

namespace nba2dpa {
int times(int a, int b) {
  static auto logger = spdlog::stdout_color_mt("times_logger");
  logger->debug("times called with ({}, {})", a, b);
  int result = a * b;
  logger->debug("Result is {}", result);
  return result;
}
}  // namespace nba2dpa

int main(int argc, char *argv[]) { std::cout << "hello world!" << std::endl; }
