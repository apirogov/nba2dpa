#pragma once

#include <utility>
#include <iostream>
#include "types.hh"
#include "det.hh"

#include <spdlog/spdlog.h>

namespace nbautils {

BA::uptr parse_ba(std::string const& filename, std::shared_ptr<spdlog::logger> log=nullptr);

void print_hoa_pa(PA const& aut, ostream &out = cout);

}  // namespace nbautils
