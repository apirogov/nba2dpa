#pragma once

#include <utility>
#include <string>
#include <iostream>
#include "ba.hh"
#include "det.hh"

#include <spdlog/spdlog.h>

namespace nbautils {

std::vector<BA::uptr> parse_hoa_ba(std::string const& filename, std::shared_ptr<spdlog::logger> log=nullptr);

std::string sym_to_edgelabel(sym_t s, std::vector<std::string> const& aps, bool as_aps=false);

void print_hoa_pa(PA const& aut, ostream &out = cout);

}  // namespace nbautils
