#pragma once

#include <utility>
#include "types.hh"

namespace nbautils {

// this should be "optional" type
BA::uptr parse_ba(std::string const& filename);

// TODO: HOA printer
}  // namespace nbautils
