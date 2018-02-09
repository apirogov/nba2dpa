#pragma once

#include "types.hh"
#include <utility>

namespace nbautils {

// this should be "optional" type
BA::uptr parse_ba(std::string const& filename);

// TODO: HOA printer
}
