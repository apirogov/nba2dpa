#pragma once

#include "types.hh"
#include <utility>

namespace nbautils {
// this should be "optional" type
std::pair<BA::ptr,ParsedMeta> parse_ba(std::string const& filename);

// TODO: HOA printer
}
