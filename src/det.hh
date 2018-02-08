#pragma once

#include "types.hh"
#include <memory>
#include <set>
#include <vector>
#include <queue>

namespace nbautils {

PA::ptr determinize(BA const& aut, SCCInfo const& scci);

}
