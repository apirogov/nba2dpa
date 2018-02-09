#pragma once

#include "types.hh"
#include "scc.hh"
#include "level.hh"
#include <memory>
#include <set>
#include <vector>
#include <queue>

namespace nbautils {

using pa_tag_store = naive_bimap<Level,state_t>;
using PA = SWA<priority_t,Level,pa_tag_store>;

PA::uptr determinize(LevelConfig const& lu);

}
