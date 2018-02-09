#pragma once

#include <memory>
#include <queue>
#include <set>
#include <vector>
#include "level.hh"
#include "scc.hh"
#include "types.hh"

namespace nbautils {

using pa_tag_store = naive_bimap<Level, state_t>;
using PA = SWA<priority_t, Level, pa_tag_store>;

PA::uptr determinize(LevelConfig const& lu);

}  // namespace nbautils
