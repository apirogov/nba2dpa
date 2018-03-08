#pragma once

#include <memory>
#include <queue>
#include <set>
#include <vector>
#include "level.hh"
#include "common/scc.hh"
#include "swa.hh"

namespace nbautils {

using PA = SWA<Level>;

PA::uptr determinize(LevelConfig const& lu);

PA::uptr determinize(LevelConfig const& lu, PS const& psa, SCCDat<state_t> const& psai);

}  // namespace nbautils
