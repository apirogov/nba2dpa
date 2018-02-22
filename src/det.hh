#pragma once

#include <memory>
#include <queue>
#include <set>
#include <vector>
#include "level.hh"
#include "scc.hh"
#include "types.hh"

namespace nbautils {

using PA = SWA<Acceptance::PARITY, Level>;

PA::uptr determinize(LevelConfig const& lu);

PA::uptr determinize(LevelConfig const& lu, PS<Acceptance::BUCHI> const& psa, SCCInfo const& psai);

}  // namespace nbautils
