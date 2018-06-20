#pragma once

#include <memory>
#include <queue>
#include <set>
#include <vector>
#include "detstate.hh"
#include "common/scc.hh"
#include "aut.hh"

namespace nbautils {

using PA = Aut<DetState>;

PA determinize(DetConf const& lu);

PA determinize(DetConf const& lu, PS const& psa, SCCDat const& psai);

}  // namespace nbautils
