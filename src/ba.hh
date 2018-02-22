#pragma once
#include <vector>
#include <set>
#include "types.hh"
#include "scc.hh"

namespace nbautils {
  using namespace std;

// this represents a Büchi automaton where
// states are labelled with (the original) state ids
using BA = SWA<Acceptance::BUCHI, state_t>;

vector<state_t> get_accepting_sinks(BA const& aut);

set<scc_t> get_dead_sccs(BA const& aut, SCCInfo const& scci);

size_t trim_ba(BA& ba, SCCInfo& scci, set<scc_t> const& dead = {});

}
