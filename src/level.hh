#pragma once

#include <vector>

#include "nba.hh"
#include "relorder.hh"

namespace nbautils {

// encodes a set of states as a tuple of disjoint subsets
// augmented with an "importance" ordering on the components.
// Can also be interpreted as tree
class StateSieve {
  typedef RelOrder::ord_t ord_t;           // type used for relative order of sets
  typedef nbautils::NBA::state_t state_t;  // state id type of original automaton
  typedef std::vector<state_t>
      state_set;  // how a state set is represented. a sorted vector

  state_set accs;  // we keep states from accepting sccs separately

  std::vector<state_set> tups;  // tuple of indexed disjoint sets
  std::vector<ord_t> tupo;      // importance order assigned to sets

  state_set nccs;  // we keep states from non-acc. sccs separately

  // accs, nccs and all tups disjoint

  // sorted, calculated in the end for == with other levels
  state_set const powerset;  // = union of accs tups and nccs

  // priority of this sieve,
  // assigned once during creation depending on predecessor
  int priority = 0;

  // create new level from unstructured state set
  StateSieve(state_set qs);
};

}  // namespace nbautils
