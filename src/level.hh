#pragma once

#include <vector>
#include <bitset>

#include "types.hh"
#include "relorder.hh"

namespace nbautils {

enum LevelUpdateMode { MUELLERSCHUPP, SAFRA, num };

struct LevelConfig {
	typedef std::shared_ptr<LevelConfig> ptr;

	//mandatory
	BA::ptr aut = nullptr;
	SCCInfo::ptr auti = nullptr;

	//with defaults
	bool sep_rej = false;
	bool sep_acc = false;
	LevelUpdateMode update = MUELLERSCHUPP;

	//optional context information to refine separation
	//if remains nullptr, means that no context is used
	BAPS::ptr ctx = nullptr;
	SCCInfo::ptr ctxi = nullptr;
};

// encodes a set of states as a tuple of disjoint subsets
// augmented with an "importance" ordering on the components.
// Can also be interpreted as tree
struct Level {
  typedef RelOrder::ord_t ord_t; // type used for relative order of sets
  typedef nbautils::small_state_t state_t; // state id type of original automaton
  // typedef std::vector<bool> hash_t; // hash to compare for powerset equality
  typedef std::bitset<256> hash_t; // hash to compare for powerset equality

  std::vector<std::vector<state_t>> tups; // tuple of indexed disjoint sets
  std::vector<ord_t> tupo;     // importance order assigned to sets
  // if we keep states from accepting sccs separately, they are last in these vectors

  // sorted, calculated in the end for == with other levels
  hash_t powerset; // = union of accs tups and nccs
  priority_t prio = 0;
};

Level make_level(LevelConfig const& lvc, std::vector<Level::state_t> const& qs);

Level succ_level(LevelConfig const& lvc, Level l, sym_t x);

} // namespace nbautils

