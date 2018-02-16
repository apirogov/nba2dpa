#pragma once

#include <bitset>
#include <string>
#include <vector>
#include <ostream>

#include "ps.hh"
#include "relorder.hh"
#include "scc.hh"
#include "types.hh"

namespace nbautils {
using namespace nbautils;

using priority_t = unsigned int;

enum LevelUpdateMode { MUELLERSCHUPP, SAFRA, num };

struct LevelConfig {
  typedef std::unique_ptr<LevelConfig> uptr;
  typedef std::shared_ptr<LevelConfig> sptr;

  bool debug = false;

  // mandatory
  BA const* aut = nullptr;
  SCCInfo const* auti = nullptr;

  // with defaults
  bool sep_rej = false;
  bool sep_acc = false;
  LevelUpdateMode update = MUELLERSCHUPP;

  // optional context information to refine separation
  // if remains nullptr, means that no context is used
  BAPP const* ctx = nullptr;
  SCCInfo const* ctxi = nullptr;

  //TODO: parametrize rank_to_prio mapping? or leave level agnostic?
};

// encodes a set of states as a tuple of disjoint subsets
// augmented with an "importance" ordering on the components.
// Can also be interpreted as tree
struct Level {
  typedef RelOrder::ord_t ord_t;            // type used for relative order of sets
  typedef nbautils::small_state_t state_t;  // state id type of original automaton
  // typedef std::vector<bool> hash_t; // hash to compare for powerset equality
  typedef std::bitset<256> hash_t;  // hash to compare for powerset equality

  std::vector<std::vector<state_t>> tups;  // tuple of indexed disjoint sets
  std::vector<ord_t> tupo;                 // importance order assigned to sets
  // if we keep states from accepting sccs separately, they are last in these vectors

  // sorted, calculated in the end for == with other levels
  hash_t powerset;  // = union of accs tups and nccs
  priority_t prio = 0; //according to min even condition, values from 0 to 2*nbastates+1

  Level();
  Level(LevelConfig const& lvc, std::vector<Level::state_t> const& qs);
  Level succ(LevelConfig const& lvc, sym_t x) const;

  bool operator<(Level const& other) const;
  vector<state_t> states() const;

  string to_string() const;
};

std::ostream& operator<<(std::ostream& os,Level const& l);


}  // namespace nbautils
