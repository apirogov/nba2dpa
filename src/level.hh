#pragma once

#include <bitset>
#include <string>
#include <vector>
#include <ostream>

#include "common/scc.hh"
#include "common/algo.hh"
#include "common/relorder.hh"
#include "ps.hh"

namespace nbautils {
using namespace nbautils;


enum class LevelUpdateMode { MUELLERSCHUPP, SAFRA, num };

struct LevelConfig {
  typedef std::unique_ptr<LevelConfig> uptr;
  typedef std::shared_ptr<LevelConfig> sptr;

  bool debug = false;

  // mandatory
  SWA<string> const* aut = nullptr;
  SCCDat<state_t> const* aut_scc = nullptr;
  BaSccClassification const* aut_cl = nullptr;

  // with defaults
  vector<small_state_t> accsinks; //if non-empty, will be used
  bool sep_rej = false;
  bool sep_acc = false;
  bool sep_acc_cyc = false;
  LevelUpdateMode update = LevelUpdateMode::MUELLERSCHUPP;

  bool pure = false; //move accepting states always into leaves after merges
  bool optguarded = false; //reuse trees if priority is guarded

  // optional context information to refine separation
  // if remains nullptr, means that no context is used
  PP const* ctx = nullptr;
  SCCDat<state_t> const* ctx_scc = nullptr;
  BaSccClassification const* ctx_cl = nullptr;
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

  bool same_tree(Level const& other) const;
  bool operator==(Level const& other) const;
  bool operator<(Level const& other) const;
  vector<state_t> states() const;

  string to_string() const;
};

std::ostream& operator<<(std::ostream& os,Level const& l);


}  // namespace nbautils
