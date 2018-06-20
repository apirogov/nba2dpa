#pragma once

#include <bitset>
#include <string>
#include <vector>
#include <ostream>

#include "common/scc.hh"
#include "preproc.hh"
#include "common/relorder.hh"
#include "ps.hh"

namespace nbautils {
using namespace nbautils;

enum class UpdateMode { MUELLERSCHUPP, SAFRA, FULLMERGE, num };

//no opt:
//nscc_states empty
//ascc_states empty
//mscc_states contains everything in one

//sep nscc:
//nscc_states contains nscc states
//mscc_states contains rest in one

//sep (nscc+)ascc:
//(nscc_states filled)
//ascc_states filled
//mscc_states contains rest in one

//sep (nscc+ascc+)mscc:
//(nscc_states filled
//ascc_states filled)
//mscc_states is a vector with different sccs

//optimizing detsccs:
//as above, but have det_sccs among msccs marked extra

struct DetConf {
  bool debug = false;

  // mandatory SCC infos that are used with various heuristics
  adj_mat aut_mat;            //adj matrix
  nba_bitset aut_states = 0; //all used states
  nba_bitset aut_acc = 0;    //accepting states

  //invariants: nscc_states + ascc_states + dscc_states + mscc_states = aut_states
  //            nscc_states , ascc_states , dscc_states , mscc_states all pw disj
  //            union of asccs_states = ascc_states
  //            asccs_states pw disj
  //            msccs_states pw disj

  nba_bitset nscc_states = 0; //nsccs, merged together

  nba_bitset ascc_states = 0; //asccs, merged together
  vector<nba_bitset> asccs_states; //either each ascc sep. or contains exactly ascc_states

  nba_bitset dscc_states = 0; //msccs that are deterministic, merged together

  vector<nba_bitset> msccs_states; //either each scc sep. or all remaining mscc together

  //heuristics and optimisations:
  nba_bitset aut_asinks = 0;  //if non-empty, will be used to stop early
  Context ctx; //if non-empty, context used for seperation refinement

  UpdateMode update = UpdateMode::MUELLERSCHUPP; // kind of merge
  bool weaksat = false;   //already saturate if some single SCC moved down
  bool puretrees = false; //move accepting states always into leaves after merges

  bool sep_rej = false;
  bool sep_acc = false;
  bool sep_acc_cyc = false;
  bool sep_mix = false;
  bool opt_det = false;
};

std::ostream& operator<<(std::ostream& os, DetConf const& dc);

// encodes a set of states as a tuple of disjoint subsets
// augmented with an "importance" ordering on the components.
// Can also be interpreted as tree / forest
struct DetState {
  nba_bitset powerset=0; // = all states present in this state (useful for lang. eq. comparison)

  nba_bitset nsccs=0; //storage for (relatively) nonacc. SCC states

  //storage for (relatively) acc. SCC state breakpoint construction + assigned rank
  nba_bitset asccs_buf=0;
  nba_bitset asccs=0;
  pri_t asccs_pri=0;

  //deterministic mixed SCC(s) - nodes not expanded, other saturation condition
  vector<pair<nba_bitset,pri_t>> dsccs;

  //(remaining) mixed SCC(s) - MS tuple / Safra tree(s)
  vector<vector<pair<nba_bitset,pri_t>>> msccs;

  DetState();
  DetState(DetConf const& dc, nba_bitset const& qs);

  //given a state and symbol returns successor and edge priority
  pair<DetState, pri_t> succ(DetConf const& dc, sym_t x) const;

  bool operator==(DetState const& other) const;
  // bool operator<(DetState const& other) const;

  string to_string() const;
};

std::ostream& operator<<(std::ostream& os, DetState const& s);

}  // namespace nbautils
