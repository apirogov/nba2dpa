#pragma once

#include <bitset>
#include <string>
#include <vector>
#include <ostream>

#include "aut.hh"
#include "common/scc.hh"
#include "preproc.hh"
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

struct DetConfSets {
  //invariants: nscc_states + ascc_states + dscc_states + mscc_states = aut_states
  //            nscc_states , ascc_states , dscc_states , mscc_states all pw disj
  //            union of asccs_states = ascc_states
  //            asccs_states pw disj
  //            msccs_states pw disj

  nba_bitset nscc_states = 0; //nsccs, merged together

  nba_bitset ascc_states = 0; //asccs, merged together
  vector<nba_bitset> asccs_states; //either each ascc sep. or contains exactly ascc_states

  vector<nba_bitset> dsccs_states; //deterministic msccs, separately

  vector<nba_bitset> msccs_states; //either each scc sep. or all remaining mscc together
};

struct DetConf {
  bool debug = false;

  // mandatory SCC infos that are used with various heuristics
  adj_mat aut_mat;           //adj matrix
  nba_bitset aut_states = 0; //all used states
  nba_bitset aut_acc = 0;    //accepting states

  //heuristics and optimisations:
  nba_bitset aut_asinks = 0;  //if non-empty, will be used to stop early
  map<unsigned,nba_bitset> impl_mask; //to store implication relation
  map<unsigned,nba_bitset> impl_pruning_mask; //to store implication relation disregarding SCC relationship
  Context ctx;                //if non-empty, context used for seperation refinement
  int maxsets = 1;

  //these must be filled
  DetConfSets sets;

  UpdateMode update = UpdateMode::MUELLERSCHUPP; // kind of merge
  bool puretrees = false; //move accepting states always into leaves after merges

  bool sep_rej = false;
  bool sep_acc = false;
  bool sep_acc_cyc = false;
  bool sep_mix = false;
  bool opt_det = false;
  bool opt_suc = false;
};

DetConfSets calc_detconfsets(DetConf const& dc, SCCDat const& scci,
    BASccAClass const& sccAcc, set<unsigned> const& sccDet);

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
  vector<ranked_slice> dsccs;

  //(remaining) mixed SCC(s) - MS tuple / Safra tree(s)
  vector<ranked_slice> msccs;

  DetState();
  DetState(DetConf const& dc, nba_bitset const& qs);

  //given a state and symbol returns successor and edge priority
  pair<DetState, pri_t> succ(DetConf const& dc, sym_t x) const;

  bool operator==(DetState const& other) const;
  bool operator!=(DetState const& other) const;
  // bool operator<(DetState const& other) const;

  string to_string() const;

  tree_history to_tree_history() const;
  bool tuples_finer_or_equal(DetState const&) const;
};

std::ostream& operator<<(std::ostream& os, DetState const& s);

}  // namespace nbautils

namespace std {
using namespace nbautils;
template <>
    struct hash<DetState> {
        size_t operator()(DetState const& k) const {
            // Compute individual hash values for first, second and third
            // http://stackoverflow.com/a/1646913/126995
            size_t res = 17;
            res = res * 31 + hash<nba_bitset>()(k.powerset);
            res = res * 31 + hash<nba_bitset>()(k.nsccs);
            res = res * 31 + hash<nba_bitset>()(k.asccs_buf);
            res = res * 31 + hash<nba_bitset>()(k.asccs);
            res = res * 31 + hash<pri_t>()(k.asccs_pri);
            for (auto const& dscc : k.dsccs)
              for (auto const& it : dscc) {
                res = res * 31 + hash<nba_bitset>()(it.first);
                res = res * 31 + hash<pri_t>()(it.second);
              }
            for (auto const& mscc : k.msccs)
              for (auto const& it : mscc) {
                res = res * 31 + hash<nba_bitset>()(it.first);
                res = res * 31 + hash<pri_t>()(it.second);
              }
            return res;
        }
    };
}
