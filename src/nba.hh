#pragma once

#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include <spot/twa/twagraph.hh>
#include <spot/twaalgos/sccinfo.hh>

namespace nbautils {
spot::twa_graph_ptr spot_nba_from_file(std::string const &file);

// we use a quite primitive representation of an NBA
// TODO: keep ap string to sym bit mapping preserved
class NBA {
 public:
  typedef unsigned state_t;
  typedef unsigned sym_t;

  std::shared_ptr<spot::scc_info> scci;  // scc information

  std::vector<spot::formula>
      aps;  // ordered atomic propositions in lsb order, unused come last

  size_t num_syms;           // number of actual symbols 0..(2^AP - 1)
  state_t init;              // initial state
  std::set<state_t> states;  // set of all states
  std::set<state_t> acc;     // subset of accepting states

  // list of successors per symbol
  std::map<state_t, std::map<sym_t, std::vector<state_t>>> adj;

  // collect x-successors from all states in qs
  std::vector<state_t> powersucc(std::vector<state_t> const &qs, sym_t const &x);

  // currently we can only get our NBA from spot
  NBA(spot::twa_graph_ptr aut);
};

}  // namespace nbautils
