#include "scc.hh"
#include "ba.hh"

namespace nbautils {
using namespace std;

//returns sorted list of accepting sinks (acc. states with self-loop for each sym)
vector<small_state_t> get_accepting_sinks(BA const& aut) {
  vector<small_state_t> ret;
  for (auto const& v : aut.states()) {
    auto outsyms = aut.outsyms(v);
    //must be accepting and have successors for each symbol
    bool accsink = aut.has_accs(v) && outsyms.size() == (size_t)aut.num_syms();
    for (auto i=0; i<aut.num_syms(); i++) {
        // must have self-loop for each symbol
        auto xsucs = aut.succ(v, i);
        if (find(cbegin(xsucs), cend(xsucs), v) == cend(xsucs))
          accsink = false;
    }
    if (accsink)
      ret.push_back(v);
  }
  return ret;
}

void mark_dead_sccs(BA const& aut, SCCInfo const& scci, map<scc_t,bool>& dead, scc_t num) {
  if (map_has_key(dead, num))  // done already
    return;

  // std::cout << "scc " << num << std::endl;
  auto sucsccs = succ_sccs(aut, scci, num);
  // mark children first
  for (auto sucscc : sucsccs) mark_dead_sccs(aut, scci, dead, sucscc);

  // if we are rejecting and trivial, assume we're dead
  bool isdead = contains(scci.rejecting, num) || contains(scci.trivial, num);
  // check children and try to falsify
  for (auto sucscc : sucsccs) isdead = isdead && dead[sucscc];
  // if still dead, we're really dead.
  dead[num] = isdead;
}

// run dfs that marks dead sccs
set<scc_t> get_dead_sccs(BA const& aut, SCCInfo const& scci) {
  map<scc_t, bool> dead;
  for (auto &v : aut.get_init())
    mark_dead_sccs(aut, scci, dead, scci.scc.at(v));

  set<scc_t> ret;
  for (auto const& it : dead)
    if (it.second)
      ret.emplace(it.first);
  return ret;
}


// use SCC info to purge unreachable and dead states from automaton and scc info
// (but keeps initial states, even when they are dead)
size_t trim_ba(BA& ba, SCCInfo& scci, set<scc_t> const& dead) {
  set<state_t> erase;

  // collect unreachable
  swap(scci.unreachable, erase);

  if (erase.empty() && dead.empty())
    return 0; //nothing to do

  // collect dead states
  for (auto const& s : ba.states()) {
    if (map_has_key(scci.scc, s)) {  // has assigned scc
      auto scit = scci.scc.at(s);
      if (contains(dead,scit) && !contains(ba.get_init(), s)) {  // is a dead scc
        erase.emplace(s);
        scci.scc.erase(s);
      }
    }
  }
  // unmark these states
  for (auto& scc : dead) {
    scci.accepting.erase(scc);
    scci.rejecting.erase(scc);
    scci.trivial.erase(scc);
    scci.sccrep.erase(scc);
    scci.sccsz.erase(scc);
  }

  vector<state_t> erasevec(begin(erase), end(erase));
  ba.remove_states(erasevec);
  return erasevec.size();
}

}
