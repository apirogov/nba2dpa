#include "ba.hh"

#include "swa.hh"
#include "common/scc.hh"

namespace nbautils {
using namespace std;

// find and purge unreachable and dead states from automaton
// (but keeps initial states, even when they are dead)
size_t trim_ba(SWA<string>& ba) {
  auto const ba_st = ba.states();
  auto const ba_suc = swa_succ(ba);
  auto const ba_acc = swa_ba_acc(ba);

  auto const scci = get_sccs(ba_st, ba_suc);
  auto const scc_suc = [&](auto num) { return succ_sccs(*scci, num, ba_suc); };

  auto const unreach = unreachable_states(ba_st, ba.get_init().front(), ba_suc);

  auto const triv = trivial_sccs(*scci, ba_suc);
  auto const cl = ba_classify_sccs(*scci, ba_acc);
  auto const deadscc = ba_get_dead_sccs(scci->sccs.size(), cl->rejecting, triv, scc_suc);

  // collect dead states
  vector<state_t> dead;
  for (auto const& s : ba_st) {
      auto scit = scci->scc_of.at(s);
      if (contains(deadscc,scit) && !contains(ba.get_init(), s)) {  // is a dead scc
        dead.push_back(s);
      }
  }

  auto const erase = set_merge(unreach, dead);
  ba.remove_states(erase);
  return erase.size();
}

}
