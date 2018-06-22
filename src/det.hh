#pragma once

#include <queue>
#include <set>
#include <vector>
#include "detstate.hh"
#include "common/scc.hh"
#include "aut.hh"

namespace nbautils {

using PA = Aut<DetState>;

// BFS-based determinization with supplied level update config
PA determinize(auto const& nba, DetConf const& dc, nba_bitset const& startset, auto const& pred) {
  assert(nba.is_buchi());
  // create automaton with same letters etc
  state_t const myinit = 0;
  auto pa = PA(false, nba.get_name(), nba.get_aps(), myinit);
  pa.set_patype(PAType::MIN_EVEN);
  pa.tag_to_str = [](ostream& out, auto const& ds){ out << ds; };
  pa.tag.put(DetState(dc, startset), myinit); // initial state tag

  int numvis=0;
  bfs(myinit, [&](auto const& st, auto const& visit, auto const&) {
    // get inner states of current ps state
    auto const cur = pa.tag.geti(st);

    // cout << "visit " << curlevel.to_string() << endl;
    ++numvis;
    if (numvis % 5000 == 0) //progress indicator
      cerr << numvis << endl;

    for (auto const i : pa.syms()) {
      // calculate successor level
      DetState suclevel;
      pri_t sucpri;
      tie(suclevel, sucpri) = cur.succ(dc, i);
      // cout << "suc " << suclevel.to_string() << endl;

      if (suclevel.powerset == 0) //is an empty set -> invalid successor
        continue;

      if (!pred(suclevel)) //predicate not satisfied -> don't explore this node
        continue;

      //check whether there is already a state in the graph with this label
      auto const sucst = pa.tag.put_or_get(suclevel, pa.num_states());

      //if this is a new successor, add it to graph and enqueue it:
      if (!pa.has_state(sucst))
        pa.add_state(sucst);
      // create edge
      pa.add_edge(st, i, sucst, sucpri);
      // schedule for bfs
      visit(sucst);
    }
  });

  // cerr << "determinized to " << pa->num_states() << " states" << endl;
  return pa;
}

// start with initial state of NBA, explore by DFS completely
PA determinize(auto const& nba, DetConf const& dc) {
  return determinize(nba, dc, nba_bitset(1<<nba.get_init()), const_true);
}

}  // namespace nbautils
