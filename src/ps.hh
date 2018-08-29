#pragma once

#include "aut.hh"
#include "common/util.hh"

#include <memory>
#include <queue>
#include <set>
#include <vector>

namespace nbautils {
using namespace std;

using ps_tag = nba_bitset;
// 2^A for some A
using PS = Aut<ps_tag>;

// BA -> 2^BA (as reachable from initial state)
PS powerset_construction(auto const& nba, adj_mat const& mat, nba_bitset const& sinks=0, map<unsigned,nba_bitset> const& impls={}) {
  assert(nba.is_buchi());

  // create aut, add initial state, associate with initial states in original aut
  state_t const myinit = 0;
  auto ps = PS(true, nba.get_name(), nba.get_aps(), myinit);
  ps.tag_to_str = [](ostream& out, ps_tag const& t){ out << pretty_bitset(t); };
  nba_bitset initset = 0;
  initset[nba.get_init()] = 1; // 1<<x does not work as expected
  ps.tag.put(initset, myinit);

  bfs(myinit, [&](auto const& st, auto const& visit, auto const&) {
    // get inner states of current ps state
    auto const curset = ps.tag.geti(st);
    // cerr << "succs of " << pretty_bitset(curset) << endl;

    // calculate successors and add to graph
    for (auto const i : ps.syms()) {
      auto const sucset = powersucc(mat, curset, i, sinks, impls);
      if (sucset == 0)
        continue;

      auto const sucst = ps.tag.put_or_get(sucset, ps.num_states());

      if (!ps.has_state(sucst))
        ps.add_state(sucst);
      // add edge
      ps.add_edge(st,i,sucst);
      // schedule bfs visit of successor
      visit(sucst);
    }
  });

  return ps;
}

using pp_tag = pair<nba_bitset, state_t>;
// 2^AxA for some A
using PP = Aut<pp_tag>;
}

namespace std {
using namespace nbautils;
template <>
    struct hash<pp_tag> {
        size_t operator()(pp_tag const& k) const {
            // Compute individual hash values for first, second and third
            // http://stackoverflow.com/a/1646913/126995
            size_t res = 17;
            res = res * 31 + hash<nba_bitset>()(k.first);
            res = res * 31 + hash<state_t>()(k.second);
            return res;
        }
    };
}

namespace nbautils {

// BA -> 2^BA x BA, returns basically blown-up original nondet automaton with context annot
PP powerset_product(auto const& nba, adj_mat const& mat, nba_bitset const& sinks=0, map<unsigned,nba_bitset> const& impls={}) {
  assert(nba.is_buchi());

  // create aut, add initial state, associate with initial states in original aut
  state_t const myinit = 0;
  auto ps = PP(true, nba.get_name(), nba.get_aps(), myinit);
  ps.tag_to_str = [](ostream& out, pp_tag const& t){
    out << "(" << t.second << ", " << pretty_bitset(t.first) << ")";
  };

  auto const bainit = nba.get_init(); //take unique initial state of original system
  // add initial state, keep acceptance of initial state
  if (nba.state_buchi_accepting(bainit))
    ps.set_pri(myinit, nba.get_pri(bainit));
  // associate with initial states in original aut
  nba_bitset initset = 0;
  initset[bainit] = 1; // 1<<x does not work as expected
  ps.tag.put(make_pair(initset, bainit), myinit);

  bfs(myinit, [&](auto const& st, auto const& visit, auto const&) {
    // get inner states of current ps state
    auto const curtag = ps.tag.geti(st);
    auto const curstate = curtag.second;
    auto const& curset = curtag.first;
    // calculate successors
    for (auto const i : ps.syms()) {
      // calc successors of powerset
      auto const sucset = powersucc(mat, curset, i, sinks, impls);
      auto suctag = make_pair(sucset, 0);

      // for each successor of pointed state add successors
      for (auto const& q : nba.succ(curstate, i)) {
        suctag.second = q;

        // get or create vertex corresponding to this product tag
        // and add an edge to it
        auto const sucst = ps.tag.put_or_get(suctag, ps.num_states());
        if (!ps.has_state(sucst))
          ps.add_state(sucst);
        // add edge
        ps.add_edge(st,i,sucst);

        // keep acceptance of pointed state
        if (nba.state_buchi_accepting(q))
          ps.set_pri(sucst, nba.get_pri(q));

        visit(sucst); // schedule for bfs discovery
      }
    }
  });

  return ps;
}

}  // namespace nbautils
