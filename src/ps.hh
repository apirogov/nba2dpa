#pragma once

#include "swa.hh"
#include "common/util.hh"

#include <memory>
#include <queue>
#include <set>
#include <vector>

namespace nbautils {
using namespace std;

using ps_tag = vector<small_state_t>;
template <Acceptance A>
using PS = SWA<A, ps_tag>;
// 2^A for some A
using BAPS = PS<Acceptance::BUCHI>;

using pp_tag = pair<vector<small_state_t>, small_state_t>;
template <Acceptance A>
using PP = SWA<A, pp_tag>;
// 2^AxA for some A
using BAPP = PP<Acceptance::BUCHI>;

// BA -> 2^BA (as reachable from initial state)
template <Acceptance A, typename T>
typename PS<A>::uptr powerset_construction(SWA<A, T> const& ks, vector<small_state_t> const& sinks={}) {
  // create aut, add initial state, associate with initial states in original aut
  state_t const myinit = 0;
  auto pks = std::make_unique<PS<A>>(ks.get_name(), ks.get_aps(), vector<state_t>{myinit});
  pks->tag_to_str = [](auto const& vec){ return "{" + seq_to_str(vec) +"}"; };
  pks->tag->put(to_small_state_t(ks.get_init()), myinit);

  bfs(myinit, [&](auto const& st, auto const& visit, auto const&) {
    // get inner states of current ps state
    auto const curset = pks->tag->geti(st);
    // calculate successors and add to graph
    for (auto i = 0; i < pks->num_syms(); i++) {
      auto const sucset = powersucc(ks, curset, i, sinks);
      if (sucset.empty())
        continue;

      auto const sucst = pks->tag->put_or_get(sucset, pks->num_states());

      if (!pks->has_state(sucst))
        pks->add_state(sucst);
      // add edge
      pks->set_succs(st,i,{sucst});
      // schedule bfs visit of successor
      visit(sucst);
    }
  });

  return move(pks);
}

// BA -> 2^BA x BA, returns basically blown-up original nondet automaton with context annot
template <Acceptance A, typename T>
typename PP<A>::uptr powerset_product(SWA<A, T> const& ks) {
  assert(ks.get_init().size()==1); //2^BA x BA makes only sense with one initial state!

  auto ksinit = ks.get_init().front(); //take unique initial state of original system
  // add initial state, keep acceptance of initial state
  state_t myinit = 0;
  auto pks = std::make_unique<PP<A>>(ks.get_name(), ks.get_aps(), vector<state_t>{myinit});
  pks->tag_to_str = [](auto const& t){ return "({" + seq_to_str(t.first) +"}, "+to_string(t.second)+")"; };

  if (ks.has_accs(ksinit))
    pks->set_accs(myinit, ks.get_accs(ksinit));
  // associate with initial states in original aut
  pks->tag->put(make_pair(to_small_state_t({ksinit}), ksinit), myinit);

  bfs(myinit, [&](auto const& st, auto const& visit, auto const&) {
    // get inner states of current ps state
    auto const curtag = pks->tag->geti(st);
    // remove pointed state from stored set
    auto const curstate = curtag.second;
    auto const& curset = curtag.first;
    // calculate successors
    for (auto i = 0; i < pks->num_syms(); i++) {
      auto const sucset = powersucc(ks, curset, i);  // calc successors of powerset
      //NOTE: maybe here also sucset must be adjusted for accepting sinks?

      auto suctag = make_pair(move(sucset), 0);

      if (!ks.state_has_outsym(curstate,i))
        continue; //dead end

      // for each successor of pointed state add successors
      for (auto const& q : ks.succ(curstate, i)) {
        suctag.second = q;

        // get or create vertex corresponding to this product tag
        // and add an edge to it
        auto const sucst = pks->tag->put_or_get(suctag, pks->num_states());
        if (!pks->has_state(sucst))
          pks->add_state(sucst);
        // add edge
        auto tmp = set_merge(pks->succ(st,i), {sucst});
        pks->set_succs(st,i,tmp);

        // keep acceptance of pointed state
        if (ks.has_accs(q))
          pks->set_accs(sucst, ks.get_accs(q));

        visit(sucst); // schedule for bfs discovery
      }
    }
  });

  return move(pks);
}

}  // namespace nbautils
