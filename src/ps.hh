#pragma once

#include "bench.hh"
#include "types.hh"
#include "util.hh"

#include <memory>
#include <queue>
#include <set>
#include <vector>

namespace nbautils {
using namespace std;

using ps_tag = vector<small_state_t>;
// using ps_tag_store = generic_trie_bimap<ps_tag, small_state_t, state_t>;
using ps_tag_store = naive_bimap<ps_tag, state_t>;
template <Acceptance A>
using PS = SWA<A, ps_tag, ps_tag_store>;
// 2^A for some A
using BAPS = PS<Acceptance::BUCHI>;

using pp_tag = pair<vector<small_state_t>, small_state_t>;
// using pp_tag_store = generic_trie_bimap<pp_tag, small_state_t, state_t>;
using pp_tag_store = naive_bimap<pp_tag, state_t>;
template <Acceptance A>
using PP = SWA<A, pp_tag, pp_tag_store>;
// 2^AxA for some A
using BAPP = PP<Acceptance::BUCHI>;

// BA -> 2^BA
template <Acceptance A, typename T>
typename PS<A>::uptr powerset_construction(SWA<A, T> const& ks) {
  // auto const idmap = [](auto const& v) { return v; };
  // auto tag = std::make_unique<ps_tag_store>(ps_tag_store(idmap, idmap));
  auto tag = std::make_unique<ps_tag_store>(ps_tag_store());
  auto pksp = std::make_unique<PS<A>>(PS<A>());

  auto& pks = *pksp;
  pks.tag = move(tag);
  auto& pbm = pks.tag;

  pks.add_state(ks.init);
  pks.init = ks.init;
  pks.meta = ks.meta;

  // q_0 -> {q_0}
  pbm->put(std::vector<small_state_t>{(small_state_t)ks.init}, pks.init);

  bfs(pks.init, [&](auto const& st, auto const& visit) {
    auto curset = pbm->get(st);  // get inner states of current ps state
    // calculate successors
    for (auto i = 0; i < (int)pks.num_syms(); i++) {
      auto const sucset = powersucc(ks, curset, i);  // calc successors
      state_t sucst = pbm->put_or_get(sucset, pbm->size());
      pks.adj[st][i].push_back(sucst);
      visit(sucst);
    }
  });

  return move(pksp);
}

// BA -> 2^BA x BA, returns basically blown-up original automaton with context annot
template <Acceptance A, typename T>
typename PP<A>::uptr powerset_product(SWA<A, T> const& ks) {
  auto pksp = std::make_unique<PP<A>>(PP<A>());
  auto& pks = *pksp;

  /*
  auto from = [](auto const& v) {
    vector<small_state_t> v2(v.first.begin(), v.first.end());
    v2.push_back(v.second);
    return v2;
  };
  auto to = [](auto const& v) {
    vector<small_state_t> v2(v.begin(), v.end());
    auto p = v2.back();
    v2.pop_back();
    return make_pair(v, p);
  };
  pks.tag = std::make_unique<pp_tag_store>(pp_tag_store(from, to));
  */
  pks.tag = std::make_unique<pp_tag_store>(pp_tag_store());
  auto& pbm = pks.tag;

  // same letters etc
  pks.meta= ks.meta;
  // initial state is 0
  pks.add_state(0);
  pks.init = 0;
  // keep acceptance of initial state
  if (ks.has_accs(ks.init)) pks.set_accs(pks.init, ks.get_accs(ks.init));

  // q_0 -> ({q_0}, q_0)
  auto pinit = make_pair(std::vector<small_state_t>{(small_state_t)ks.init}, ks.init);
  pbm->put(pinit, pks.init);

  bfs(pks.init, [&](auto const& st, auto const& visit) {
    // get inner states of current ps state
    auto curtag = pbm->get(st);
    // remove pointed state from stored set
    auto curstate = curtag.second;
    auto& curset = curtag.first;
    // calculate successors
    for (auto i = 0; i < (int)pks.num_syms(); i++) {
      auto const sucset = powersucc(ks, curset, i);  // calc successors of powerset
      auto suctag = make_pair(move(sucset), 0);

      if (!ks.state_has_outsym(curstate,i))
        continue; //deadend

      // for each successor of pointed state add successors
      for (auto const& q : ks.succ(curstate, i)) {
        suctag.second = q;

        // get or create vertex corresponding to this product tag
        // and add an edge to it
        auto sucst = pbm->put_or_get(suctag, pbm->size());
        if (!pks.has_state(sucst))
          pks.add_state(sucst);
        pks.adj[st][i].push_back(sucst);

        // keep acceptance of pointed state
        if (ks.has_accs(q)) pks.set_accs(sucst, ks.get_accs(q));

        visit(sucst); // add for discovery
      }
    }
  });

  return move(pksp);
}

}  // namespace nbautils
