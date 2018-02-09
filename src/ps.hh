#pragma once

#include "bench.hh"
#include "types.hh"

#include <memory>
#include <queue>
#include <set>
#include <vector>

#include <spdlog/spdlog.h>
namespace spd = spdlog;

namespace nbautils {
using namespace std;

using ps_tag = vector<small_state_t>;
using ps_tag_store = generic_trie_bimap<ps_tag, small_state_t, state_t>;
template <typename L>
using PS = SWA<L, ps_tag, ps_tag_store>;
// 2^A for some A
using BAPS = PS<bool>;

using pp_tag = pair<vector<small_state_t>, small_state_t>;
// using pp_tag_store = generic_trie_bimap<pp_tag, small_state_t, state_t>;
using pp_tag_store = naive_bimap<pp_tag, state_t>;
template <typename L>
using PP = SWA<L, pp_tag, pp_tag_store>;
// 2^AxA for some A
using BAPP = PP<bool>;

// BA -> 2^BA
template <typename L, typename T>
typename PS<L>::uptr powerset_construction(SWA<L, T> const& ks) {
  spd::get("log")->debug("powerset_construction(aut)");
  auto starttime = get_time();

  auto pksp = std::make_unique<PS<L>>(PS<L>());
  auto& pks = *pksp;

  auto idmap = [](auto const& v) { return v; };
  pks.tag = std::make_unique<ps_tag_store>(ps_tag_store(idmap, idmap));
  auto& pbm = pks.tag;

  // same letters
  pks.meta = ks.meta;
  // initial state is 0
  pks.init = 0;
  // q_0 -> {q_0}
  pbm->put(std::vector<small_state_t>{(small_state_t)ks.init}, pks.init);

  // explore powerset by bfs
  std::queue<state_t> bfsq;
  bfsq.push(pks.init);
  while (!bfsq.empty()) {
    auto const st = bfsq.front();
    bfsq.pop();
    if (pks.has_state(st)) continue;  // have visited this one

    auto curset = pbm->get(st);  // get inner states of current ps state
    // calculate successors
    for (auto i = 0; i < pks.meta.num_syms; i++) {
      auto sucset = powersucc(ks, curset, i);  // calc successors
      state_t sucst = pbm->put_or_get(sucset, pbm->size());
      pks.adj[st][i].push_back(sucst);
      bfsq.push(sucst);
    }
  }

  spd::get("log")->debug("powerset_construction completed ({:.4f} s)",
                         get_secs_since(starttime));
  return move(pksp);
}

// BA -> 2^BA x BA, returns basically blown-up original automaton with context annot
template <typename L, typename T>
typename PP<L>::uptr powerset_product(SWA<L, T> const& ks) {
  spd::get("log")->debug("powerset_product(aut)");
  auto starttime = get_time();

  auto pksp = std::make_unique<PP<L>>(PP<L>());
  auto& pks = *pksp;

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
  // pks.tag = std::make_unique<pp_tag_store>(pp_tag_store(from, to));
  pks.tag = std::make_unique<pp_tag_store>(pp_tag_store());
  auto& pbm = pks.tag;

  // same letters
  pks.meta.num_syms = ks.meta.num_syms;
  // initial state is 0
  pks.init = 0;
  // keep acceptance of initial state
  if (ks.has_acc(ks.init)) pks.acc[pks.init] = ks.acc.at(ks.init);

  // q_0 -> ({q_0}, q_0)
  auto pinit = make_pair(std::vector<small_state_t>{(small_state_t)ks.init}, ks.init);
  pbm->put(pinit, pks.init);

  // explore powerset by bfs
  std::queue<state_t> bfsq;
  bfsq.push(pks.init);
  while (!bfsq.empty()) {
    auto const st = bfsq.front();
    bfsq.pop();
    if (pks.has_state(st)) continue;  // have visited this one

    // get inner states of current ps state
    auto curtag = pbm->get(st);
    // remove pointed state from stored set
    auto curstate = curtag.second;
    auto& curset = curtag.first;
    // calculate successors
    for (auto i = 0; i < pks.meta.num_syms; i++) {
      auto sucset = powersucc(ks, curset, i);  // calc successors of powerset
      auto suctag = make_pair(move(sucset), 0);

      // for each successor of pointed state add successors
      for (auto const& q : ks.succ(curstate, i)) {
        // sucset.push_back(q); //add new pointed state
        suctag.second = q;

        // get or create vertex corresponding to this product tag
        // and add an edge to it
        auto sucst = pbm->put_or_get(suctag, pbm->size());
        pks.adj[st][i].push_back(sucst);

        // keep acceptance of pointed state
        if (ks.has_acc(q)) pks.acc[sucst] = ks.acc.at(q);

        // add for discovery
        bfsq.push(sucst);

        // sucset.pop_back(); //remove pointed state
      }
    }
  }

  spd::get("log")->debug("powerset_product completed ({:.4f} s)",
                         get_secs_since(starttime));
  return move(pksp);
}

}  // namespace nbautils
