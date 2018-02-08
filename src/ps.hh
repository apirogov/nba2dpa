#pragma once

#include "types.hh"
#include "debug.hh"
#include <memory>
#include <set>
#include <vector>
#include <queue>

#include <spdlog/spdlog.h>
namespace spd = spdlog;

namespace nbautils {

// BA -> 2^BA
template<typename L,typename T>
typename PS<L>::ptr powerset_construction(SWA<L,T> const& ks) {
  spd::get("log")->debug("powerset_construction(aut)");
  auto starttime = get_time();

  auto pksp = std::make_shared<PS<L>>(PS<L>());
  auto &pks = *pksp;

  pks.tag = std::make_shared<ps_tag>(ps_tag());
  auto& pbm = pks.tag;

  // same letters
  pks.num_syms = ks.num_syms;
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
    if (pks.adj.find(st) != end(pks.adj))
      continue; // have visited this one

    // calculate successors
    for (auto i = 0; i < pks.num_syms; i++) {
      auto curset = pbm->get(st); // get inner states of current ps state
      auto sucset = powersucc(ks, curset, i); // calc successors
      auto sucst = pbm->put_or_get(sucset, pbm->size());
      pks.adj[st][i].push_back(sucst);
      bfsq.push(sucst);
    }
  }

  spd::get("log")->debug("powerset_construction completed ({:.4f} s)"
		  , get_secs_since(starttime));
  return pksp;
}

// BA -> BA x 2^BA, returns basically blown-up original automaton with context annot
template<typename L, typename T>
typename PS<L>::ptr powerset_product(SWA<L,T> const& ks) {
  spd::get("log")->debug("powerset_product(aut)");
  auto starttime = get_time();

  auto pksp = std::make_shared<PS<L>>(PS<L>());
  auto &pks = *pksp;

  pks.tag = std::make_shared<ps_tag>(ps_tag());
  auto& pbm = pks.tag;

  // same letters
  pks.num_syms = ks.num_syms;
  // initial state is 0
  pks.init = 0;
  //keep acceptance of initial state
  if (ks.acc.find(ks.init) != end(ks.acc))
    pks.acc[pks.init] = ks.acc.at(ks.init);

  // q_0 -> ({q_0}, q_0)
  pbm->put(std::vector<small_state_t>{(small_state_t)ks.init,(small_state_t)ks.init}, pks.init);

  // explore powerset by bfs
  std::queue<state_t> bfsq;
  bfsq.push(pks.init);
  while (!bfsq.empty()) {
    auto const st = bfsq.front();
    bfsq.pop();
    if (pks.adj.find(st) != end(pks.adj))
      continue; // have visited this one

    // calculate successors
    for (auto i = 0; i < pks.num_syms; i++) {
      auto curset = pbm->get(st); // get inner states of current ps state
      auto curstate = curset.back();
      curset.pop_back(); //remove pointed state from stored set
      auto sucset = powersucc(ks, curset, i); // calc successors of powerset

      // for each successor of pointed state add successors
      for (auto q : succ(ks, curstate, i)) {
        sucset.push_back(q); //add new pointed state

        auto sucst = pbm->put_or_get(sucset, pbm->size());
        pks.adj[st][i].push_back(sucst);

        //keep acceptance of pointed state
        if (ks.acc.find(q) != end(ks.acc))
          pks.acc[sucst] = ks.acc.at(q);

        bfsq.push(sucst);

        sucset.pop_back(); //remove pointed state
      }
    }
  }


  spd::get("log")->debug("powerset_product completed ({:.4f} s)"
		  , get_secs_since(starttime));
  return pksp;
}

}
