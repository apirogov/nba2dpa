#pragma once

#include "types.hh"
#include "debug.hh"
#include <stack>
#include <set>
#include <utility>
#include <iostream>
#include <algorithm>
#include <vector>
#include <queue>

#include <spdlog/spdlog.h>
namespace spd = spdlog;

namespace nbautils {

/*
template<typename L,typename T>
std::vector<state_t> scc_states(SWA<L,T> const& aut, SCCInfo const& scci, scc_t const& num) {
  std::queue<state_t> bfsq;
  std::set<state_t> visited;
  bfsq.push(scci.sccrep.at(num));
  while (!bfsq.empty()) {
    auto const st = bfsq.front();
    bfsq.pop();
    if (visited.find(st) != end(visited))
      continue; // have visited this one
    for (auto i = 0; i < aut.num_syms; i++) {
      auto const sucs = succ(aut,st,i);
      for (auto const& sucst : sucs) {
        if (scci.scc.at(sucst)==num) {
          bfsq.push(sucst);
        }
      }
    }
  }
  std::vector<nbautils::scc_t> ret;
  copy(begin(visited),end(visited),back_inserter(ret));
  return ret;
}
*/

template<typename L,typename T>
std::vector<state_t> succ_sccs(SWA<L,T> const& aut, SCCInfo const& scci, scc_t const& num) {
  std::queue<state_t> bfsq;
  std::set<state_t> visited;
  std::set<nbautils::scc_t> sucsccs;

  bfsq.push(scci.sccrep.at(num));
  while (!bfsq.empty()) {
    auto const st = bfsq.front();
    bfsq.pop();
    if (visited.find(st) != end(visited))
      continue; // have visited this one
    visited.emplace(st);

    for (auto i = 0; i < aut.num_syms; i++) {
      auto sucs = succ(aut,st,i);
      for (auto& sucst : sucs) {
        auto sucscc = scci.scc.at(sucst);
      if (sucscc==num)
        bfsq.push(sucst);
      else
        sucsccs.emplace(sucscc);
      }
    }
  }
  std::vector<nbautils::scc_t> ret;
  copy(begin(sucsccs),end(sucsccs),back_inserter(ret));
  return ret;
}

template<typename L,typename T>
void mark_dead_sccs(SWA<L,T> const& aut, SCCInfo& scci, scc_t num) {
  if (scci.dead.find(num) != end(scci.dead)) //done already
    return;

  // std::cout << "scc " << num << std::endl;
  auto sucsccs = succ_sccs(aut, scci, num);
  //mark children first
  for (auto sucscc : sucsccs)
    mark_dead_sccs(aut, scci, sucscc);

  //if we are rejecting and trivial, assume we're dead
  bool dead = (scci.rejecting.find(num) != end(scci.rejecting)
                  || scci.trivial.find(num) != end(scci.trivial));
  //check children and try to falsify
  for (auto sucscc : sucsccs)
    dead = dead && scci.dead[sucscc];
  //if still dead, we're really dead.

  scci.dead[num] = dead;
}


//https://en.wikipedia.org/wiki/Path-based_strong_component_algorithm with extensions
template<typename L,typename T>
SCCInfo::ptr get_scc_info(SWA<L,T> const& aut, bool analyse_acc = true) {
  spd::get("log")->debug("get_scc_info(aut,{})",analyse_acc);
  auto starttime = get_time();

  auto sccip = std::make_shared<SCCInfo>(SCCInfo());
  auto &scci = *sccip;

  std::stack<state_t> call; // dfs call stack
  std::stack<state_t> reps; // scc representative stack
  std::stack<state_t> open; // not yet fully completed vertex stack

  std::map<state_t, unsigned> order; // first visit order
  int count = 0;

  // we just start with the initial state
  call.push(aut.init);
  while (!call.empty()) {
    auto v = call.top();
    if (order.find(v) == order.end()) { // dfs just "called" with current node
      order[v] = count++; //assign visit order

      reps.push(v); //SCC representative candidate (for now)
      open.push(v); //this node is "pending" (not completely discovered from here)

      for (auto i = 0; i < aut.num_syms; i++) {
        auto const& vesit = aut.adj.find(v);
        if (vesit == end(aut.adj))
          continue;
        auto const& ves = vesit->second;
        auto const& vsucsit = ves.find(i);
        if (vsucsit == end(ves))
          continue;
        auto const& vsucs = vsucsit->second;

        for (auto w : vsucs) { // process edges with any label

          if (order.find(w) == order.end()) {
            call.push(w); //recursively explore nodes that have not been visited yet
          } else if (scci.scc.find(w) == scci.scc.end()) {
            // if already visited, but not with assigned scc, we have found a loop
            // -> drop candidates, keep oldest on this loop as SCC representative
            while (order[reps.top()] > order[w])
              reps.pop();
          }

        }
      }

    } else { // returned from recursive calls
      // is still rep. -> we found an SCC
      if (reps.top() == v) {
        // assume SCC is accepting and rejecting, then falsify
        bool accscc = true;
        bool rejscc = true;
        //count how many state this scc has to detect triv. SCC
        int scc_sz = 0;
        auto scc_num = scci.sccrep.size(); //size corresponds to next number

        // drop states up to the current state, they are done and part of the SCC
        state_t tmp;
        do {
          tmp = open.top();
          open.pop();

          //if acceptance info given, detect acc./rej. sccs
          if (analyse_acc) {
            auto actmp = aut.acc.find(tmp) != end(aut.acc);
            accscc = accscc && actmp;
            rejscc = rejscc && !actmp;
          }

          //assign current state current SCC
          scci.scc[tmp] = scc_num;

          scc_sz++;
        } while (tmp != v);

        // now lets get the extra information we want
        if (analyse_acc) {

          // mark completed SCC as accepting/rejecting
          if (accscc)
            scci.accepting.emplace(scc_num);
          if (rejscc)
            scci.rejecting.emplace(scc_num);

          // trivial scc with no self-loop?
          bool trivacc = scc_sz == 1;
          bool noselfloop = true;
          for (auto i = 0; i < aut.num_syms; i++) {
            auto const vsucs = succ(aut, v, i);
            if (find(begin(vsucs),end(vsucs), tmp)!=end(vsucs))
              noselfloop = false;
          }
          if (trivacc && noselfloop)
            scci.trivial.emplace(scc_num);
        }

        // store the representative state (from which the rest can be obtained easily)
        scci.sccrep[scc_num] = reps.top();

        // current SCC is done
        reps.pop();
      }

      call.pop(); // current node done -> dfs "returns" to previous caller
    }
  }
  // states still without assigned scc are unreachable
  for (auto &it : aut.adj)
    if (scci.scc.find(it.first) == end(scci.scc))
      scci.unreachable.emplace(it.first);

  // run dfs that marks dead sccs
  mark_dead_sccs(aut, scci, scci.scc.at(aut.init));

  spd::get("log")->debug("get_scc_info completed ({:.4f} s)"
		  , get_secs_since(starttime));

  return sccip;
}

// use SCC info to purge unreachable and dead states from automaton and scc info
template <typename L, typename T>
size_t trim_ba(SWA<L,T>& ba, SCCInfo& scci) {
  std::set<state_t> erase;

  //collect unreachable
  for (auto &it : scci.unreachable) {
      erase.emplace(it);
  }
  scci.unreachable.clear();

  //collect dead
  for (auto &it : ba.adj) {
    if (scci.scc.find(it.first) != end(scci.scc)) { //has assigned scc
        auto scit = scci.scc.at(it.first);
        if (scci.dead.at(scit)) { //is a dead scc
          erase.emplace(it.first);
          scci.scc.erase(it.first);
        }
    }
  }
  for (auto &kv : scci.dead) {
    if (kv.second) {
      scci.accepting.erase(kv.first);
      scci.rejecting.erase(kv.first);
      scci.trivial.erase(kv.first);
      scci.sccrep.erase(kv.first);
    }
  }
  scci.dead.clear();

  std::vector<state_t> erasevec;
  copy(begin(erase),end(erase),back_inserter(erasevec));
  erase.clear();

  kill_states(ba, erasevec);
  return erasevec.size();
}

}

