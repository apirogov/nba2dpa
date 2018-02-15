#pragma once

#include "bench.hh"
#include "types.hh"

#include <algorithm>
#include <iostream>
#include <queue>
#include <set>
#include <stack>
#include <unordered_set>
#include <utility>
#include <vector>

#include <spdlog/spdlog.h>
namespace spd = spdlog;

namespace nbautils {
using namespace std;

using scc_t = state_t;
using scc_flag = unordered_set<scc_t>;

struct SCCInfo {
  typedef std::shared_ptr<SCCInfo> sptr;
  typedef std::unique_ptr<SCCInfo> uptr;

  map<state_t, scc_t> scc;   // tags each state with scc num
  set<state_t> unreachable;  // tags a state as unreachable from initial state

  std::map<scc_t, state_t> sccrep;  // representative state
  scc_flag accepting;               // tags an scc as fully accepting
  scc_flag rejecting;               // tags an scc as fully rejecting
  scc_flag trivial;                 // tags an scc as trivial (single state, no loop)
  std::map<scc_t, bool> dead;  // tags an scc as dead (means that all reachable states are
                               // rejecting)
};

/*
template<typename L,typename T>
std::vector<state_t> scc_states(SWA<L,T> const& aut, SCCInfo const& scci, scc_t const&
num) { std::queue<state_t> bfsq; std::set<state_t> visited;
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

template <typename L, typename T, typename S>
std::vector<state_t> succ_sccs(SWA<L, T, S> const& aut, SCCInfo const& scci,
                               scc_t const& num) {
  std::queue<state_t> bfsq;
  std::set<state_t> visited;
  std::set<nbautils::scc_t> sucsccs;

  bfsq.push(scci.sccrep.at(num));
  while (!bfsq.empty()) {
    auto const st = bfsq.front();
    bfsq.pop();
    if (visited.find(st) != end(visited)) continue;  // have visited this one
    visited.emplace(st);

    for (auto& sucst : aut.succ(st)) {
      auto sucscc = scci.scc.at(sucst);
      if (sucscc == num)
        bfsq.push(sucst);
      else
        sucsccs.emplace(sucscc);
    }
  }
  std::vector<nbautils::scc_t> ret;
  copy(begin(sucsccs), end(sucsccs), back_inserter(ret));
  return ret;
}

template <typename L, typename T, typename S>
void mark_dead_sccs(SWA<L, T, S> const& aut, SCCInfo& scci, scc_t num) {
  if (scci.dead.find(num) != end(scci.dead))  // done already
    return;

  // std::cout << "scc " << num << std::endl;
  auto sucsccs = succ_sccs(aut, scci, num);
  // mark children first
  for (auto sucscc : sucsccs) mark_dead_sccs(aut, scci, sucscc);

  // if we are rejecting and trivial, assume we're dead
  bool dead = (scci.rejecting.find(num) != end(scci.rejecting) ||
               scci.trivial.find(num) != end(scci.trivial));
  // check children and try to falsify
  for (auto sucscc : sucsccs) dead = dead && scci.dead[sucscc];
  // if still dead, we're really dead.

  scci.dead[num] = dead;
}

// https://en.wikipedia.org/wiki/Path-based_strong_component_algorithm with extensions
template <typename L, typename T, typename S>
SCCInfo::uptr get_scc_info(SWA<L, T, S> const& aut, bool analyse_acc = true) {
  spd::get("log")->debug("get_scc_info(aut,{})", analyse_acc);
  auto starttime = get_time();

  auto sccip = std::make_unique<SCCInfo>(SCCInfo());
  auto& scci = *sccip;

  std::stack<state_t> call;  // dfs call stack
  std::stack<state_t> reps;  // scc representative stack
  std::stack<state_t> open;  // not yet fully completed vertex stack

  std::map<state_t, unsigned> order;  // first visit order
  int count = 0;

  // we just start with the initial state
  call.push(aut.init);
  while (!call.empty()) {
    auto v = call.top();
    if (order.find(v) == order.end()) {  // dfs just "called" with current node
      order[v] = count++;                // assign visit order

      reps.push(v);  // SCC representative candidate (for now)
      open.push(v);  // this node is "pending" (not completely discovered from here)

      for (auto w : aut.succ(v)) {  // process edges with any label
        if (order.find(w) == order.end()) {
          call.push(w);  // recursively explore nodes that have not been visited yet
        } else if (scci.scc.find(w) == scci.scc.end()) {
          // if already visited, but not with assigned scc, we have found a loop
          // -> drop candidates, keep oldest on this loop as SCC representative
          while (order[reps.top()] > order[w]) reps.pop();
        }
      }

    } else {  // returned from recursive calls
      // is still rep. -> we found an SCC
      if (reps.top() == v) {
        // assume SCC is accepting and rejecting, then falsify
        bool accscc = true;
        bool rejscc = true;
        // count how many state this scc has to detect triv. SCC
        int scc_sz = 0;
        auto scc_num = scci.sccrep.size();  // size corresponds to next number

        // drop states up to the current state, they are done and part of the SCC
        state_t tmp;
        do {
          tmp = open.top();
          open.pop();

          // if acceptance info given, detect acc./rej. sccs
          if (analyse_acc) {
            auto actmp = aut.acc.find(tmp) != end(aut.acc);
            accscc = accscc && actmp;
            rejscc = rejscc && !actmp;
          }

          // assign current state current SCC
          scci.scc[tmp] = scc_num;

          scc_sz++;
        } while (tmp != v);

        // now lets get the extra information we want
        if (analyse_acc) {
          // mark completed SCC as accepting/rejecting
          if (accscc) scci.accepting.emplace(scc_num);
          if (rejscc) scci.rejecting.emplace(scc_num);

          // trivial scc with no self-loop?
          bool trivacc = scc_sz == 1;
          bool noselfloop = true;
          auto vsucs = aut.succ(v);
          if (find(begin(vsucs), end(vsucs), tmp) != end(vsucs)) noselfloop = false;
          if (trivacc && noselfloop) scci.trivial.emplace(scc_num);
        }

        // store the representative state (from which the rest can be obtained easily)
        scci.sccrep[scc_num] = reps.top();

        // current SCC is done
        reps.pop();
      }

      call.pop();  // current node done -> dfs "returns" to previous caller
    }
  }
  // states still without assigned scc are unreachable
  for (auto& it : aut.adj)
    if (scci.scc.find(it.first) == end(scci.scc)) scci.unreachable.emplace(it.first);

  // run dfs that marks dead sccs
  mark_dead_sccs(aut, scci, scci.scc.at(aut.init));

  spd::get("log")->debug("get_scc_info completed ({:.4f} s)", get_secs_since(starttime));

  return move(sccip);
}

// use SCC info to purge unreachable and dead states from automaton and scc info
template <typename L, typename T, typename S>
size_t trim_ba(SWA<L, T, S>& ba, SCCInfo& scci) {
  set<state_t> erase;

  // collect unreachable
  move(begin(scci.unreachable), end(scci.unreachable), inserter(erase,end(erase)));

  // collect dead states
  for (auto& it : ba.adj) {
    if (scci.scc.find(it.first) != end(scci.scc)) {  // has assigned scc
      auto scit = scci.scc.at(it.first);
      if (scci.dead.at(scit)) {  // is a dead scc
        erase.emplace(it.first);
        scci.scc.erase(it.first);
      }
    }
  }
  // unmark these states
  for (auto& kv : scci.dead) {
    if (kv.second) {
      scci.accepting.erase(kv.first);
      scci.rejecting.erase(kv.first);
      scci.trivial.erase(kv.first);
      scci.sccrep.erase(kv.first);
    }
  }
  scci.dead.clear();

  vector<state_t> erasevec;
  move(begin(erase), end(erase), back_inserter(erasevec));

  ba.remove_states(erasevec);
  return erasevec.size();
}

}  // namespace nbautils
