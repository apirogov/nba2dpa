#pragma once

#include "swa.hh"
#include "common/algo.hh"

#include <cassert>
#include <algorithm>
#include <iostream>
#include <queue>
#include <set>
#include <stack>
#include <unordered_set>
#include <utility>
#include <vector>

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
  std::map<scc_t, size_t> sccsz;    // size of SCC
  inline size_t num_sccs() { return sccrep.size(); }

  scc_flag accepting;               // tags an scc as fully accepting
  scc_flag rejecting;               // tags an scc as fully rejecting
  scc_flag trivial;                 // tags an scc as trivial (single state, no loop)
};

//states in given scc
template<Acceptance A,typename T,typename S>
std::vector<state_t> scc_states(SWA<A,T,S> const& aut, SCCInfo const& scci, scc_t const& num) {
  std::vector<nbautils::scc_t> ret;
  bfs(scci.sccrep.at(num), [&](auto const& st, auto const& visit, auto const&) {
    ret.push_back(st);
    for (auto const& sucst : aut.succ(st)) {
      if (scci.scc.at(sucst)==num) {
        visit(sucst);
      }
    }
  });
  sort(begin(ret),end(ret));
  return move(ret);
}

//successor sccs of given scc (by bfs of that scc)
template <Acceptance A, typename T, typename S>
std::vector<state_t> succ_sccs(SWA<A, T, S> const& aut, SCCInfo const& scci,
                               scc_t const& num) {
  std::set<nbautils::scc_t> sucsccs;
  bfs(scci.sccrep.at(num), [&](auto const& st, auto const& visit, auto const&) {
    for (auto const& sucst : aut.succ(st)) {
      auto sucscc = scci.scc.at(sucst);
      if (sucscc == num)
        visit(sucst);
      else
        sucsccs.emplace(sucscc);
    }
  });
  return vector<nbautils::scc_t>(cbegin(sucsccs), cend(sucsccs));
}

// https://en.wikipedia.org/wiki/Path-based_strong_component_algorithm with extensions
template <Acceptance A, typename T>
SCCInfo::uptr get_scc_info(SWA<A, T> const& aut, bool analyse_acc = false) {
  assert(!analyse_acc || A==Acceptance::BUCHI); //extended SCC analysis only for BA

  auto sccip = std::make_unique<SCCInfo>(SCCInfo());
  auto& scci = *sccip;

  std::stack<state_t> call;  // dfs call stack
  std::stack<state_t> reps;  // scc representative stack
  std::stack<state_t> open;  // not yet fully completed vertex stack

  std::map<state_t, unsigned> order;  // first visit order
  int count = 0;

  // we just start with the initial state(s)
  for (auto &v : aut.get_init())
    call.push(v);

  while (!call.empty()) {
    auto v = call.top();
    if (!map_has_key(order,v)) {  // dfs just "called" with current node
      order[v] = count++;                // assign visit order

      reps.push(v);  // SCC representative candidate (for now)
      open.push(v);  // this node is "pending" (not completely discovered from here)

      for (auto w : aut.succ(v)) {  // process edges with any label
        if (!map_has_key(order, w)) {
          call.push(w);  // recursively explore nodes that have not been visited yet
        } else if (!map_has_key(scci.scc, w)) {
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
        int scc_sz = 0; // count how many state this scc has
        auto scc_num = scci.sccrep.size();  // size of scc set corresponds to next number

        // drop states up to the current state, they are done and part of the SCC
        state_t tmp;
        do {
          tmp = open.top();
          open.pop();

          // if acceptance info given, detect acc./rej. sccs
          if (analyse_acc) {
            auto actmp = aut.has_accs(tmp);
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
          bool const trivacc = scc_sz == 1;
          bool const noselfloop = !contains(aut.succ(v), tmp);
          if (trivacc && noselfloop) scci.trivial.emplace(scc_num);
        }

        // store the representative state (from which the rest can be obtained easily) + size
        scci.sccrep[scc_num] = reps.top();
        scci.sccsz[scc_num] = scc_sz;

        // current SCC is done
        reps.pop();
      }

      call.pop();  // current node done -> dfs "returns" to previous caller
    }
  }
  // states still without assigned scc are unreachable
  for (auto const& it : aut.states()) //TODO: iterator
    if (!map_has_key(scci.scc, it)) scci.unreachable.emplace(it);

  return move(sccip);
}

}  // namespace nbautils
