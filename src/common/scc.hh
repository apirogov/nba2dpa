#pragma once

#include <stack>
#include <set>
#include <map>

namespace nbautils {

using namespace std;
using namespace nbautils;

struct SCCDat {
  std::map<state_t, unsigned> scc_of;            //state to scc
  std::map<unsigned, std::vector<state_t>> sccs; //scc to states
};

// https://en.wikipedia.org/wiki/Path-based_strong_component_algorithm with extensions
// takes list of all states of graph we want to have an scc for
// a function that supplies successors of a state
// performs an SCC DFS traversal.
// returns list of SCCs such that later SCCs can not reach earlier SCCs
template <typename F>
SCCDat get_sccs(std::vector<state_t> const& states, F get_succs, bool store_partitions=true) {
  SCCDat ret;

  std::stack<state_t> call;  // dfs call stack
  std::stack<state_t> reps;  // scc representative stack
  std::stack<state_t> open;  // not yet fully completed vertex stack

  int count = 0;
  std::map<state_t, unsigned> order;  // first visit order

  // schedule all states to be called
  for (auto const& v : states)
    call.push(v);

  while (!call.empty()) {
    auto const v = call.top();

    if (!map_has_key(order,v)) {  // dfs just "called" with current node
      // std::cout << "discover " << v << std::endl;

      order[v] = count++;         // assign visit order
      reps.push(v);  // SCC representative candidate (for now)
      open.push(v);  // this node is "pending" (not completely discovered from here)

      for (auto const w : get_succs(v)) {  // process edges with any label
        if (!map_has_key(order, w)) {
          call.push(w);  // recursively explore nodes that have not been visited yet
        } else if (!map_has_key(ret.scc_of, w)) {
          // if already visited, but not with assigned scc, we have found a loop
          // -> drop candidates, keep oldest on this loop as SCC representative
          while (order.at(reps.top()) > order.at(w)) reps.pop();
        }
      }

    } else {
      if (map_has_key(ret.scc_of, v)) {
        //this node is already completed and uselessly visited
        call.pop();
        continue;
      }

      // returned from recursive calls
      // std::cout << "return to " << v << std::endl;

      // is still rep. -> we found an SCC
      if (reps.top() == v) {

        std::vector<state_t> scc_states;
        // drop states up to the current state, they are done and part of the SCC
        state_t tmp;
        unsigned const curnum = ret.sccs.size();
        do {
          tmp = open.top();
          open.pop();
          ret.scc_of[tmp] = curnum;

          if (store_partitions)
            scc_states.push_back(tmp);
        } while (tmp != v);
        if (store_partitions) {
          sort(begin(scc_states), end(scc_states));
          ret.sccs[curnum] = {};
          swap(ret.sccs[curnum], scc_states);
        }

        // bool go_on = scc_handler(*ret);
        // if (!go_on) //can be used to abort on-the-fly scc search
        //   return ret;

        // current SCC is done
        reps.pop();
      }

      // std::cout << "done " << v << std::endl;
      call.pop();  // current node done -> dfs "returns" to previous caller
    }
  }

  return ret;
}

// takes SCC information and an SCC, returns successor SCCs
template <typename F>
std::set<unsigned> succ_sccs(F const& succ, SCCDat const& scci, unsigned const& num) {
  std::set<unsigned> sucsccs;
  for (auto const st : scci.sccs.at(num)) {
    for (auto const sucst : succ(st)) {
      auto const sucscc = scci.scc_of.at(sucst);
      if (sucscc != num)
        sucsccs.emplace(sucscc);
    }
  };
  return sucsccs;
}

// return trivial SCCs (no edge leads back to SCC)
template <typename F>
std::set<unsigned> trivial_sccs(F const& succ, SCCDat const& scci) {
  set<unsigned> ret;
  for (auto const& scc : scci.sccs) {
    auto const& sts = scc.second;
    // single state with no self-loop?
    bool const trivacc = sts.size() == 1;
    bool const noselfloop = !sorted_contains(succ(sts.front()), sts.front());
    if (trivacc && noselfloop)
      ret.emplace(scc.first);
  }
  return ret;
}

}
