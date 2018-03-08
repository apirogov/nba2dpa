#pragma once

#include "common/util.hh"
#include "common/graph.hh"

#include <functional>
#include <memory>
#include <queue>
#include <stack>
#include <set>
#include <map>

namespace nbautils {

using namespace std;
using namespace nbautils;

using succ_scc_fun = succ_fun<unsigned>;

template <typename Node>
struct SCCDat {
  using uptr = std::unique_ptr<SCCDat<Node>>;

  std::vector<std::vector<Node>> sccs;
  std::map<Node, unsigned> scc_of;
};

// https://en.wikipedia.org/wiki/Path-based_strong_component_algorithm with extensions
// takes list of all states of graph we want to have an scc for
// a function that supplies successors of a state
// a function that is called for each new found SCC and returns whether search has to go on
// performs an SCC DFS traversal.
// returns list of SCCs such that later SCCs can not reach earlier SCCs
template <typename Node, typename F, typename G>
typename SCCDat<Node>::uptr get_sccs(std::vector<Node> const& states, F get_succs, G scc_handler) {
  auto ret = make_unique<SCCDat<Node>>();

  std::stack<Node> call;  // dfs call stack
  std::stack<Node> reps;  // scc representative stack
  std::stack<Node> open;  // not yet fully completed vertex stack

  int count = 0;
  std::map<Node, unsigned> order;  // first visit order

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

      auto const sucs = get_succs(v);

      for (auto const w : sucs) {  // process edges with any label
        if (!map_has_key(order, w)) {
          call.push(w);  // recursively explore nodes that have not been visited yet
        } else if (!map_has_key(ret->scc_of, w)) {
          // if already visited, but not with assigned scc, we have found a loop
          // -> drop candidates, keep oldest on this loop as SCC representative
          while (order.at(reps.top()) > order.at(w)) reps.pop();
        }
      }

    } else {
      if (map_has_key(ret->scc_of, v)) {
        //this node is already completed and uselessly visited
        call.pop();
        continue;
      }

      // returned from recursive calls
      // std::cout << "return to " << v << std::endl;

      // is still rep. -> we found an SCC
      if (reps.top() == v) {
        // drop states up to the current state, they are done and part of the SCC
        std::vector<Node> scc_states;
        Node tmp;
        do {
          tmp = open.top();
          open.pop();
          ret->scc_of[tmp] = ret->sccs.size();
          scc_states.push_back(tmp);
        } while (tmp != v);

        ret->sccs.push_back({});
        swap(ret->sccs.back(), scc_states);

        bool go_on = scc_handler(*ret);
        if (!go_on) //can be used to abort on-the-fly scc search
          return ret;

        // current SCC is done
        reps.pop();
      }

      // std::cout << "done " << v << std::endl;
      call.pop();  // current node done -> dfs "returns" to previous caller
    }
  }

  return std::move(ret);
}

template <typename Node, typename F>
typename SCCDat<Node>::uptr get_sccs(std::vector<Node> const& states, F get_succs) {
  return get_sccs(states, get_succs, const_true);
}

// takes SCC information and an SCC, returns successor SCCs
template <typename Node>
std::vector<unsigned> succ_sccs(SCCDat<Node> const& scci, unsigned const& num,
    succ_fun<Node> const& get_suc_st) {
  std::set<unsigned> sucsccs;
  bfs(scci.sccs.at(num).front(), [&](auto const& st, auto const& visit, auto const&) {
    for (auto const& sucst : get_suc_st(st)) {
      auto sucscc = scci.scc_of.at(sucst);
      if (sucscc == num)
        visit(sucst);
      else
        sucsccs.emplace(sucscc);
    }
  });
  return vector<unsigned>(cbegin(sucsccs), cend(sucsccs));
}

template <typename Node>
std::set<unsigned> trivial_sccs(SCCDat<Node> const& scci, succ_fun<Node> const& get_succs) {
  set<unsigned> ret;
  unsigned i=0;
  for (auto const& scc : scci.sccs) {
    // single state with no self-loop?
    bool const trivacc = scc.size() == 1;
    bool const noselfloop = !contains(get_succs(scc.front()), scc.front());
    if (trivacc && noselfloop)
      ret.emplace(i);

    ++i;
  }
  return ret;
}

// given sccs and map of states to semigroup together with semigroup operation
// and default value, return map from scc to semigroup element
// (used to get accepting, rejecting sccs of NBA and dominant SCC priority of DPAs)
template <typename Node, typename V>
map<unsigned, V> fold_sccs(SCCDat<Node> const& scci,
    V def, function<V(Node)> lift, function<V(V,V)> op) {
  map<unsigned, V> ret;
  unsigned i=0;
  for (auto const& scc : scci.sccs) {
    ret[i] = vec_fold_mapped(scc, def, lift, op);
    ++i;
  }
  return ret;
}

}
