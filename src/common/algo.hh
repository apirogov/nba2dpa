#pragma once

#include <iostream>
#include "common/util.hh"
#include "common/part_refinement.hh"
#include "common/scc.hh"

#include <range/v3/all.hpp>

#include <functional>
#include <queue>
#include <stack>
#include <set>
#include <map>

namespace nbautils {
using namespace std;
using namespace nbautils;

template <typename Node,typename Sym>
using outsym_fun = function<vector<Sym>(Node)> const&;

template <typename Node,typename Sym>
using succ_sym_fun = function<vector<Node>(Node,Sym)> const&;

struct BaSccClassification {
  set<unsigned> accepting;  // tags an scc as fully accepting
  set<unsigned> rejecting;  // tags an scc as fully rejecting
};

// classify sccs as accepting or rejecting
template <typename Node>
BaSccClassification ba_classify_sccs(SCCDat<Node> const& scci,
    function<bool(Node)> const& is_acc) {
  BaSccClassification ret;

  auto const conj = [](bool a, bool b){return a&&b;};
  ret.accepting = mapbool_to_set(fold_sccs<Node,bool>(scci, true, is_acc, conj));
  ret.rejecting = mapbool_to_set(fold_sccs<Node,bool>(scci, true,
        [&is_acc](Node v){ return !is_acc(v); }, conj));

  return ret;
}

//returns sorted list of accepting sinks (acc. states with self-loop for each sym)
template <typename Node, typename Sym>
vector<Node> ba_get_accepting_sinks(vector<Node> const& states, int num_syms,
    function<bool(Node)> is_acc,
    outsym_fun<Node,Sym> get_outsyms,
    succ_sym_fun<Node,Sym> get_xsuccs) {
  vector<Node> ret;
  for (auto const& v : states) {
    auto outsyms = get_outsyms(v);
    //must be accepting and have successors for each symbol
    bool accsink = is_acc(v) && (int)outsyms.size() == num_syms;
    for (auto i=0; i<num_syms; i++) {
        // must have self-loop for each symbol
        auto xsucs = get_xsuccs(v, i);
        if (!contains(xsucs, v))
          accsink = false;
    }
    if (accsink)
      ret.push_back(v);
  }
  return ret;
}

inline void mark_dead_sccs(set<unsigned> const& rejecting, set<unsigned> const& trivial,
    succ_scc_fun const& get_suc_sccs, map<unsigned,bool>& dead, unsigned num) {
  if (map_has_key(dead, num))  // done already
    return;

  // std::cout << "scc " << num << std::endl;
  auto sucsccs = get_suc_sccs(num);
  // mark children first
  for (auto sucnum : sucsccs)
    mark_dead_sccs(rejecting, trivial, get_suc_sccs, dead, sucnum);

  // if we are rejecting and trivial, assume we're dead
  bool isdead = contains(rejecting, num) || contains(trivial, num);
  // check children and try to falsify
  for (auto sucscc : sucsccs) isdead = isdead && dead[sucscc];
  // if still dead, we're really dead.
  dead[num] = isdead;
}

// run dfs that marks dead sccs (assuming sccs 0..num_sccs-1 exist)
inline set<unsigned> ba_get_dead_sccs(int num_sccs,
    set<unsigned> const& rejecting, set<unsigned> const& trivial,
    succ_scc_fun const& get_succs) {
  map<unsigned, bool> dead;
  for (int i=0; i<num_sccs; i++)
    mark_dead_sccs(rejecting, trivial, get_succs, dead, i);

  set<unsigned> ret;
  for (auto const& it : dead)
    if (it.second)
      ret.emplace(it.first);
  return ret;
}


// ----------------------------------------------------------------------------

template <typename T>
int max_chain(function<int(T)> const& oldpri, map<T,int>& newpri,
    vector<T> const& p, succ_fun<T> const& get_succs) {
  if (p.empty()) //by definition, empty set has no chain
    return 0;

  int maxlen = 0;
  // cout << "max_chain " << seq_to_str(p) << endl;

  // maximal essential subsets = non-trivial SCCs in restricted graph
  succ_fun<T> succs_in_p = [&](auto v) { return set_intersect(get_succs(v), p); };
  auto scci = get_sccs(p, succs_in_p, const_true);
  auto triv = trivial_sccs(scci, succs_in_p);

  // lift priority map to state sets
  auto strongest_oldpri = fold_sccs<T,int>(scci, 0, oldpri, [](auto a, auto b){return max(a,b);});

  int i=0;
  for (auto const& scc : scci.sccs) {
    // cout << "scc " << seq_to_str(scc) << endl;

    if (contains(triv,i)) {
      // cout << "trivial, skip" << endl;

      ++i;
      continue; //non-essential
    }

    auto const scc_pri = strongest_oldpri.at(i);

    // cout << "scc pri " << scc_pri << endl;

    auto deriv_scc = vec_filter(scc, [&](auto v){return oldpri(v) < scc_pri;});
    vec_to_set(deriv_scc);
    auto const not_deriv_scc = vec_filter(scc, [&](auto v){return oldpri(v) >= scc_pri;});

    int m = 0;
    if (scc_pri > 0) {
      m = max_chain(oldpri, newpri, deriv_scc, succs_in_p);

      if ((scc_pri - m) % 2 == 1) //alternation -> requires new priority
        m++;
    }

    for (auto s : not_deriv_scc) {
      newpri[s] = m;
    }

    maxlen = max(maxlen, m);

    ++i;
  }
  // cout << "max_chain " << seq_to_str(p)  << " done" << endl;
  return maxlen;
}

// takes priorities (with max odd parity!), minimizes them
// takes all states, normal successor function and old priority map function
// priorities must be from a max odd acceptance
// returns new state to priority map
// Paper: "Computing the Rabin Index of a parity automaton"
template <typename T>
map<T, int> pa_minimize_priorities(vector<T> const& states,
    succ_fun<T> const& get_succs, function<int(T)> const& oldpri) {
  map<T, int> primap;
  for (auto const s : states)
    primap[s] = 0;
  max_chain(oldpri, primap, states, get_succs);
  return primap;
}


template <typename Node,typename Sym>
using det_succ_sym_fun = function<Node(Node,Sym)> const&;
template <typename Node,typename Val>
using node_prop_fun = function<Val(Node)> const&;

// https://en.wikipedia.org/wiki/DFA_minimization
//
// hopcroft algorithm to calculate equivalence classes that can be merged
// without changing the language / output behaviour
// requires complete, deterministic automaton
template <typename T, typename S, typename C>
vector<vector<T>> dfa_equivalent_states(vector<T> const& states, node_prop_fun<T,C> color,
    int num_syms, det_succ_sym_fun<T,S> get_xsucc) {
  // prepare initial partitions = different colors
  vector<T> sorted = states;
  ranges::sort(sorted, [&](T a, T b){ return color(a) < color(b); });
  vector<vector<T>> startsets = sorted
    | ranges::view::group_by([&](T a, T b){ return color(a) == color(b); });

  // cerr << "starter sets:" << endl;
  // for (auto  s : startsets) {
  //   cerr << seq_to_str(s) << endl;
  // }
  // cerr << endl;

  PartitionRefiner<T> p(startsets);
  auto w=p.get_set_ids();

  // cerr << "starting refinement" << endl;

  while (!w.empty()) {
    auto const a = w.back();
    // cerr << "separator " << seq_to_str(p.get_elements_of(a)) << endl;
    w.pop_back();

    for (auto i=0; i<num_syms; i++) {
      auto const succ_in_a = [&](T st){ return p.get_set_of(get_xsucc(st, i)) == a; };

      for (auto y : p.get_set_ids()) {
        auto const z = p.separate(y, succ_in_a);
        if (z) { //separation happened
          if (contains(w, y)) //if y is in w, its symbol still is, and we need the other set
            w.push_back(*z);
          else { //if y is not in w, take smaller part as separator
            if (p.get_set_size(y) <= p.get_set_size(*z)) {
              w.push_back(y);
            } else {
              w.push_back(*z);
            }
          }
        }
      }
    }
  }

  return p.get_refined_sets();
}

}
