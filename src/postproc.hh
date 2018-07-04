#pragma once

#include <iostream>
#include "common/util.hh"
#include "common/part_refinement.hh"
#include "common/scc.hh"

// #include <range/v3/all.hpp>

#include <functional>
#include <memory>
#include <set>
#include <map>

namespace nbautils {
using namespace std;
using namespace nbautils;

template <typename Node,typename Sym>
using outsym_fun = function<vector<Sym>(Node)>;

template <typename Node,typename Sym>
using succ_sym_fun = function<vector<Node>(Node,Sym)>;

template <typename Node>
using succ_fun = function<vector<Node>(Node)>;

// ----------------------------------------------------------------------------

template <typename Range, typename T>
int max_chain(function<int(T)> const& oldpri, map<T,int>& newpri,
    Range const& p, succ_fun<T> const& get_succs) {
  if (p.empty()) //by definition, empty set has no chain
    return 0;

  int maxlen = 0;
  // cout << "max_chain " << seq_to_str(p) << endl;

  // maximal essential subsets = non-trivial SCCs in restricted graph
  succ_fun<T> succs_in_p = [&](auto v) {
    auto const tmp = get_succs(v);
    vector<T> ret;
    ranges::set_intersection(ranges::view::all(tmp), p, ranges::back_inserter(ret));
    return ret;
  };
  auto scci = get_sccs(p, succs_in_p);
  auto triv = trivial_sccs(succs_in_p, scci);

  // lift priority map to state sets
  map<unsigned, int> strongest_oldpri;
  for (auto const& it : scci.sccs) {
    int maxp = -1;
    for (auto const& st : it.second) {
      maxp = max(maxp, oldpri(st));
    }
    strongest_oldpri[it.first] = maxp;
  }

  for (auto const& it : scci.sccs) {
    int const i = it.first;
    vector<state_t> const& scc = it.second;
    // cout << "scc " << seq_to_str(scc) << endl;

    if (contains(triv,it.first)) {
      continue;
      // cout << "trivial, skip" << endl;
    }

    pri_t const scc_pri = strongest_oldpri[i];
    // cout << "scc pri " << scc_pri << endl;

    /*
    auto deriv_scc = vec_filter(scc, [&](state_t v){return oldpri(v) < scc_pri;});
    vec_to_set(deriv_scc);
    auto const not_deriv_scc = vec_filter(scc, [&](state_t v){return oldpri(v) >= scc_pri;});
    */
    auto const deriv_scc = ranges::view::all(scc)
      | ranges::view::filter([&](state_t v){return oldpri(v) < scc_pri;}) | ranges::to_vector;
    auto const not_deriv_scc = ranges::view::all(scc)
      | ranges::view::filter([&](state_t v){return oldpri(v) >= scc_pri;}) | ranges::to_vector;

    int m = 0;
    if (scc_pri > 0) {
      m = max_chain(oldpri, newpri, deriv_scc, succs_in_p);

      if ((scc_pri - m) % 2 == 1) //alternation -> requires new priority
        m++;
    }

    for (auto const s : not_deriv_scc) {
      newpri[s] = m;
    }

    maxlen = max(maxlen, m);
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

// ----------------------------------------------------------------------------

template <typename Node,typename Sym>
using det_succ_sym_fun = function<Node(Node,Sym)>;
template <typename Node,typename Val>
using node_prop_fun = function<Val(Node)>;

// https://en.wikipedia.org/wiki/DFA_minimization
//
// hopcroft algorithm to calculate equivalence classes that can be merged
// without changing the language / output behaviour
// requires complete, deterministic automaton
template <typename T, typename S, typename C>
vector<vector<T>> dfa_equivalent_states(vector<T> const& states, node_prop_fun<T,C> color,
    int num_syms, det_succ_sym_fun<T,S> get_xsucc) {
  vector<T> sorted = states;

  // prepare initial partitions = different colors
  // ranges::sort(sorted, [&](T a, T b){ return color(a) < color(b); });
  sort(begin(sorted), end(sorted), [&](T a, T b){ return color(a) < color(b); });

  // vector<vector<T>> startsets = sorted
  //   | ranges::view::group_by([&](T a, T b){ return color(a) == color(b); });
  vector<vector<T>> startsets = group_by(sorted, [&](T a, T b){ return color(a) == color(b); });

  // cerr << "starter sets:" << endl;
  // for (auto  s : startsets) {
  //   cerr << seq_to_str(s) << endl;
  // }

  PartitionRefiner<T> p(startsets);

  // set seems not to make a difference :/
  /*
  auto const wvec=p.get_set_ids();
  auto const sym_set_cmp = [](auto const& a, auto const& b) {
        return lexicographical_compare(a->first, a->second, b->first, b->second);
      };
  auto w=set<typename PartitionRefiner<T>::sym_set, decltype(sym_set_cmp)>(cbegin(wvec), cend(wvec), sym_set_cmp);
  */

  auto w=p.get_set_ids();

  while (!w.empty()) {
    auto const a = w.back(); w.pop_back();
    // auto const a = *w.begin(); w.erase(w.begin());

    auto const sepset = p.get_elements_of(a); //need to take a copy, as it is modified in loop

    for (auto i=0; i<num_syms; i++) {
      auto const succ_in_a = [&](T st){ return sorted_contains(sepset, get_xsucc(st, i)); };

      for (auto y : p.get_set_ids()) {
        //TODO: try to precalculate sep X subset of Y set instead of using predicate?
        auto const z = p.separate(y, succ_in_a);

        if (z) { //separation happened
          //TODO: w is vector, this is slow
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
