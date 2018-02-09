#pragma once

#include "interfaces.hh"

#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace nbautils {
using namespace std;

// type for "small" automata (as the input should be)
using small_state_t = uint8_t;
constexpr size_t max_nba_states = 1 << (8 * sizeof(small_state_t));

// type for alphabet (<=8 aps with all <=2^8=256 combinations)
using sym_t = uint8_t;
constexpr size_t max_nba_syms = 1 << (8 * sizeof(sym_t));

// type for "big" automata
using state_t = uint32_t;

struct SWAMeta {
  // automaton name
  string name;
  // atomic props
  vector<string> aps;

  size_t num_syms = 0;
};

// acceptance-labelled KS with tagged nodes
template <typename L, typename T,
          class tag_storage = generic_trie_bimap<T, small_state_t, state_t>>
struct SWA {
  typedef unique_ptr<SWA<L, T, tag_storage>> uptr;
  typedef shared_ptr<SWA<L, T, tag_storage>> sptr;

  SWAMeta meta;

  state_t init;  // initial state

  // list of successors per symbol
  map<state_t, map<sym_t, vector<state_t>>> adj;
  // state acceptance label
  map<state_t, L> acc;
  // shared_ptr<bimap<state_t,T>> tag;
  unique_ptr<bimap<T, state_t, tag_storage>> tag;

  bool has_state(state_t const& s) const { return adj.find(s) != end(adj); }

  bool has_acc(state_t const& s) const { return acc.find(s) != end(acc); }

  vector<state_t> succ(state_t const& p, sym_t const& x) const {
    if (!has_state(p)) return {};
    auto const& edges = adj.at(p);
    if (edges.find(x) == end(edges)) return {};
    return edges.at(x);
  }

  vector<sym_t> outsyms(state_t const& p) const {
    if (adj.find(p) == end(adj)) return {};
    auto const& edges = adj.at(p);
    vector<sym_t> ret;
    for (auto const& it : edges) ret.push_back(it.first);
    return ret;
  }

  vector<state_t> succ(state_t const& p) const {
    if (adj.find(p) == end(adj)) return {};
    auto const& edges = adj.at(p);
    // collect successors for any symbol
    set<state_t> s;
    for (auto const& it : edges) {
      copy(begin(it.second), end(it.second), inserter(s, end(s)));
    }
    vector<state_t> ret;
    copy(begin(s), end(s), back_inserter(ret));
    return ret;
  }

  // given set in sorted(!!) vector, kill states
  void remove_states(vector<state_t> const& tokill) {
    // TODO: refuse to kill initial
    for (auto const& it : tokill) {
      auto el = adj.find(it);
      if (el != end(adj)) adj.erase(el);  // kill state + edges
    }
    for (auto const& it : tokill) {
      auto el = acc.find(it);
      if (el != end(acc)) acc.erase(el);  // kill acc. mark
    }

    // kill states from edge targets
    for (auto& it : adj) {
      for (auto& jt : it.second) {
        vector<state_t> res;
        set_difference(begin(jt.second), end(jt.second), begin(tokill), end(tokill),
                       back_inserter(res));
        jt.second = res;
      }
    }
  }

  // given set in sorted(!!) vector, merge states into one
  // TODO: somehow make eqclass const but still without copy?
  void merge_states(SWA<L, T>& aut, vector<state_t>& eqclass) {
    if (eqclass.size() < 2) return;

    // get representative
    state_t rep = eqclass.back();

    // add outgoing edges from all class members to representative
    for (auto const& st : eqclass) {
      for (auto& symmap : aut.adj.at(st)) {
        auto& sym = symmap.first;
        auto& sucs = symmap.second;

        vector<state_t> res;
        set_union(begin(sucs), end(sucs), begin(aut.adj.at(rep)[sym]),
                  end(aut.adj.at(rep)[sym]), back_inserter(res));
        aut.adj.at(rep).at(sym) = res;
      }
    }
    // add ingoing edges into all class members to representative
    for (auto& it : aut.adj) {
      auto& st = it.first;
      for (auto& symmap : aut.adj.at(st)) {
        auto& sym = symmap.first;
        auto& sucs = symmap.second;
        vector<state_t> res;
        set_intersection(begin(sucs), end(sucs), begin(eqclass), end(eqclass),
                         back_inserter(res));
        if (res.size() > 0) {
          res = vector<state_t>{};
          auto singl = std::vector<state_t>{rep};
          set_union(begin(sucs), end(sucs), begin(singl), end(singl), back_inserter(res));
          sucs = res;
        }
      }
    }

    // remove representative from kill list
    eqclass.pop_back();
    // kill the others!!
    aut.remove_states(eqclass);
  }
};

template <typename L, typename T, typename S>
vector<small_state_t> powersucc(SWA<L, T> const& ks, std::vector<S> const& ps,
                                sym_t const& x) {
  std::set<state_t> suc;
  for (auto const& p : ps) {
    auto qs = ks.succ(p, x);
    std::copy(qs.begin(), qs.end(), std::inserter(suc, suc.end()));
  }
  return vector<S>(suc.begin(), suc.end());
}

using BA = SWA<bool, bool>;

}  // namespace nbautils
