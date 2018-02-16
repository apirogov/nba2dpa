#pragma once

#include "interfaces.hh"
#include "util.hh"

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
constexpr int max_nba_syms = 1 << (8 * sizeof(sym_t));

// type for "big" automata
using state_t = uint32_t;

// type of acceptance sets
using acc_t = unsigned;

struct SWAMeta {
  // automaton name
  string name;
  // atomic props
  vector<string> aps;
};

enum Acceptance { BUCHI, PARITY };

// acceptance-labelled KS with tagged nodes, constant size alphabet
template <Acceptance A,typename T,class tag_storage = generic_trie_bimap<T, small_state_t, state_t>>
struct SWA {
  typedef unique_ptr<SWA<A, T, tag_storage>> uptr;
  typedef shared_ptr<SWA<A, T, tag_storage>> sptr;
  typedef unique_ptr<bimap<T, state_t, tag_storage>> tag_ptr;

  SWA() {}

  SWA(vector<string> const& aps, state_t const& start) : init(start) {
      add_state(start);
      meta.aps = aps;
    }

  SWAMeta meta;

  state_t init;  // initial state
  set<acc_t> accsets; //all used acceptance sets

  // state acceptance marks
  map<state_t, vector<acc_t>> acc;

  // list of successors per symbol
  map<state_t, map<sym_t, vector<state_t>>> adj;
  // shared_ptr<bimap<state_t,T>> tag;
  tag_ptr tag;

  //TODO: state iterator wrapping map iterator?

  vector<string> const& aps() const { return meta.aps; }
  inline size_t num_syms() const { return 1 << meta.aps.size(); }

  size_t num_states() const { return adj.size(); }
  vector<state_t> states() const { return map_get_keys(adj); }
  bool has_state(state_t const& s) const { return map_has_key(adj, s); }

  // add a state if not yet existing. return true if added
  bool add_state(state_t const& s) {
    if (!has_state(s)) {
      adj[s] = {};
      return true;
    }
    return false;
  }

  // bool has_acc(state_t const& s, acc_t a) const { return acc.at(s).find(a) != acc.at(s).end(); }
  // bool add_acc(state_t const& s, acc_t a) const { return acc.at(s).find(a) != acc.at(s).end(); }
  // bool remove_acc(state_t const& s, acc_t a) const { return acc.at(s).find(a) != acc.at(s).end(); }

  bool has_accs(state_t const& s) const { return map_has_key(acc, s); }
  vector<acc_t> const& get_accs(state_t const& s) const { return acc.at(s); }
  void set_accs(state_t const& s, vector<acc_t> const& a) {
    if (!has_state(s))
      throw runtime_error("ERROR: Adding acceptance to non-exinting state!");

    acc[s] = a;
    copy(begin(a),end(a), inserter(accsets, end(accsets)));
  }

  vector<sym_t> outsyms(state_t const& p) const {
    if (!has_state(p))
      throw runtime_error("ERROR: callung outsyms on non-exinting state!");

    auto const& edges = adj.at(p);
    return map_get_keys(edges);
  }


  bool state_has_outsym(state_t const& s, sym_t const& x) const {
    auto const& edges = adj.at(s);
    return map_has_key(edges, x);
  }

  vector<state_t> succ(state_t const& p, sym_t const& x) const {
    if (!has_state(p)) return {};
    auto const& edges = adj.at(p);
    if (!map_has_key(edges, x)) return {};
    return edges.at(x);
  }

  vector<state_t> const& succ_raw(state_t const& p, sym_t const& x) const { return adj.at(p).at(x); }

  vector<state_t> succ(state_t const& p) const {
    if (!has_state(p)) return {};
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
      if (map_has_key(adj,it))
        adj.erase(adj.find(it));  // kill state + edges
    }
    for (auto const& it : tokill) {
      if (map_has_key(acc,it))
        acc.erase(acc.find(it));  // kill acc. mark
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
  void merge_states(SWA<A, T>& aut, vector<state_t>& eqclass) {
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

template <Acceptance A, typename T, typename S>
vector<small_state_t> powersucc(SWA<A, T> const& ks, std::vector<S> const& ps,
                                sym_t const& x) {
  std::set<state_t> suc;
  for (auto const& p : ps) {
    if (!ks.has_state(p))
      continue;
    if (!ks.state_has_outsym(p,x))
      continue;
    auto const& qs = ks.succ(p, x);
    std::copy(qs.begin(), qs.end(), std::inserter(suc, suc.end()));
  }
  return vector<S>(suc.begin(), suc.end());
}

using BA = SWA<BUCHI, Unit>;

}  // namespace nbautils
