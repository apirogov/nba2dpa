#pragma once

#include "common/bimap.hh"
#include "common/util.hh"

#include <cassert>
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
using acc_t = uint32_t;

inline vector<small_state_t> to_small_state_t(vector<state_t> const& v) { return vector<small_state_t>(cbegin(v),cend(v)); }
// inline vector<small_state_t> to_big_state_t(vector<small_state_t> const& v) { return vector<state_t>(cbegin(v),cend(v)); }

enum class Acceptance { DYNAMIC, BUCHI, PARITY };

// acceptance-labelled KS with tagged nodes, fixed alphabet
template <Acceptance A,typename T,class tag_storage = naive_bimap<T, state_t>>
struct SWA {
  using uptr =  unique_ptr<SWA<A, T, tag_storage>>;
  using sptr =  shared_ptr<SWA<A, T, tag_storage>>;
  using tag_container =  bimap<T, state_t, tag_storage>;
  using tag_ptr = unique_ptr<tag_container>;

  private:
  bool normalized = true;
  string name;

  const Acceptance acond;   //immutable acceptance condition type

  vector<string> aps; //atomic props can be set just once

  // initial state(s)
  vector<state_t> init;
  //all used acceptance sets
  set<acc_t> accsets;
  // state acceptance marks
  map<state_t, vector<acc_t>> acc;

  // generic graph stuff

  public:
  string const& get_name() const { return name; }
  void set_name(string const& n) { name = n; }
  set<acc_t> const& get_accsets() const { return accsets; }

  // list of successors per symbol
  map<state_t, map<sym_t, vector<state_t>>> adj;
  // node tags
  tag_ptr tag;

  //TODO: state iterator wrapping map iterator? all sym iterator? succ sym iterator?

  SWA(string const& title="", vector<string> const& ap={}, vector<state_t> const& initial={})
    : name(title), acond(A), aps(ap), tag(make_unique<tag_storage>()) { set_init(initial); }

  void set_aps(vector<string> const& ap) {
    assert(aps.empty()); // Can set APs only once!
    aps = ap;
  }

  int num_syms() const { return 1 << aps.size(); }
  vector<string> const& get_aps() const { return aps; }

  size_t num_states() const { return adj.size(); }
  vector<state_t> states() const { return map_get_keys(adj); }
  bool has_state(state_t const& s) const { return map_has_key(adj, s); }

  // add a new state (must have unused id)
  void add_state(state_t const& s) {
    assert(!has_state(s));
    // if (!has_state(s)) {
      adj[s] = {};
      // return true;
    // }
    // return false;
  }

  void set_init(vector<state_t> const& initial) {
    assert(is_set_vec(initial));

    init = initial;
    //add states if they do not exist yet
    for (auto &v : init)
      if (!has_state(v))
        add_state(v);
  }
  vector<state_t> const& get_init() const { return init; }
  bool is_init(state_t const& s){ return contains(init, s); }

  // bool has_acc(state_t const& s, acc_t a) const { return acc.at(s).find(a) != acc.at(s).end(); }
  // bool add_acc(state_t const& s, acc_t a) const { return acc.at(s).find(a) != acc.at(s).end(); }
  // bool remove_acc(state_t const& s, acc_t a) const { return acc.at(s).find(a) != acc.at(s).end(); }

  bool has_accs(state_t const& s) const { return map_has_key(acc, s); }
  vector<acc_t> get_accs(state_t const& s) const { if (!has_accs(s)) return {}; else return acc.at(s); }
  void set_accs(state_t const& s, vector<acc_t> const& a) {
    assert(has_state(s)); // otherwise adding acceptance to non-exinting state
    assert(is_set_vec(a));

    acc[s] = a;
    copy(cbegin(a),cend(a), inserter(accsets, end(accsets)));
  }

  vector<sym_t> outsyms(state_t const& p) const {
    assert(has_state(p));
    return map_get_keys(adj.at(p));
  }

  bool state_has_outsym(state_t const& p, sym_t const& x) const {
    return map_has_key(adj.at(p), x);
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

  //renumber all states continuously starting from offset
  void normalize(state_t offset=0) {
    //calculate state renumbering
    map<state_t, state_t> m;
    state_t i=offset;
    for (auto const& it : adj) {
      m[it.first] = i++;
    }
    auto map_states = [&m](auto &s){
      if (!map_has_key(m,s))
        cerr << s << " not in state map!" << endl;
      return m.at(s);
    };

    map<state_t, vector<acc_t>> newacc;
    tag_ptr newtag = make_unique<tag_storage>();
    map<state_t, map<sym_t, vector<state_t>>> newadj;
    for (auto const& it : adj) {
      auto const olds = it.first;
      auto const news = m.at(olds);
      // cout << olds << " -> " << news << endl;

      // cout << "move acc" << endl;
      newacc[news] = acc[olds];

      // cout << "rewrite tag" << endl;
      if (tag->hasi(olds)) {
        // cout << " get tag" << endl;
        auto t = tag->geti(olds);
        // cout << " write tag" << endl;
        newtag->put(t, news);
      }

      // cout << "rewrite edges" << endl;
      auto tmp = it.second;
      for (auto& jt : tmp) {
        auto &v = jt.second;
        transform(begin(v), end(v), begin(v), map_states);
      }
      newadj[news] = tmp;
      // cout << "state done" << endl;

    }
    //replace old stuff
    transform(begin(init), end(init), begin(init), map_states);
    acc = newacc;
    tag = move(newtag);
    adj = newadj;
    normalized = true;
  }

  //stupidly paste another automaton (ignore initial status)
  void insert(SWA<A,T,tag_storage> const& other) {
    assert(set_intersect(states(), other.states()).empty()); // no common states

    for (auto const& it : other.adj) {
      auto st = it.first;
      // cout << it.second.size() << " ";
      adj[st] = it.second;
      // cout << adj.at(st).size() << endl;
      set_accs(st, other.get_accs(st));
      tag->put(other.tag->geti(st), st);
    }

    normalized = false;
  }

  // given set in sorted(!!) vector, kill states
  void remove_states(vector<state_t> const& tokill) {
    for (auto const& it : tokill) {
      if (map_has_key(adj,it))
        adj.erase(adj.find(it));  // kill state + edges
    }
    for (auto const& it : tokill) {
      if (map_has_key(acc,it))
        acc.erase(acc.find(it));  // kill acc. mark
    }
    for (auto const& it : tokill) {
      if (tag->hasi(it))
        tag->erase(it);  // kill tag
    }

    // kill states from edge targets
    for (auto& it : adj) {
      for (auto& jt : it.second) {
        auto res = set_diff(jt.second, tokill);
        jt.second = res;
      }
    }

    // kill from initial
    auto res = set_diff(init, tokill);
    init = res;

    normalized = false;
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

//common stuff working on any automaton

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

template <Acceptance A, typename T>
bool is_deterministic(SWA<A,T> const& aut) {
  for (auto const& p : aut.states()) {
    for (auto const& s : aut.outsyms(p))
      if (aut.succ(p,s).size() > 1) {
        cerr << (int) p << "-" << (int) s << " > ";
        for (auto q : aut.succ(p,s)) {
          cerr << (int) q << ",";
        }
        cerr << endl;

        return false;
      }
  }
  return true;
}

template <Acceptance A, typename T>
bool is_colored(SWA<A,T> const& aut) {
  for (auto const& p : aut.states()) {
    //TODO
  }
  return true;
}

}  // namespace nbautils
