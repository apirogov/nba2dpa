#pragma once

#include "common/bimap.hh"
#include "common/util.hh"

#include <boost/variant.hpp>

#include <cassert>
#include <functional>
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

inline vector<small_state_t> to_small_state_t(vector<state_t> const& v) {
  return vector<small_state_t>(cbegin(v),cend(v));
}
// inline vector<small_state_t> to_big_state_t(vector<small_state_t> const& v) { return vector<state_t>(cbegin(v),cend(v)); }


enum class Acceptance { DYNAMIC, BUCHI, PARITY };
enum class PAType { MIN_EVEN, MIN_ODD, MAX_EVEN, MAX_ODD };
inline constexpr PAType opposite_parity(PAType pt) {
  switch (pt) {
    case PAType::MIN_EVEN: return PAType::MIN_ODD;
    case PAType::MIN_ODD:  return PAType::MIN_EVEN;
    case PAType::MAX_EVEN: return PAType::MAX_ODD;
    case PAType::MAX_ODD:  return PAType::MAX_EVEN;
  }
}

inline constexpr PAType opposite_polarity(PAType pt) {
  switch (pt) {
    case PAType::MIN_EVEN: return PAType::MAX_EVEN;
    case PAType::MIN_ODD:  return PAType::MAX_ODD;
    case PAType::MAX_EVEN: return PAType::MIN_EVEN;
    case PAType::MAX_ODD:  return PAType::MIN_ODD;
  }
}

inline constexpr bool pa_acc_is_min(PAType a) {
  return a==PAType::MIN_EVEN || a==PAType::MIN_ODD;
}
inline constexpr bool pa_acc_is_even(PAType a) {
  return a==PAType::MIN_EVEN || a==PAType::MAX_EVEN;
}
inline constexpr bool pa_acc_is_max(PAType a) {
  return !pa_acc_is_min(a);
}
inline constexpr bool pa_acc_is_odd(PAType a) {
  return !pa_acc_is_even(a);
}

inline constexpr bool same_parity(PAType a, PAType b) {
  return (pa_acc_is_even(a) == pa_acc_is_even(b));
}
inline constexpr bool same_parity(acc_t a, acc_t b) {
  return (a%2==0)==(b%2==0);
}

inline constexpr bool same_polarity(PAType a, PAType b) {
  return (pa_acc_is_min(a) == pa_acc_is_min(b));
}

// acceptance-labelled KS with tagged nodes, fixed alphabet
template <Acceptance A,typename T,class tag_storage = naive_bimap<T, state_t>>
struct SWA {
  using uptr =  unique_ptr<SWA<A, T, tag_storage>>;
  using sptr =  shared_ptr<SWA<A, T, tag_storage>>;
  using tag_container =  bimap<T, state_t, tag_storage>;
  using tag_ptr = unique_ptr<tag_container>;
  using tag_printer = function<std::string(const T&)>;

  private:
  bool normalized = true;
  string name;

  const Acceptance acond;   //immutable acceptance condition type
  PAType patype;

  vector<string> aps; //atomic props can be set just once

  // initial state(s)
  vector<state_t> init;
  //all used acceptance sets
  map<acc_t,int> accsets;
  // state acceptance marks
  map<state_t, vector<acc_t>> acc;

  // generic graph stuff:

  // list of successors per symbol
  map<state_t, map<sym_t, vector<state_t>>> adj;

  public:
  string const& get_name() const { return name; }
  void set_name(string const& n) { name = n; }

  vector<acc_t> get_accsets() const { return map_get_keys(accsets); }
  PAType get_patype() const {
    assert(A==Acceptance::PARITY);
    return patype;
  }
  void set_patype(PAType t) {
    assert(A==Acceptance::PARITY);
    patype = t;
  }

  // node tags
  tag_ptr tag;
  tag_printer tag_to_str;

  //TODO: state iterator wrapping map iterator instead of states()?
  //all sym iterator for beauty? wrapper outsym iterator instead of outsyms()?

  SWA(string const& title="", vector<string> const& ap={}, vector<state_t> const& initial={})
    : name(title), acond(A), aps(ap), tag(make_unique<tag_storage>()),
      tag_to_str([](T const&){ return "<?>"; })
    {
      for (auto const& v : initial)
        add_state(v);
      set_init(initial);
    }

  string print_tag(state_t const& s) const {
    if (tag->hasi(s))
      return tag_to_str(tag->geti(s));
    return "";
  }

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

    if (s != num_states()) {
      normalized = false;
    }

    // if (!has_state(s)) {
    adj[s] = {};
      // return true;
    // }
    // return false;
  }

  void set_init(vector<state_t> const& initial) {
    assert(is_set_vec(initial));
    for (auto &v : initial)
      assert(has_state(v));

    init = initial;
  }
  vector<state_t> const& get_init() const { return init; }
  bool is_init(state_t const& s){ return contains(init, s); }

  // bool has_acc(state_t const& s, acc_t a) const { return acc.at(s).find(a) != acc.at(s).end(); }
  // bool add_acc(state_t const& s, acc_t a) const { return acc.at(s).find(a) != acc.at(s).end(); }
  // bool remove_acc(state_t const& s, acc_t a) const { return acc.at(s).find(a) != acc.at(s).end(); }

  bool has_accs(state_t const& s) const {
    assert(has_state(s));
    return map_has_key(acc, s);
  }

  vector<acc_t> get_accs(state_t const& s) const { if (!has_accs(s)) return {}; else return acc.at(s); }
  void set_accs(state_t const& s, vector<acc_t> const& a) {
    assert(has_state(s)); // otherwise adding acceptance to non-exinting state
    assert(is_set_vec(a));

    if (A==Acceptance::BUCHI)
      assert(a.empty() || a==vector<acc_t>{0});
    else if (A==Acceptance::PARITY)
      assert(a.size()<2);

    // remove count of old, if any
    if (has_accs(s))
      for (auto old : acc.at(s))
        if (map_has_key(accsets, old)) {
          accsets.at(old)--;
          if (accsets.at(old) == 0)
            accsets.erase(accsets.find(old));
        }

    // overwrite and increase counts of new, if any provided
    if (!a.empty()) {
      acc[s] = a;
      for (auto as : acc.at(s))
        accsets[as]++;
    } else { //remove if empty
      acc.erase(acc.find(s));
    }
  }

  // --------------------------------------------

  vector<sym_t> outsyms(state_t const& p) const {
    assert(has_state(p));
    return map_get_keys(adj.at(p));
  }


  bool state_has_outsym(state_t const& p, sym_t const& x) const {
    return map_has_key(adj.at(p), x);
  }

  void set_succs(state_t const& p, sym_t const& x, vector<state_t> const& v) {
    assert(has_state(p));
    assert(x < num_syms());
    assert(is_set_vec(v));
#ifndef NDEBUG
    for (auto q : v)
      assert(has_state(q));
#endif

    if (!v.empty())
      adj.at(p)[x] = v;
    else // setting empty vector remoes this outsym
      if (map_has_key(adj.at(p), x))
          adj.at(p).erase(adj.at(p).find(x));
  }

  vector<state_t> const& succ_raw(state_t const& p, sym_t const& x) const { return adj.at(p).at(x); }

  // x-label edge successors
  vector<state_t> succ(state_t const& p, sym_t const& x) const {
    // if (!has_state(p)) return {};
    assert(has_state(p));
    auto const& edges = adj.at(p);
    if (!map_has_key(edges, x)) return {};
    return edges.at(x);
  }

  // all successors (independent of symbol)
  vector<state_t> succ(state_t const& p) const {
    // if (!has_state(p)) return {};
    assert(has_state(p));
    auto const& edges = adj.at(p);
    // collect successors for any symbol
    set<state_t> s;
    for (auto const& it : edges) {
      copy(begin(it.second), end(it.second), inserter(s, end(s)));
    }
    return vector<state_t>(begin(s), end(s));
  }

  // --------------------------------------------

  // given set in sorted(!!) vector, kill states
  void remove_states(vector<state_t> const& tokill) {
    assert(is_set_vec(tokill));
    for (auto s : tokill)
      assert(has_state(s));

    for (auto const& it : tokill) {
      if (has_accs(it))
        set_accs(it, {});  // kill acc. mark
    }
    for (auto const& it : tokill) {
      if (tag->hasi(it))
        tag->erase(it);  // kill tag
    }
    for (auto const& it : tokill) {
      if (map_has_key(adj,it))
        adj.erase(adj.find(it));  // kill state + edges
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

  //stupidly paste another automaton (ignore initial status)
  void insert(SWA<A,T,tag_storage> const& other) {
    assert(set_intersect(states(), other.states()).empty()); // no common states
    assert(get_aps() == other.get_aps()); //same alphabet and acceptance
    assert(acond == other.acond); //same alphabet and acceptance

    for (auto const& it : other.adj) {
      auto st = it.first;
      // cout << it.second.size() << " ";
      adj[st] = it.second;
      // cout << adj.at(st).size() << endl;
      if (other.has_accs(st))
        set_accs(st, other.get_accs(st));
      if (other.tag->hasi(st))
        tag->put(other.tag->geti(st), st);
    }

    normalized = false;
  }

  // given set in sorted(!!) vector, merge states into one rep
  // TODO: what about initial states?
  void merge_states(vector<state_t> const& others, state_t rep) {
    assert(is_set_vec(others));

    if (others.empty()) return; //nothing to do

    assert(!contains(others,rep));
    assert(has_state(rep));
    for (auto q : others)
      assert(has_state(q));

    // add outgoing edges from all class members to representative
    for (auto const st : others) {
      for (auto const sym : outsyms(st)) {
        auto const sucs = succ(st, sym);
        auto const res = set_merge(sucs, succ(rep, sym));
        set_succs(rep, sym, res);
      }
    }

    // add ingoing edges into all class members to representative
    for (auto const st : states()) {
      for (auto const sym : outsyms(st)) {
        auto const sucs = succ(st, sym);
        auto const res = set_intersect(sucs, others);
        set_succs(st, sym, set_merge(sucs, {rep}));
      }
    }

    // finally kill the others!!
    remove_states(others);
  }

  //renumber all states continuously starting from offset
  //TODO: maybe don't do anything for id-mapped nodes? is this optimization realistic?
  map<state_t, state_t> normalize(state_t offset=0) {
    //calculate state renumbering
    map<state_t, state_t> m;
    state_t i=offset;
    bool needs_renumbering = false;
    for (auto const& it : adj) {
      m[it.first] = i++;
      if (m.at(it.first) != it.first)
        needs_renumbering = true;
    }
    if (!needs_renumbering) { //everything is fine already
      normalized = true;
      return {}; //empty map = did nothing
    }

    auto map_states = [&m](auto const& s){ return m.at(s); };

    map<state_t, vector<acc_t>> newacc;
    tag_ptr newtag = make_unique<tag_storage>();
    map<state_t, map<sym_t, vector<state_t>>> newadj;
    for (auto const& it : adj) {
      auto const olds = it.first;
      auto const news = m.at(olds);
      // cout << olds << " -> " << news << endl;

      // cout << "move acc" << endl;
      if (map_has_key(acc, olds))
        newacc[news] = acc.at(olds);

      // cout << "rewrite tag" << endl;
      if (tag->hasi(olds)) {
        // cout << " get tag" << endl;
        auto t = tag->geti(olds);
        // cout << " write tag" << t << " " << news << endl;
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
    return m;
  }

};

//common stuff working on any automaton

template <Acceptance A, typename T, typename S>
vector<small_state_t> powersucc(SWA<A, T> const& ks, std::vector<S> const& ps,
                                sym_t const& x, std::vector<small_state_t> const& asinks={}) {
  assert(is_set_vec(asinks));
  assert(is_set_vec(ps));

  std::set<small_state_t> suc;
  for (auto const& p : ps) {
    if (!ks.has_state(p))
      continue;
    if (!ks.state_has_outsym(p,x))
      continue;
    auto const& qs = to_small_state_t(ks.succ(p, x));

    //an edge leads into accepting sink state -> return accepting sink PS
    if (!asinks.empty() && !set_intersect(qs, asinks).empty())
      return vector<small_state_t>(cbegin(asinks), cend(asinks));

    std::copy(qs.cbegin(), qs.cend(), std::inserter(suc, suc.end()));
  }
  return vector<small_state_t>(suc.cbegin(), suc.cend());
}

template <Acceptance A, typename T>
bool is_deterministic(SWA<A,T> const& aut) {
  if (aut.get_init().size() > 1)
    return false;

  for (auto const& p : aut.states()) {
    for (auto const& s : aut.outsyms(p))
      if (aut.succ(p,s).size() > 1) {
        // cerr << (int) p << "-" << (int) s << " > ";
        // for (auto q : aut.succ(p,s)) {
        //   cerr << (int) q << ",";
        // }
        // cerr << endl;
        return false;
      }
  }
  return true;
}

// each state is part of exactly one accset
template <Acceptance A, typename T>
bool is_colored(SWA<A,T> const& aut) {
  for (auto const& p : aut.states()) {
    if (!aut.has_accs(p) || aut.get_accs(p).size() != 1)
      return false;
  }
  return true;
}

}  // namespace nbautils
