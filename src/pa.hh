#pragma once

#include <vector>
#include <cassert>
#include <functional>
#include <unordered_map>
#include "aut.hh"
#include "graph.hh"
#include "common/scc.hh"
#include "common/util.hh"
#include "common/parity.hh"
#include "common/part_refinement.hh"

#include <spdlog/spdlog.h>

namespace nbautils {
using namespace std;
using namespace nbautils;

//map over the priorities (using a function obtained with other functions provided here)
template<typename T>
void transform_priorities(Aut<T> &aut, function<pri_t(pri_t)> const& pf) {
  assert(aut.is_colored());
  for (auto const p : aut.states())
    for (auto const x : aut.state_outsyms(p))
      for (auto& es : aut.succ_edges(p,x))
        aut.mod_edge(p, x, es.first,  pf(es.second));
}

//complement PA by flipping parity of states
template<typename T>
void complement_pa(Aut<T> &aut) {
  transform_priorities(aut, [](int p){return p+1;});
}

//switch between different parity conditions by changing priorities
template<typename T>
void change_patype(Aut<T> &aut, PAType pt) {
  auto const f = priority_transformer(aut.get_patype(), pt, aut.pri_bounds());
  transform_priorities(aut, f);
  aut.set_patype(pt);
}

// ----------------------------------------------------------------------------

using EdgeNode = tuple<state_t, sym_t, state_t, pri_t>;

struct PAProdState {
  state_t a;
  state_t b;

  //representation of <bool = original automaton (false = left, true = right), important bad, less important good neighbor = int>
  vector<pair<bool,int>> priord;


  // construct an initial parity product state from initial state and parity bounds
  PAProdState(state_t l, int lmin, int lmax, state_t r, int rmin, int rmax);
  PAProdState(PAProdState const& other);
  PAProdState();

  // get new state from current with given new component states (and their prio) and adapted priorities
  pair<PAProdState, int> succ(state_t l, int pl, state_t r, int pr, bool fulldown=true) const;

  string to_string() const;
  bool operator<(PAProdState const& other) const;
  bool operator==(PAProdState const& other) const;
};
}

namespace std {
using namespace nbautils;

template <>
struct hash<PAProdState> {
  size_t operator()(PAProdState const& k) const {
    // Compute individual hash values for first, second and third
    // http://stackoverflow.com/a/1646913/126995
    size_t res = 17;
    res = res * 31 + hash<state_t>()( k.a );
    res = res * 31 + hash<state_t>()( k.b );
    for (auto const& it : k.priord) {
      res = res * 31 + (hash<int>()( it.second << it.first ));
    }
    return res;
  }
};

template <>
struct hash<EdgeNode> {
  size_t operator()(EdgeNode const& k) const {
    // Compute individual hash values for first, second and third
    // http://stackoverflow.com/a/1646913/126995
    size_t res = 17;
    res = res * 31 + hash<state_t>()( get<0>(k) );
    res = res * 31 + hash<sym_t>()( get<1>(k) );
    res = res * 31 + hash<state_t>()( get<2>(k) );
    res = res * 31 + hash<pri_t>()( get<3>(k) );
    return res;
  }
};

template <>
struct hash<pair<state_t, state_t>> {
  size_t operator()(pair<state_t,state_t> const& k) const {
    // Compute individual hash values for first, second and third
    // http://stackoverflow.com/a/1646913/126995
    size_t res = 17;
    res = res * 31 + hash<state_t>()( get<0>(k) );
    res = res * 31 + hash<sym_t>()( get<1>(k) );
    return res;
  }
};

}

namespace nbautils {
using namespace std;
using namespace nbautils;

//return empty PA accepting nothing
template<typename T>
Aut<T> new_empty_pa(vector<string> const& aps) {
  Aut<T> pa(true, "Empty", aps, 0);
  pa.set_patype(PAType::MIN_EVEN);
  pa.set_pri(0, 1);
  for (auto const i : pa.syms()) {
    pa.add_edge(0, i, 0, -1);
  }
  return pa;
}

//this is the naive and straight-forward emptiness check in DTPA
//assumes that all states are reachable
//returns accepting (sub)scc from which a run can be easily constructed
//returns an edge with a good priority to build a run from
template<typename T, typename F>
pair<vector<state_t>,EdgeNode> find_acc_pa_scc_ext(Aut<T> const& aut, F extra_predicate) {
  assert(!aut.is_sba());

  auto const stronger = stronger_op_f(aut.get_patype());
  auto const stronger_of = stronger_priority_f(aut.get_patype());
  auto const pris = aut.pri_bounds();
  int const weakest = stronger(pris.first, pris.second) ? pris.second : pris.first;
  int const lowprio = pa_acc_is_max(aut.get_patype()) ? weakest-2 : weakest+2;

  // for each good priority
  for (auto const p : aut.pris()) {
    if (!good_priority(aut.get_patype(), p))
      continue;
    // cerr << "check " << p << endl;

    // try to find non-trivial SCC that never needs to visit any stronger priority
    auto const restricted_succ = [&](auto const& s){
      vector<state_t> sucs;
      // cerr << "suc " << s << endl;
      for (auto const x : aut.state_outsyms(s))
        for (auto const it : aut.succ_edges(s,x))
          if (!stronger(it.second, p))
            sucs.push_back(it.first);
      sucs |= ranges::action::sort | ranges::action::unique;
      return sucs;
    };

    // get SCCs s.t. stronger priorities are forbidden
    auto const scci = get_sccs(aut.states(), restricted_succ);
    auto const triv = trivial_sccs(restricted_succ, scci);
    if (scci.sccs.size() == triv.size())
      continue;

    // return any non-trivial good subscc where
    // an acc cycle can be found with current priority
    for (auto const& it : scci.sccs) {
      // cerr << seq_to_str(it.second) << endl;
      if (contains(triv, it.first))
        continue;

      int sccprio = lowprio;
      EdgeNode en;
      for (auto const s : it.second)
        for (auto const x : aut.state_outsyms(s))
          for (auto const e : aut.succ_edges(s,x)) {
            if (scci.scc_of.at(e.first) != scci.scc_of.at(s))
              continue;
            if (stronger(e.second, p))
              continue;
            int newprio = stronger_of(sccprio, e.second);
            if (newprio != sccprio) {
              en = make_tuple(s, x, e.first, e.second);
            }
            sccprio = newprio;
          }
      if (sccprio == p) {
        // cerr << "nonempty: " << seq_to_str(it.second) << endl;
        //
        if (extra_predicate(it.second)) // <- for BA/DPA inclusion check
          return make_pair(it.second, en);
      }
    }
  }
  return {};
}

template<typename T>
pair<vector<state_t>,EdgeNode> find_acc_pa_scc(Aut<T> const& aut) {
  return find_acc_pa_scc_ext(aut, const_true);
}

//empty := no accepting subscc
template<typename T>
bool pa_is_empty(Aut<T> const& aut) {
  return find_acc_pa_scc(aut).first.empty();
}

//acc run := prefix into accepting subscc + cycle in accepting subscc
template<typename T>
pair<vector<state_t>,vector<state_t>> get_acc_pa_run(Aut<T> const& aut) {
  auto ret = find_acc_pa_scc(aut);
  auto& scc = ret.first;
  if (scc.empty())
    return {};
  EdgeNode const& edge = ret.second;

  // cout << "getting pref" << endl;
  vector<state_t> const fin = find_path_from_to(aut, aut.get_init(), get<0>(edge));
  // cout << "getting cyc" << endl;
  vector<state_t> const cyc = find_path_from_to(aut, get<2>(edge), get<0>(edge));
  return make_pair(move(fin),move(cyc));
}

//TODO: get run, reconstruct word by choosing edge to known successor on run
// pair<vector<sym_t>,vector<sym_t>> get_acc_pa_word(SWA<T,S> const& aut) {
// }

// ----------------------------------------------------------------------------

using PAP = Aut<PAProdState>;

//TODO: same topo stuff as with determinization, using "raw" product as base
template<typename A, typename B>
PAP pa_union(Aut<A> const& aut_a, Aut<B> const& aut_b) {
  assert(aut_a.get_patype() == PAType::MIN_EVEN);
  assert(aut_b.get_patype() == PAType::MIN_EVEN);
  assert(aut_a.get_aps() == aut_b.get_aps());
  assert(is_colored(aut_a));
  assert(is_colored(aut_b));
  assert(is_deterministic(aut_a));
  assert(is_deterministic(aut_b));
  assert(!aut_a.is_sba());
  assert(!aut_b.is_sba());

  state_t const myinit = 0;
  auto pa = Aut<PAProdState>(false, "PA Product (unnamed)", aut_a.get_aps(), myinit);
  pa.set_patype(PAType::MIN_EVEN);
  pa.tag_to_str = [](ostream& out, auto const& t){ out << t.to_string(); };

  state_t const l    = aut_a.get_init();
  int const lmin = aut_a.pri_bounds().first;
  int const lmax = aut_a.pri_bounds().second;

  state_t const r    = aut_b.get_init();
  int const rmin = aut_b.pri_bounds().first;
  int const rmax = aut_b.pri_bounds().second;

  PAProdState const inittag(l, lmin, lmax, r, rmin, rmax);
  pa.tag.put(inittag, myinit);

  // int numvis=0;
  bfs(myinit, [&](auto const& st, auto const& visit, auto const&) {
    // get inner states of current state
    auto const curst = pa.tag.geti(st);

    // ++numvis;
    // if (numvis % 100 == 0) //progress indicator
    //   cerr << numvis << endl;

    vector<sym_t> syms;
    ranges::set_intersection( aut_a.state_outsyms(curst.a), aut_b.state_outsyms(curst.b)
                            , ranges::back_inserter(syms) );

    for (auto const i : syms) {
      for (auto const ea : aut_a.succ_edges(curst.a,i))
      for (auto const eb : aut_b.succ_edges(curst.b,i)) {
        auto const tmp = curst.succ(ea.first, ea.second,
                                    eb.first, eb.second);
        auto const& sucprod = tmp.first;
        auto const prio = tmp.second;

        //check whether there is already a state in the graph with this label
        auto const sucst = pa.tag.put_or_get(sucprod, pa.num_states());

        //if this is a new successor, add it to graph and enqueue it:
        if (!pa.has_state(sucst))
          pa.add_state(sucst);
        // create edge
        pa.add_edge(st, i, sucst, prio);
        // schedule for bfs
        visit(sucst);
      }
    }
  });

  return move(pa);
}

//check language inclusion by emptiness test of corresp. intersection via union
template<typename A, typename B>
bool dpa_inclusion(Aut<A> const& a, Aut<B> const& b) {
  auto aut_a(a);
  auto aut_b(b);

  aut_a.make_colored();
  complement_pa(aut_a);
  if (aut_a.is_sba())
    aut_a.to_tba();

  aut_b.make_colored();
  if (aut_b.is_sba())
    aut_b.to_tba();

  auto prodpa = pa_union(aut_a, aut_b);
  complement_pa(prodpa);
  // print_aut(prodpa);

  return pa_is_empty(prodpa);
}

template<typename A, typename B>
bool dpa_equivalence(Aut<A> const& aut1, Aut<B> const& aut2) {
  return dpa_inclusion(aut1, aut2) && dpa_inclusion(aut2, aut1);
}

// ----------------------------------------------------------------------------

//product of SBA and TDPA. edge priorities = incremented prios of DPA edges
template<typename A, typename B>
Aut<pair<state_t, state_t>> ba_dpacomp_prod(Aut<A> const& ba, Aut<B> const& dpa) {
  assert(ba.get_aps() == dpa.get_aps());
  assert(is_deterministic(dpa));
  assert(ba.is_sba());
  assert(!dpa.is_sba());

  state_t const myinit = 0;
  auto pa = Aut<pair<state_t, state_t>>(false, "BA/DPAcomp product (unnamed)", dpa.get_aps(), myinit);
  pa.set_patype(PAType::MIN_EVEN);
  pa.tag_to_str = [](ostream& out, auto const& t){ out << t.first << "," << t.second; };
  pa.tag.put(make_pair(ba.get_init(),dpa.get_init()), myinit);

  // int numvis=0;
  bfs(myinit, [&](auto const& st, auto const& visit, auto const&) {
    auto const curst = pa.tag.geti(st);

    vector<sym_t> syms;
    ranges::set_intersection( ba.state_outsyms(curst.first), dpa.state_outsyms(curst.second)
                            , ranges::back_inserter(syms) );

    for (auto const i : syms) {
      for (auto const ea : ba.succ_edges(curst.first,i))
      for (auto const eb : dpa.succ_edges(curst.second,i)) {
        auto const sucprod = make_pair(ea.first, eb.first);
        auto const prio = eb.second + 1; //inc for complement

        //check whether there is already a state in the graph with this label
        auto const sucst = pa.tag.put_or_get(sucprod, pa.num_states());

        //if this is a new successor, add it to graph and enqueue it:
        if (!pa.has_state(sucst))
          pa.add_state(sucst);
        pa.add_edge(st, i, sucst, prio); // create edge
        visit(sucst); // schedule for bfs
      }
    }
  });
  return move(pa);
}

//check that given SNBA is included (language-wise) in TDPA
//this is a custom emptiness check in a mixed product
//useful when the DPA is obtained by iterated underapproximation
template<typename A, typename B>
bool ba_dpa_inclusion(Aut<A> const& ba, Aut<B> const& dpa) {
  auto const ppa = ba_dpacomp_prod(ba, dpa);
  return find_acc_pa_scc_ext(ppa, [&](vector<state_t> const& scc) {
      for (auto const st : scc) {
        auto const pst = ppa.tag.geti(st);
        if (ba.state_buchi_accepting(pst.first))
          return true;
      }
      return false;
    }).first.empty();
}

// ----------------------------------------------------------------------------

template <typename Range, typename SuccFun>
int max_chain(unordered_map<EdgeNode,int>& newpri, Range const& p, SuccFun const& get_succs, int curmax) {
  if (p.empty()) //by definition, empty set has no chain
    return 0;

  // edges we're allowed to consider with given forbidden priority bound
  auto const allowed_edge = [&curmax](EdgeNode const& e){return get<3>(e) < curmax; };
  // edges we're allowed to consider due to given state set
  auto const target_in_p = [&](EdgeNode const& e){ return ranges::find(p, get<2>(e))!=ranges::end(p); };
  // returns only allowed successors for the SCC search
  auto const restricted_succs = [&](state_t v) {
    // important: requires that successors pre-sorted correctly!
    return get_succs(v)
         | ranges::view::take_while(allowed_edge)
         | ranges::view::filter(target_in_p)
         | ranges::view::transform([](EdgeNode const& e){ return get<2>(e);})
         ;
  };

  auto const scci = get_sccs(p, restricted_succs);
  auto const triv = trivial_sccs(restricted_succs, scci);

  int maxlen = 0;
  for (auto const& it : scci.sccs) {
    auto const& scc = it.second;

    if (contains(triv, it.first)) //skip trivial SCCs
      continue;

    // get maximal outgoing priority in SCC
    auto const target_in_scc = [&](EdgeNode const& e){return contains(scc,get<2>(e));};
    auto const max_pri_edge = [&](state_t s){
      auto const tmp = get_succs(s);
      return tmp.empty() ? 0 : ranges::max(tmp | ranges::view::take_while(allowed_edge)
                             | ranges::view::filter(target_in_scc)
                             | ranges::view::transform([](EdgeNode const& e){ return get<3>(e); }) );
    };
    int const scc_pri = ranges::max(scc | ranges::view::transform(max_pri_edge));

    int m = 0;
    if (scc_pri > 0) {
      m = max_chain(newpri, ranges::view::bounded(scc), get_succs, scc_pri);

      if ((scc_pri - m) % 2 == 1) //parity alternation -> requires new priority
        m++;
    }

    //get edges that don't belong to the "derivative" of current SCC, these can be assigned a priority now
    for (auto const s : scc) {
      auto nonderiv_edges = get_succs(s)
        | ranges::view::drop_while([&](EdgeNode const& e){return get<3>(e) < scc_pri;})
        | ranges::view::filter(target_in_scc);
      for (auto const& e : ranges::view::bounded(nonderiv_edges))
        newpri[e] = m;

      //assign unset derivative edges the new prio of current SCC
      //(as they apparently don't have cycles for any smaller restriction)
      //NOTE: can be left out. any priority from 0 up to current SCC prio is valid!
      //      may lead to different minimization results...
      auto unset_edges = get_succs(s)
        | ranges::view::take_while([&](EdgeNode const& e){return get<3>(e) < scc_pri;})
        | ranges::view::filter(target_in_scc);
      for (auto const& e : ranges::view::bounded(unset_edges))
        if (!map_has_key(newpri, e))
          newpri[e] = m;
    }

    maxlen = max(maxlen, m);
  }

  return maxlen;
}

// takes all states, outgoing edge (EdgeNode) function and maximal max-odd prio
// priorities on edges must be from a max odd acceptance (!!!)
// returns minimized edge -> priority mapping
// Paper: "Computing the Rabin Index of a parity automaton"
template <typename Range, typename SuccFun>
unordered_map<EdgeNode, int> pa_minimize_priorities(Range const& states, SuccFun const& get_succs, int maxoldpri) {
  unordered_map<EdgeNode, int> primap;
  max_chain(primap, states, get_succs, maxoldpri+1);
  return primap;
}

// take a DPA, minimize number of used priorities
template<typename T>
bool minimize_priorities(Aut<T>& aut, shared_ptr<spdlog::logger> log = nullptr) {
  assert(aut.is_colored());

  //priority function (more convenient to work with max odd here)
  auto const to_max_odd = priority_transformer(aut.get_patype(), PAType::MAX_ODD, aut.pri_bounds());

  if (log)
    log->info("preparing edge graph...");

  //calc list of outgoing edges, sorted by max-odd prio
  unordered_map<state_t,vector<EdgeNode>> esucs;
  bool has_edges = false;
  for (auto const p : aut.states()) {
    esucs[p] = {}; //init empty
    for (sym_t const& x : aut.state_outsyms(p))
      for (auto const& es : aut.succ_edges(p,x)) {
        // cerr << es.second << " mapped to " << to_max_odd(es.second) << endl;
        esucs[p].push_back(make_tuple(p, x, es.first, to_max_odd(es.second)));
        has_edges = true;
      }
    //sort successors by max odd prio (to hopefully speedup restriction)
    ranges::action::sort(esucs[p], [](auto const a, auto const b){ return get<3>(a) < get<3>(b); });
  }

  //catch edge case
  if (!has_edges) {
    if (log)
      log->info("Automaton has no edges! Nothing to do!");
    return true;
  }

  //corresponding successor-edge function
  auto const sucs = [&esucs](state_t const v){ return ranges::view::bounded(esucs.at(v)); };

  if (log)
    log->info("calculating new priorities...");

  //calculate priority map (old edge pri -> new edge pri)
  auto const strongest = pa_acc_is_min(aut.get_patype()) ? aut.pris().front() : aut.pris().back();
  auto const primap = pa_minimize_priorities(aut.states(), sucs, to_max_odd(strongest));

  if (log)
    log->info("applying new priority map...");

  //calculate conversion back from max odd to original
  auto const mmel = ranges::minmax(ranges::view::values(primap));
  auto const from_max_odd = priority_transformer(PAType::MAX_ODD, aut.get_patype(), mmel);
  //map over priorities, transforming obtained to original acc. type
  for (auto const& it : esucs) {
    for (auto const& e : it.second) {
      auto const max_odd_prio = map_has_key(primap, e) ? primap.at(e) : 0;
      auto const new_edge_prio = from_max_odd(max_odd_prio);
      // cerr << max_odd_prio << " unmapped to " << new_edge_prio << endl;
      aut.mod_edge(get<0>(e), get<1>(e), get<2>(e), new_edge_prio);
    }
  }

  return true;
}

// https://en.wikipedia.org/wiki/DFA_minimization
//
// hopcroft algorithm to calculate equivalence classes that can be merged
// without changing the language / output behaviour
// requires complete, deterministic automaton with colored edges
template<typename T>
vector<vector<state_t>> get_equiv_states(Aut<T> const& aut) {
  // assign each combination of output priorities per symbol a number
  // -> for efficiency, this is the "color" of each state
  map<state_t, int> clr;
  int i=0;
  map<vector<pri_t>, int> clrs;
  for (auto const s : aut.states()) {
    auto const cvec = aut.syms()
      | ranges::view::transform([s,&aut](sym_t x){ return cbegin(aut.succ_edges(s,x))->second; })
      | ranges::to_vector;
    if (!clrs[cvec])
      clrs[cvec] = ++i;
    clr[s] = clrs[cvec];
  }
  auto const color = [&](state_t const s){ return clr.at(s); };

  // obtain adj matrix for big speedup (per sym, succ of each state)
  vector<vector<state_t>> mat(aut.num_syms(),vector<state_t>(aut.num_states(), -1));
  for (state_t const p : aut.states()) {
    for (sym_t const x : aut.state_outsyms(p)) {
      for (auto const& es : aut.succ_edges(p,x)) {
        mat[x][p] = es.first;
      }
    }
  }

  // partition states by initial color (= behaviour profile)
  vector<state_t> states = aut.states();
  states |= ranges::action::sort([&color](auto a, auto b){ return color(a) < color(b); });
  vector<vector<state_t>> const startsets = ranges::view::group_by(states,
    [&color](auto a, auto b){ return color(a) == color(b); });

  PartitionRefiner<state_t> p(startsets);

  // set seems not to make a difference :/
  /*
  auto const wvec=p.get_set_ids();
  auto const sym_set_cmp = [](auto const& a, auto const& b) {
        // return lexicographical_compare(a->first, a->second, b->first, b->second);
        return &(a->first) < &(b->first);
      };
  auto w=set<PartitionRefiner<state_t>::sym_set, decltype(sym_set_cmp)>(cbegin(wvec), cend(wvec), sym_set_cmp);
  */
  auto w=p.get_set_ids();

  vector<state_t> sepset;
  sepset.reserve(aut.states().size());
  while (!w.empty()) {
    auto const a = w.back(); w.pop_back();
    // auto const a = *w.begin(); w.erase(w.begin());

    // auto const sepset = p.get_elements_of(a); //need to take a copy, as it is modified in loop
    p.get_elements_of(a, sepset); //need to take a copy, as it is modified in loop

    for (auto const i : aut.syms()) {
      // auto const succ_in_a = [&](state_t st){ return sorted_contains(sepset, cbegin(aut.succ_edges(st, i))->first); };
      auto const succ_in_a = [&](state_t st){ return sorted_contains(sepset, mat[i][st]); };

      for (auto& y : p.get_set_ids()) {
      // for (auto it = begin(p.get_sets()); it!=end(p.get_sets()); ++it) {
        // PartitionRefiner<state_t>::sym_set y = it;

        //TODO: try to precalculate sep X subset of Y set instead of using predicate?
        auto const z = p.separate(y, succ_in_a);

        if (z) { //separation happened
          // ++it; //need to skip new set

          //TODO: w is vector, this is slow?
          if (contains(w, y)) //if y is in w, its symbol still is, and we need the other set
            w.push_back(*z);
            // w.emplace(*z);
          else { //if y is not in w, take smaller part as separator
            if (p.get_set_size(y) <= p.get_set_size(*z)) {
              w.push_back(y);
              // w.emplace(y);
            } else {
              w.push_back(*z);
              // w.emplace(*z);
            }
          }
        }
      }
    }
  }

  return p.get_refined_sets();
}

template<typename T>
bool minimize_pa(Aut<T>& pa, shared_ptr<spdlog::logger> log = nullptr) {
  assert(pa.is_complete());
  assert(pa.is_colored());

  if (log)
    log->info("Calculating equivalent states...");
  auto const equiv = get_equiv_states(pa);
  // for (auto const& eq : equiv)
  //   cerr << seq_to_str(eq) << endl;
  // cerr << "--" << endl;

  if (log)
    log->info("Merging equivalent states...");
  pa.quotient(equiv);
  // cerr << "quotiented" << endl;

  auto const unreach = unreachable_states(pa, pa.get_init());
  pa.remove_states(unreach);
  // cerr << "removed " << unreach.size() << " unreachable" << endl;

  // remove rejecting sink (which is unique after min.) in complete automaton
  // unless it is the initial state (i.e. automaton with empty language)
  state_t rejsink = -1;
  for (auto const v : pa.states()) {
    //must be rej. self-loop
    bool rsink = true;
    for (auto const i : pa.state_outsyms(v)) {
      for (auto const es : pa.succ_edges(v, i)) {
        if (es.first != v || (es.second % 2) == 0)
          rsink = false;
      }
    }
    if (rsink) {
      rejsink = v;
      break;
    }
  }
  if (rejsink >= 0 && rejsink != pa.get_init())
    pa.remove_states({rejsink});

  pa.normalize();
  // cerr << "normalized" << endl;

  return true;
}


}
