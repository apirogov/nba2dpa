#pragma once

#include <vector>
#include <cassert>
#include <functional>
#include <unordered_map>
#include "aut.hh"
#include "common/scc.hh"
#include "common/util.hh"
#include "common/parity.hh"
#include "common/part_refinement.hh"

#include <spdlog/spdlog.h>

namespace nbautils {
using namespace std;
using namespace nbautils;

using EdgeNode = tuple<state_t, sym_t, state_t, pri_t>;

struct PAProdState {
  state_t a;
  state_t b;
  int prio;

  //representation of <bool = original automaton (false = left, true = right), important bad, less important good neighbor = int>
  vector<pair<bool,int>> priord;


  // construct an initial parity product state from initial state and parity bounds
  PAProdState(state_t l, int lmin, int lmax, state_t r, int rmin, int rmax);
  PAProdState(PAProdState const& other);
  PAProdState();

  // get new state from current with given new component states (and their prio) and adapted priorities
  PAProdState succ(state_t l, int pl, state_t r, int pr, bool fulldown=true) const;

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
    res = res * 31 + hash<int>()( k.prio );
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

//returns accepting reachable (sub)scc from which a run can be easily constructed
/*
template<typename T>
vector<state_t> find_acc_pa_scc(Aut<T> const& aut) {
  assert(aut.get_patype() == PAType::MIN_EVEN);

  auto const st = aut.states();
  function<unsigned(state_t)> const get_pri = [&](state_t v){return aut.get_accs(v).front(); };

  for (auto const p : aut.get_accsets()) {
    if (p%2==1)
      continue;
    // cout << "check " << p << endl;

    // for each good priority, in descending importance order,
    // try to find non-trivial SCC that never needs to visit any stronger priority
    succ_fun<state_t> restricted_succ = [&](auto const& s){
      auto ret = aut.succ(s);
      ret.erase(remove_if(ret.begin(), ret.end(),
            [&aut,p](state_t q){ return aut.get_accs(q).front()<p; }), ret.end());
      return ret;
    };

    // get SCCs s.t. lower priorities are forbidden
    auto const scci = get_sccs(st, restricted_succ);
    auto const best_pri = fold_sccs<state_t,unsigned>(*scci, aut.get_accsets().back(), get_pri,
                            [](auto a, auto b){return min(a,b);});
    auto const triv = trivial_sccs(*scci, restricted_succ);
    if (scci->sccs.size() == triv.size())
      continue;

    // return any non-trivial good subscc where an acc cycle can be found
    // with current priority
    for (int i=0; i<(int)scci->sccs.size(); i++)
      if (!contains(triv, i) && best_pri.at(i)==p)
        return scci->sccs.at(i);
  }
  return {};
}

//empty := no accepting subscc
template<typename T, template <typename... Args> class S>
bool pa_is_empty(SWA<T,S> const& aut) {
  return find_acc_pa_scc(aut).empty();
}

//acc run := prefix into accepting subscc + cycle in accepting subscc
template<typename T, template <typename... Args> class S>
pair<vector<state_t>,vector<state_t>> get_acc_pa_run(SWA<T,S> const& aut) {
  auto scc = find_acc_pa_scc(aut);
  if (scc.empty())
    return {};

  sort(begin(scc), end(scc));
  auto const succ = swa_succ(aut);
  succ_fun<state_t> const restricted_succ = [&](auto const& p){
    return set_intersect(succ(p), scc);
  };

  //get a good important state as cycle pivot
  state_t pivot = scc.front();
  for (auto s : scc)
    if (aut.get_accs(pivot).front() > aut.get_accs(s).front() )
      pivot = s;

  // cout << "getting pref" << endl;
  auto const ret = find_path_from_to(aut.get_init().front(), pivot, succ);
  // cout << "getting cyc" << endl;
  auto const cyc = find_path_from_to(pivot, pivot, restricted_succ);
  return make_pair(move(ret),move(cyc));
}

//TODO: get run, reconstruct word by choosing edge to known successor on run
// pair<vector<sym_t>,vector<sym_t>> get_acc_pa_word(SWA<T,S> const& aut) {
// }

using PAP = SWA<PAProdState,naive_unordered_bimap>;
template<typename A, typename B, template <typename... Args1> class SA, template <typename... Args2> class SB>
PAP::uptr pa_union(SWA<A, SA> const& aut_a, SWA<B, SB> const& aut_b) {
  assert(aut_a.acond == Acceptance::PARITY);
  assert(aut_b.acond == Acceptance::PARITY);
  assert(aut_a.get_patype() == PAType::MIN_EVEN);
  assert(aut_b.get_patype() == PAType::MIN_EVEN);
  assert(aut_a.get_aps() == aut_b.get_aps());
  assert(is_deterministic(aut_a));
  assert(is_deterministic(aut_b));
  assert(is_colored(aut_a));
  assert(is_colored(aut_b));

  auto pa = std::make_unique<PAP>(Acceptance::PARITY, "", aut_a.get_aps());
  pa->set_patype(PAType::MIN_EVEN);
  pa->tag_to_str = [](PAProdState const& s){ return s.to_string(); };

  state_t l    = aut_a.get_init().front();
  int lmin = aut_a.get_accsets().front();
  int lmax = aut_a.get_accsets().back();

  state_t r    = aut_b.get_init().front();
  int rmin = aut_b.get_accsets().front();
  int rmax = aut_b.get_accsets().back();
  PAProdState inittag(l, lmin, lmax, r, rmin, rmax);

  state_t myinit = 0;
  pa->add_state(myinit);
  pa->set_init({myinit});
  pa->set_accs(myinit, {(unsigned)inittag.prio});
  pa->tag->put(inittag, myinit);

  int numvis=0;
  bfs(myinit, [&](auto const& st, auto const& visit, auto const&) {
    // get inner states of current state
    auto const curst = pa->tag->geti(st);

    ++numvis;
    if (numvis % 100 == 0) //progress indicator
      cerr << numvis << endl;

    for (auto i = 0; i < pa->num_syms(); i++) {
      auto const suca = aut_a.succ(curst.a, i).front();
      auto const sucb = aut_b.succ(curst.b, i).front();
      auto const sucprod = curst.succ(suca, aut_a.get_accs(suca).front(),
                                      sucb, aut_b.get_accs(sucb).front());

      //check whether there is already a state in the graph with this label
      auto const sucst = pa->tag->put_or_get(sucprod, pa->num_states());

      //if this is a new successor, add it to graph and enqueue it:
      if (!pa->has_state(sucst))
        pa->add_state(sucst);
      // assign priority according to resulting level
      if (!pa->has_accs(sucst))
        pa->set_accs(sucst,{(unsigned)sucprod.prio});
      // create edge
      pa->set_succs(st, i, {sucst});
      // schedule for bfs
      visit(sucst);
    }
  });
  //TODO: same topo stuff as with safra, using "raw" product

  return move(pa);
}

template<typename T, typename T2>
// vector<vector<state_t>> pa_equiv(SWA<T,S> const& aut1, SWA<T,S> const& aut2) {
bool pa_equivalent(Aut<T> const& aut1, Aut<T2> const& aut2) {
  auto aut_a(aut1);
  auto aut_b(aut2);
  complement_pa(aut_b);

  auto a_times_not_b = pa_union(aut_a, aut_b);
  complement_pa(*a_times_not_b);

  if (!pa_is_empty(*a_times_not_b))
    return false;

  complement_pa(aut_a);
  complement_pa(aut_b);
  auto not_a_times_b = pa_union(aut_a, aut_b);
  complement_pa(*not_a_times_b);

  if (!pa_is_empty(*not_a_times_b))
    return false;

  return true;
}
*/

// ----------------------------------------------------------------------------

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
