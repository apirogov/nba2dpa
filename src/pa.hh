#include <vector>
#include <cassert>
#include <functional>
#include <unordered_map>
#include "swa.hh"
#include "io.hh"
#include "common/util.hh"
#include "common/acceptance.hh"
#include "common/algo.hh"

namespace nbautils {
using namespace std;
using namespace nbautils;

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
  PAProdState succ(state_t l, int pl, state_t r, int pr, bool intersect) const;

  string to_string() const;
  bool operator<(PAProdState const& other) const;
  bool operator==(PAProdState const& other) const;
};
}

namespace std {
using namespace nbautils;
template <>
    struct hash<PAProdState>
    {
        size_t operator()(PAProdState const& k) const
        {
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
}

namespace nbautils {
using namespace std;
using namespace nbautils;

//return empty PA accepting nothing
template<typename T,template <typename... Args> class S>
typename SWA<T,S>::uptr empty_pa(vector<string> const& aps) {
  auto pa = std::make_unique<SWA<T,S>>(Acceptance::PARITY, "Empty", aps);
  pa->set_patype(PAType::MIN_EVEN);
  pa->add_state(0);
  pa->set_init({0});
  pa->set_accs(0, {1});
  for (auto i=0; i<pa->num_syms(); i++) {
    pa->set_succs(0, i, {0});
  }
  return move(pa);
}

//returns accepting reachable (sub)scc from which a run can be easily constructed
template<typename T, template <typename... Args> class S>
vector<state_t> find_acc_pa_scc(SWA<T,S> const& aut) {
  assert(aut.acond == Acceptance::PARITY);
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

//TODO: PA emptiness
//TODO: lang. containment check
using PAP = SWA<PAProdState,naive_unordered_bimap>;
template<typename A, typename B, template <typename... Args1> class SA, template <typename... Args2> class SB>
PAP::uptr pa_prod(SWA<A, SA> const& aut_a, SWA<B, SB> const& aut_b, bool intersect) {
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

  int l    = aut_a.get_init().front();
  // int pl   = aut_a.get_accs(l).front();
  int lmin = aut_a.get_accsets().front();
  int lmax = aut_a.get_accsets().back();

  int r    = aut_b.get_init().front();
  // int pr   = aut_b.get_accs(r).front();
  int rmin = aut_b.get_accsets().front();
  int rmax = aut_b.get_accsets().back();
  PAProdState inittag(l, lmin, lmax, r, rmin, rmax);
  // inittag = inittag.succ(l, pl, r, pr, intersect);
  //TODO: the initial state is often in a trivial SCC and unneeded!
  //
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
                                      sucb, aut_b.get_accs(sucb).front(), intersect);

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

  return move(pa);
}

//TODO: something is not right. either intersection or this is broken
template<typename T, template <typename... Args> class S>
vector<vector<state_t>> pa_equiv_states(SWA<T,S> const& aut) {
  assert(aut.acond == Acceptance::PARITY);

  //TODO: check each color separately, not each with each?

  auto const sts = aut.states();
  map<state_t, state_t> eq;
  for (auto const st : sts)
    eq[st] = st;

  for (int i=0; i<(int)sts.size(); ++i) {
    for (int j=i+1; j<(int)sts.size(); ++j) {
      // cout << sts[i] << "," << sts[j] << endl;
      if (aut.get_accs(sts[i]) != aut.get_accs(sts[j]))
        continue; //different priorities can not be equivalent! TODO?
      if (eq.at(sts[j]) != sts[j])
        continue; //this one is already equiv to someone

      auto aut_a(aut);
      auto aut_b(aut);
      aut_a.set_init({sts[i]});
      aut_b.set_init({sts[j]});
      complement_pa(aut_b);

      auto a_times_not_b = pa_prod(aut_a, aut_b, true);

      // cout << "a x not b:" << endl;
      // print_hoa(*a_times_not_b);
      // vector<state_t> pref;
      // vector<state_t> cyc;
      // tie(pref, cyc) = get_acc_pa_run(*a_times_not_b);
      // cout << seq_to_str(pref) << " | " << seq_to_str(cyc) << endl;
      // auto scc = find_acc_pa_scc(*a_times_not_b);
      // cout << seq_to_str(scc) << endl;

      if (!pa_is_empty(*a_times_not_b))
        continue;

      complement_pa(aut_a);
      complement_pa(aut_b);
      auto not_a_times_b = pa_prod(aut_a, aut_b, true);

      // cout << "not a x b:" << endl;
      // print_hoa(*not_a_times_b);

      if (!pa_is_empty(*not_a_times_b))
        continue;

      // we're here so both are empty so states are equivalent
      // cout << sts[i] << " ~ " << sts[j] << endl;
      eq[sts[j]] = eq[sts[i]];
    }
  }
  return group_by(sts, [&](state_t a, state_t b){ return eq.at(a)==eq.at(b); });
}

//interpret a BA as min even PA by setting missing priorities to 1
template<typename T, template <typename... Args> class S>
void ba_to_pa(SWA<T, S> &aut) {
  assert(aut.acond == Acceptance::BUCHI);
  for (auto const s : aut.states()) {
    if (!aut.has_accs(s))
      aut.set_accs(s, {1});
  }

  aut.acond = Acceptance::PARITY;
  aut.set_patype(PAType::MIN_EVEN);
}

//map over the priorities (using a function obtained with other functions provided here)
template<typename T, template <typename... Args> class S>
void transform_priorities(SWA<T,S> &aut, function<acc_t(acc_t)> const& pf) {
  assert(is_colored(aut));
  assert(aut.acond == Acceptance::PARITY);

  for (auto const s : aut.states())
    aut.set_accs(s, {pf(aut.get_accs(s).front())});
}

//complement PA by flipping parity of states
template<typename T, template <typename... Args> class S>
void complement_pa(SWA<T,S> &aut) {
  transform_priorities(aut, [](int p){return p+1;});
}

//switch between different parity conditions by changing priorities
template<typename T, template <typename... Args> class S>
void change_patype(SWA<T,S> &aut, PAType pt) {
  assert(aut.acond == Acceptance::PARITY);

  auto const pris = aut.get_accsets();
  auto f = priority_transformer(aut.get_patype(), pt, pris.front(), pris.back());
  transform_priorities(aut, f);
  aut.set_patype(pt);
}

template<typename T, template <typename... Args> class S>
auto minimize_priorities(SWA<T,S>& aut) {
  assert(is_colored(aut));
  assert(aut.acond == Acceptance::PARITY);

  auto const orig_patype = aut.get_patype();

  function<vector<state_t>(state_t)> const sucs = [&](state_t v){ return aut.succ(v); };
  function<int(state_t)> const get_pri = [&](state_t v){return aut.get_accs(v).front();};

  auto const oldpris = aut.get_accsets();
  auto const to_max_odd = priority_transformer(aut.get_patype(), PAType::MAX_ODD,
                                               oldpris.front(), oldpris.back());

  function<int(state_t)> max_odd_pri = [&](state_t v){ return to_max_odd(get_pri(v)); };
  auto const primap = pa_minimize_priorities(aut.states(), sucs, max_odd_pri);

  for (auto it : primap) {
    aut.set_accs(it.first, {(acc_t)it.second});
    // cout << it.first << " : " << it.second << endl;
  }
  aut.set_patype(PAType::MAX_ODD);
  change_patype(aut, orig_patype); //transform back

  return primap;
}

// https://en.wikipedia.org/wiki/DFA_minimization
//
// hopcroft algorithm to calculate equivalence classes that can be merged
// without changing the language / output behaviour
// requires complete, deterministic automaton
template<typename T, template <typename... Args> class S>
vector<vector<state_t>> get_equiv_states(SWA<T,S> const& aut) {
  auto sorted = aut.states();
  unordered_map<state_t, int> clr;
  for (auto s : sorted)
    clr[s] = aut.get_accs(s).front();

  unordered_map<state_t, vector<state_t>> suc;
  for (auto s : sorted)
    for (auto i=0; i<aut.num_syms(); i++) {
      suc[s].push_back(aut.succ(s,i).front());
    }

  sort(begin(sorted), end(sorted), [&](auto a, auto b){ return clr.at(a) < clr.at(b); });
  auto startsets = group_by(sorted, [&](auto a, auto b){ return clr.at(a) == clr.at(b); });

  PartitionRefiner<state_t> p(startsets);

  // set seems not to make a difference :/
  auto const wvec=p.get_set_ids();
  auto const sym_set_cmp = [](auto const& a, auto const& b) {
        return lexicographical_compare(a->first, a->second, b->first, b->second);
      };
  auto w=set<PartitionRefiner<state_t>::sym_set, decltype(sym_set_cmp)>(cbegin(wvec), cend(wvec), sym_set_cmp);

  // auto w=p.get_set_ids();

  while (!w.empty()) {
    // auto const a = w.back(); w.pop_back();
    auto const a = *w.begin(); w.erase(w.begin());

    auto const sepset = p.get_elements_of(a); //need to take a copy, as it is modified in loop

    for (auto i=0; i<aut.num_syms(); i++) {
      auto const succ_in_a = [&](state_t st){ return sorted_contains(sepset, suc[st][i]); };

      for (auto& y : p.get_set_ids()) {
        //TODO: try to precalculate sep X subset of Y set instead of using predicate?
        auto const z = p.separate(y, succ_in_a);

        if (z) { //separation happened
          //TODO: w is vector, this is slow?
          if (contains(w, y)) //if y is in w, its symbol still is, and we need the other set
            // w.push_back(*z);
            w.emplace(*z);
          else { //if y is not in w, take smaller part as separator
            if (p.get_set_size(y) <= p.get_set_size(*z)) {
              // w.push_back(y);
              w.emplace(y);
            } else {
              // w.push_back(*z);
              w.emplace(*z);
            }
          }
        }
      }
    }
  }

  return p.get_refined_sets();
}

template<typename T, template <typename... Args> class S>
bool minimize_pa(SWA<T,S>& pa) {
  assert(pa.get_init().size()==1);

  // function<acc_t(state_t)> colors = [&](auto s){return pa.get_accs(s).front();};
  // function<state_t(state_t,sym_t)> xsucc = [&](auto p, auto s){return pa.succ(p,s).front();};
  // auto equiv = dfa_equivalent_states(pa.states(), colors, pa.num_syms(), xsucc);

  auto equiv = get_equiv_states(pa);
  pa.quotient(equiv);
  pa.normalize();

  return true;
}


}
