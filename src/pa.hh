#include <vector>
#include <cassert>
#include <functional>
#include "swa.hh"
#include "common/util.hh"
#include "common/acceptance.hh"
#include "common/algo.hh"

namespace nbautils {
using namespace std;
using namespace nbautils;

//TODO: interpreting BA as PA and vice versa and checking whether BA \ PA = emptyset (BA <= PA)

//TODO
struct PAProdState {
  state_t a;
  state_t b;
  int prio;
  //representation of <bool = original automaton (false = left, true = right), important bad, less important good neighbor = int>
  vector<tuple<bool,int>> priord;

  bool operator<(PAProdState const& other) const;
};

//compare lexicographically
bool PAProdState::operator<(PAProdState const& other) const {
  if (a == other.a) {
    if (b == other.b) {
      if (prio == other.prio) {
        return priord < other.priord;
      }
      return prio < other.prio;
    }
    return b < other.b;
  }
  return a < other.a;
}

//map over the priorities (using a function obtained with other functions provided here)
template<typename T>
void transform_priorities(SWA<T> &aut, function<acc_t(acc_t)> const& pf) {
  assert(is_colored(aut));
  assert(aut.acond == Acceptance::PARITY);

  for (auto const s : aut.states())
    aut.set_accs(s, {pf(aut.get_accs(s).front())});
}

//complement PA by flipping parity of states
template<typename T>
void complement_pa(SWA<T> &aut, PAType pt) {
  transform_priorities(aut, [](int p){return p+1;});
}

//switch between different parity conditions by changing priorities
template<typename T>
void change_patype(SWA<T> &aut, PAType pt) {
  assert(aut.acond == Acceptance::PARITY);

  auto const pris = aut.get_accsets();
  auto f = priority_transformer(aut.get_patype(), pt, pris.front(), pris.back());
  transform_priorities(aut, f);
  aut.set_patype(pt);
}

template<typename T>
auto minimize_priorities(SWA<T>& aut) {
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
  }
  aut.set_patype(PAType::MAX_ODD);
  change_patype(aut, orig_patype); //transform back

  return primap;
}

template<typename T>
bool minimize_pa(SWA<T>& pa) {
  assert(pa.get_init().size()==1);

  function<acc_t(state_t)> colors = [&](auto s){return pa.get_accs(s).front();};
  function<state_t(state_t,sym_t)> xsucc = [&](auto p, auto s){return pa.succ(p,s).front();};

  auto equiv = dfa_equivalent_states(pa.states(), colors, pa.num_syms(), xsucc);

  auto initial = pa.get_init().front();
  bool seenini = false;
  for (auto ecl : equiv) {
    auto rep = ecl.back();
    if (!seenini) {
      auto it = lower_bound(begin(ecl), end(ecl), initial);
      if (it != end(ecl) && *it == initial) {
        ecl.erase(it);
        rep = initial;
        seenini = true;
      } else {
        ecl.pop_back();
      }
    } else {
      ecl.pop_back();
    }

    pa.merge_states(ecl, rep);
  }

  pa.normalize();

  return true;
}


}
