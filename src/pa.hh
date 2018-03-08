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

//switch between different parity conditions
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

  // for (auto p : aut.get_accsets()) {
  //   cout << p << " -> " << to_max_odd(p) << endl;
  // }

  function<int(state_t)> max_odd_pri = [&](state_t v){ return to_max_odd(get_pri(v)); };
  auto const primap = pa_minimize_priorities(aut.states(), sucs, max_odd_pri);

  // auto const minpris = map_get_vals(primap);
  // cout << seq_to_str(minpris) << endl;

  for (auto it : primap) {
    aut.set_accs(it.first, {(acc_t)it.second});
  }
  aut.set_patype(PAType::MAX_ODD);
  change_patype(aut, orig_patype); //transform back

  return primap;
}


}
