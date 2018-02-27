#include <vector>
#include <cassert>
#include <functional>
#include "swa.hh"
#include "common/util.hh"

namespace nbautils {
using namespace std;

function<acc_t(acc_t)> flip_acc_parity(vector<acc_t> const& as);
function<acc_t(acc_t)> flip_acc_polarity(vector<acc_t> const& as);
function<acc_t(acc_t)> priority_transformer(PAType from, PAType to, vector<acc_t> const& as);

//linear time priority reduction attempt
template<typename T>
function<acc_t(acc_t)> heuristic_minimize_priorities(SWA<Acceptance::PARITY,T> const& aut) {
  assert(aut.get_init().size() == 1);
  assert(is_colored(aut));

  auto const oldpris = aut.get_accsets();
  auto const min_even_adapter = priority_transformer(aut.get_patype(), PAType::MIN_EVEN, oldpris);

  //construct map from each reachable priority to biggest smaller priority
  set<acc_t> reachpri;
  map<acc_t, acc_t> maxlesspri;
  bfs(aut.get_init().front(), [&](auto const& st, auto const& visit, auto const&) {
      auto const stp = min_even_adapter(aut.get_accs(st).front());
      reachpri.emplace(stp);
      for (auto const& suc : aut.succ(st)) {
        auto const sucp = min_even_adapter(aut.get_accs(suc).front());
        if (sucp < stp) {
          if (!map_has_key(maxlesspri, stp)) {
            maxlesspri[stp] = sucp;
          } else {
            if (sucp > maxlesspri.at(stp))
              maxlesspri[stp] = sucp;
          }
        }
        visit(suc);
      }
   });

  // for (auto it : reachpri) {
  //   cout << it << " " << endl;
  // }
  // for (auto it : maxlesspri) {
  //   cout << it.first << " -> " << it.second << endl;
  // }

  //unreachable priorities can be set to any.
  //if exists biggest smaller reachable priority + has same parity, we don't need that one
  vector<acc_t> newpris;
  vector<acc_t> keptpris;
  vector<acc_t> useless;
  acc_t currpri = oldpris.front() % 2 == 0 ? 0 : 1;
  for (auto const& p : oldpris) {
    auto const adp = min_even_adapter(p);
    if (contains(reachpri, adp)) {
      if ((!map_has_key(maxlesspri, adp) || !same_parity(adp, maxlesspri.at(adp))) && (newpris.empty() || !same_parity(adp,newpris.back()))) {
        newpris.push_back(currpri++);
        keptpris.push_back(adp);
      } else {
        useless.push_back(adp);
      }
    } else {
      useless.push_back(adp);
    }
  }

  // cout << "kept: " << seq_to_str(keptpris) << endl;
  // cout << "removed: " << seq_to_str(useless) << endl;

  map<acc_t, acc_t> minprimap;
  for (int i=0; i<(int)keptpris.size(); i++)
    minprimap[keptpris.at(i)] = newpris.at(i);

  int j=keptpris.size()-1;
  for (int i=(int)useless.size()-1; i>=0; i--) {
    auto const adp = useless.at(i);
    while (adp < keptpris.at(j) && j>0)
      --j;
    minprimap[adp] = newpris.at(j);
  }

  //transform to original acceptance
  auto reverse_adapter = priority_transformer(PAType::MIN_EVEN, aut.get_patype(), newpris);
  //return map from original acceptance to minimized original acceptance
  return [=](acc_t p){ return reverse_adapter(minprimap.at(min_even_adapter(p))); };
}

template<typename T>
function<acc_t(acc_t)> minimize_priorities(SWA<Acceptance::PARITY,T> const& aut) {
  assert(is_colored(aut));

  function<vector<state_t>(state_t)> const sucs = [&](state_t v){ return aut.succ(v); };
  function<unsigned(state_t)> const get_pri = [&](state_t v){return aut.get_accs(v).front();};

  auto const oldpris = aut.get_accsets();
  auto const max_odd_adapter = priority_transformer(aut.get_patype(), PAType::MAX_ODD, oldpris);

  function<unsigned(state_t)> get_max_odd_pri = [&](state_t v){ return max_odd_adapter(get_pri(v)); };
  auto const primap = minimize_priorities(aut.states(), sucs, get_max_odd_pri);

  //transform to original acceptance
  vector<acc_t> tmppris = vec_fmap(aut.get_accsets(), max_odd_adapter);
  vec_to_set(tmppris);

  vector<acc_t> minpris = vec_fmap(tmppris, [&](auto p){return primap.at(p);});
  vec_to_set(minpris);

  auto reverse_adapter = priority_transformer(PAType::MAX_ODD, aut.get_patype(), minpris);
  //return map from original acceptance to minimized original acceptance
  return [=](acc_t p){ return reverse_adapter(primap.at(max_odd_adapter(p))); };
}

//map over the priorities (using a function obtained with other functions provided here)
template<typename T>
void transform_priorities(SWA<Acceptance::PARITY,T> &aut, function<acc_t(acc_t)> const& pf) {
  assert(is_colored(aut));
  for (auto const s : aut.states()) {
    auto const p = aut.get_accs(s).front();
    aut.set_accs(s, {pf(p)});
  }
}

//switch between different parity conditions
template<typename T>
void change_patype(SWA<Acceptance::PARITY,T> &aut, PAType pt) {
  auto f = priority_transformer(aut.get_patype(), pt, aut.get_accsets());
  transform_priorities(aut, f);
  aut.set_patype(pt);
}

}
