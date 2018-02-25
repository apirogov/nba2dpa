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

//TODO: debug this. does not work correctly yet
template<typename T>
function<acc_t(acc_t)> minimize_priorities(SWA<Acceptance::PARITY,T> const& aut) {
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
  for (auto it : reachpri) {
    cout << it << " " << endl;
  }
  for (auto it : maxlesspri) {
    cout << it.first << " -> " << it.second << endl;
  }

  //unreachable priorities can be set to any.
  //if exists biggest smaller reachable priority + has same parity, we don't need that one
  map<acc_t, acc_t> minprimap;
  vector<acc_t> newpris;
  acc_t currpri = oldpris.front() % 2 == 0 ? 0 : 1;
  for (auto const& p : oldpris) {
    auto const adp = min_even_adapter(p);
    if (contains(reachpri, adp)) {
      if ((!map_has_key(maxlesspri, adp) || !same_parity(adp, maxlesspri.at(adp))) && (newpris.empty() || !same_parity(adp,newpris.back()))) {
        newpris.push_back(currpri++);
      }
    }
    minprimap[adp] = newpris.back();
  }
  //transform to original acceptance
  auto reverse_adapter = priority_transformer(PAType::MIN_EVEN, aut.get_patype(), newpris);
  //return map from original acceptance to minimized original acceptance
  return [=](acc_t p){ return reverse_adapter(minprimap.at(min_even_adapter(p))); };
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
