#include "level.hh"
#include "relorder.hh"
#include <bitset>
#include <iostream>
#include <algorithm>
using namespace std;
using namespace nbautils;

namespace nbautils {

Level::hash_t add_powerset_hash(BA const& ba, Level const& lv) {
  Level::hash_t ret(0);

  set<Level::state_t> pset;
  for (auto &tup : lv.tups)
    copy(begin(tup),end(tup), inserter(pset, end(pset)));

  //traverse states in order, check which are present and set bits
  int num=0;
  auto pit = pset.begin();
  for (auto const& it : ba.adj) {
    if (it.first == *pit) {
      ret.set(num);
      if (++pit == end(pset))
        break;
    }
    num++;
  }

  return ret;
}

Level make_level(LevelConfig const& lvc, std::vector<Level::state_t> const& qs) {
  Level l;
  l.tups.push_back(qs);
  l.tupo.push_back(0);
  if (lvc.sep_acc) {
    l.tups.push_back({});
    l.tupo.push_back(1);
  }
  if (lvc.sep_rej) {
    l.tups.push_back({});
  }
  l.powerset = add_powerset_hash(*lvc.aut, l);
  return l;
}

Level succ_level(LevelConfig const& lvc, Level l, sym_t x) {
  RelOrder rord(2*l.tups.size());

  auto nsccs = l.tups.back();
  if (lvc.sep_rej)
    l.tups.pop_back();

  auto asccs = l.tups.back();
  auto aord = l.tupo.back();
  if (lvc.sep_acc) {
    l.tups.pop_back();
    l.tupo.pop_back();
  }

  auto suctupo = rord.from_ranks(l.tupo);

  // calculate successor sets, prioritize left
  set<Level::state_t> used_sucs;
  vector<vector<Level::state_t>> suctups;
  for (auto &tup : l.tups) {
    auto suc = powersucc(*lvc.aut, tup, x);
    vector<Level::state_t> tmp;
    set_difference(begin(suc),end(suc),begin(used_sucs),end(used_sucs),back_inserter(tmp));
    suctups.push_back(tmp);
    copy(begin(tmp),end(tmp),inserter(used_sucs, end(used_sucs)));
  }

  Level suclvl;
  suclvl.tups = suctups;
  suclvl.tupo = rord.to_ranks(suctupo);
  suclvl.powerset = add_powerset_hash(*lvc.aut, suclvl);
  return suclvl;
}

}
