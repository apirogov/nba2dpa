#include "level.hh"
#include <iostream>
#include <algorithm>
#include <bitset>
#include <iostream>
#include "relorder.hh"
using namespace std;
using namespace nbautils;

namespace nbautils {

bool Level::operator<(Level const& other) const {
  if (tups == other.tups)
    return tupo < other.tupo;
  return tups < other.tups;

}

Level::hash_t add_powerset_hash(BA const& ba, Level const& lv) {
  Level::hash_t ret(0);

  set<Level::state_t> pset;
  for (auto& tup : lv.tups) copy(begin(tup), end(tup), inserter(pset, end(pset)));

  // traverse states in order, check which are present and set bits
  int num = 0;
  auto pit = pset.begin();
  for (auto const& it : ba.adj) {
    if (it.first == *pit) {
      ret.set(num);
      if (++pit == end(pset)) break;
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

// TODO: complete this
Level succ_level(LevelConfig const& lvc, Level l, sym_t x) {
  // cerr << "succ level" << endl;
  RelOrder rord(2 * l.tups.size());

  auto tupranks = rord.from_ranks(l.tupo);

  auto nsccs = l.tups.back();
  if (lvc.sep_rej) l.tups.pop_back();

  auto asccs = l.tups.back();
  auto aord = tupranks.back();
  if (lvc.sep_acc) {
    l.tups.pop_back();
    tupranks.pop_back();
  }


  // calculate successor sets, prioritize left
  // TODO: split by F/not F
  set<Level::state_t> used_sucs;
  vector<vector<Level::state_t>> suctups;
  for (auto& tup : l.tups) {
    auto suc = powersucc(*lvc.aut, tup, x);
    vector<Level::state_t> tmp;
    set_difference(begin(suc), end(suc), begin(used_sucs), end(used_sucs),
                   back_inserter(tmp));
    suctups.push_back(tmp);
    copy(begin(tmp), end(tmp), inserter(used_sucs, end(used_sucs)));
  }

  //pash rank of accepting component states, if any
  if (lvc.sep_acc) {
    tupranks.push_back(aord);
  }

  Level suclvl;
  suclvl.tups = suctups;
  suclvl.tupo = rord.to_ranks(tupranks);
  suclvl.powerset = add_powerset_hash(*lvc.aut, suclvl);
  // cerr << "end succ level" << endl;
  return suclvl;
}

}  // namespace nbautils
