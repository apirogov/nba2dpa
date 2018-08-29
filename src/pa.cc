#include "pa.hh"
#include <vector>
#include <string>

namespace nbautils {
using namespace std;
using namespace nbautils;

bool PAProdState::operator==(PAProdState const& other) const {
  return a==other.a && b==other.b && priord==other.priord;
}

//compare lexicographically
bool PAProdState::operator<(PAProdState const& other) const {
  if (a == other.a) {
    if (b == other.b) {
        return priord < other.priord;
    }
    return b < other.b;
  }
  return a < other.a;
}

std::ostream& operator<<(std::ostream& os, PAProdState const& s) {
  os << "(" << s.a << "," << s.b << ") | ";
  for (int i=0; i<(int)s.priord.size(); i++) {
    os << "(" << (s.priord[i].first ? "r" : "l") << "," << s.priord[i].second << ")";
    if (i!=(int)s.priord.size()-1)
      os << ",";
  }
  return os;
}

PAProdState::PAProdState() {}

// construct an initial parity product state from initial state and parity bounds
// assuming min even parity!
PAProdState::PAProdState(state_t l, int lmin, int lmax, state_t r, int rmin, int rmax) : a(l), b(r) {
  if (lmin%2==1)
    ++lmin;
  if (lmax%2==1)
    ++lmax;
  if (rmin%2==1)
    ++rmin;
  if (rmax%2==1)
    ++rmax;

  for (int i=lmin; i<=lmax; i+=2) {
    priord.push_back(make_pair(false, i));
  }
  for (int i=rmin; i<=rmax; i+=2) {
    priord.push_back(make_pair(true, i));
  }
}

PAProdState::PAProdState(PAProdState const& other) : a(other.a), b(other.b), priord(other.priord) {}

// get new state from current with given new component states (and their prio) and adapted priorities
pair<PAProdState,int> PAProdState::succ(state_t l, int pl, state_t r, int pr, bool fulldown) const {
  PAProdState s(*this);
  int prio;

  //update state pair
  s.a = l;
  s.b = r;

  auto it = s.priord.begin();
  int i=0;
  while (it != s.priord.end()) {
    int p = it->first ? pr : pl;
    if (it->second == p) { //good priority fires -> we're done
      prio = 2*(i+1); //fire good
      break;
    } else if (it->second-1 == p) { //bad priority fires -> need shifting stuff
      prio = 2*i+1;  // fire bad

      //reshuffle priority tuples (preserving relative order)

      //mode: move completely down (looks like better choice)
      auto it2 = s.priord.end();
      if (!fulldown) {
      //mode: move one pair of other automaton above the red one
        it2 = it;
        while (it2 != s.priord.end() && it2->first == it->first)
          ++it2;
        if (it2 != s.priord.end())
          ++it2;
      }

      //perform the corresponding shifting
      stable_partition(it, it2, [b=it->first](auto t){ return t.first != b; });

      break;
    }

    //nothing interesting happened so just go on
    ++it;
    ++i;
  }

  return make_pair(s,prio);
}

}
