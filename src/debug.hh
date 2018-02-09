#pragma once

#include "types.hh"
#include "io.hh"
#include "scc.hh"
#include "ps.hh"
#include "level.hh"

namespace nbautils {

void printSCCI(nbautils::SCCInfo const& scci);
void printBA(nbautils::BA const& aut, nbautils::SCCInfo const& scci);

template<typename L>
void printPS(nbautils::PS<L>& aut,nbautils::SCCInfo &scci, bool pointed=false) {
  using namespace std;
  for (auto &it : aut.adj) {
    auto s = it.first;
    cout << s;
      auto pacc = aut.acc.find(s);
      string acc = "";
      if (pacc!=end(aut.acc)) {
        cout << "(" << pacc->second << ")";
      }
    cout << " - ";

    auto tag = aut.tag->get(s);
    if (pointed) {
      nbautils::state_t ptd = tag.back();
      tag.pop_back();
      cout << "(" << ptd << ", {";
    } else {
      cout << "{";
    }
    for (nbautils::state_t si : tag)
      cout << si << ",";
    cout << "}";
    if (pointed)
      cout << ")";

	if (scci.unreachable.find(s)==end(scci.unreachable))
		cout << " -- SCC: " << scci.scc[s];

	cout << endl;

    for (auto i = 0; i < aut.num_syms; i++) {
      if (aut.adj[s][i].empty())
        continue;
      cout << "\t" << i << " -> ";
      for (auto t : aut.adj[s][i])
        cout << t << ", ";
      cout << endl;
    }
  }
}

}
