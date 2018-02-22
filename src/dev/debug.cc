#include "debug.hh"
#include "common/util.hh"
#include <iostream>

using namespace std;
using namespace nbautils;

namespace nbautils {

void printSCCI(nbautils::SCCInfo const& scci) {
  cout << "total number of SCCs: " << scci.sccrep.size() << endl;
  cout << scci.accepting.size() << " accepting SCCs: " << endl;
  for (auto i : scci.accepting) cout << i << ", ";
  cout << endl;
  cout << scci.rejecting.size() << " rejecting SCCs: " << endl;
  for (auto i : scci.rejecting) cout << i << ", ";
  cout << endl;
  cout << "Dead SCCs: ";
  for (auto it : scci.dead)
    if (it.second) cout << it.first << ", ";
  cout << "." << endl;
  cout << scci.unreachable.size() << " unreachable states." << endl;
}

void printBA(nbautils::BA const& aut, nbautils::SCCInfo const& scci = SCCInfo()) {
  cout << "Name: " << aut.get_name() << ", APs: ";
  for (auto ap : aut.get_aps()) cout << ap << " ";
  cout << endl;
  cout << "BA with " << aut.adj.size() << " states:" << endl;
  for (auto& it : aut.adj) {
    auto s = it.first;
    cout << "State: " << s;
    if (aut.has_accs(s)) {
      auto const& ac = aut.get_accs(s);
      cout << "{";
      for (auto a : ac) {
        cout << a << ",";
      }
      cout << "}";
    }
    if (!contains(scci.unreachable, s)) {
      if (map_has_key(scci.scc, s)) {
        cout << " -- SCC: " << scci.scc.at(s);
      }
    }
    cout << endl;

    for (auto i = 0; i < (int)aut.num_syms(); i++) {
      auto sucs = aut.succ(s, i);
      if (sucs.empty()) continue;
      cout << "\t" << i << " -> ";
      for (auto t : sucs) cout << t << ", ";
      cout << endl;
    }
  }
}

}  // namespace nbautils
