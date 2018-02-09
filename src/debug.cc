#include "types.hh"
#include "debug.hh"
#include <iostream>

using namespace std;
using namespace nbautils;

namespace nbautils {

void printSCCI(nbautils::SCCInfo const& scci) {
  cout << "total number of SCCs: " << scci.sccrep.size() << endl;
  cout << scci.accepting.size() << " accepting SCCs: " << endl;
  for (auto i : scci.accepting)
    cout << i << ", ";
  cout << endl;
  cout << scci.rejecting.size() << " rejecting SCCs: " << endl;
  for (auto i : scci.rejecting)
    cout << i << ", ";
  cout << endl;
  cout << "Dead SCCs: ";
  for (auto it : scci.dead)
	  if (it.second)
			cout << it.first << ", ";
  cout << "." << endl;
  cout << scci.unreachable.size() << " unreachable states." << endl;
}

void printAcc(bool b) { cout << "*"; }
void printAcc(priority_t p) { cout << "{"<<p<<"}"; }

void printBA(nbautils::BA const& aut, nbautils::SCCInfo const& scci = SCCInfo()) {
  cout << "Name: " << aut.meta.name << ", APs: ";
  for (auto ap : aut.meta.aps)
    cout << ap << " ";
  cout << endl;
  cout << "BA with " << aut.adj.size() << " states:" << endl;
  for (auto &it : aut.adj) {
    auto s = it.first;
    cout << "State: " << s;
    if (aut.acc.find(s)!=end(aut.acc))
      printAcc(aut.acc.at(s));
	if (scci.unreachable.find(s)==end(scci.unreachable)) {
	  if (scci.scc.find(s)!=end(scci.scc)) {
		cout << " -- SCC: " << scci.scc.at(s);
      }
    }
	cout << endl;

    for (auto i = 0; i < aut.meta.num_syms; i++) {
      auto sucs = aut.succ(s, i);
      if (sucs.empty())
        continue;
      cout << "\t" << i << " -> ";
      for (auto t : sucs)
        cout << t << ", ";
      cout << endl;
    }
  }
}


}
