#include <iostream>
#include "io.hh"
#include "common/algo.hh"
#include "common/part_refinement.hh"
#include "pa.hh"

using namespace std;
using namespace nbautils;

//run with p==pa==pb
bool equiv(auto const& aut, int k, int pa, int pb, state_t a, state_t b) {
  for (int i=0; i<aut.num_syms(); i++) {
    auto suca = aut.succ(a, i).front();
    auto sucb = aut.succ(b, i).front();
    auto sucpa = min(pa, (int)aut.get_accs(suca).front());
    auto sucpb = min(pb, (int)aut.get_accs(sucb).front());

    cerr << "-> " << suca << "," << sucb << " " << sucpa << ":" << sucpb << endl;

    if (k==1) {
      if (suca != sucb || sucpa != sucpb)
        return false;
    } else {
      if (suca == sucb && sucpa != sucpb)
        return false;

      cerr << suca << "," << sucb << " " << sucpa << ":" << sucpb << endl;
      if (!equiv(aut, k-1, sucpa, sucpb, suca, sucb))
        return false;
    }
  }

  //all returned true -> true
  cerr << "success " << a <<","<<b<<endl;
  return true;
}

vector<vector<state_t>> pa_k_equiv(auto const& aut, int k) {
  auto const sts = aut.states();
  map<state_t, state_t> eq;
  for (auto const st : sts)
    eq[st] = st;

  for (int i=0; i<(int)sts.size(); ++i) {
    for (int j=i+1; j<(int)sts.size(); ++j) {
      // if (aut.get_accs(sts[i]) != aut.get_accs(sts[j]))
      //   continue; //different priorities can not be k-equivalent

      // int pi = aut.get_accs(sts[i]).front();
      // int pj = aut.get_accs(sts[j]).front();
      int p = aut.get_accsets().back();
      // cerr << "---" << endl;
      // cerr << sts[i] << "," << sts[j] << " " << pi << ":" << pj << endl;
      if (!equiv(aut, k, p, p, sts[i], sts[j]))
        continue;

      eq[sts[j]] = eq[sts[i]];
    }
  }

  return group_by(sts, [&](state_t a, state_t b){ return eq.at(a)==eq.at(b); });
}

//this is for trying stuff out. do whatever you want here.
int main(int argc, char *argv[]) {
  (void)argv[0];
  bool intersect=argc>1;
  auto auts = parse_hoa("");

  if (auts.empty())
    cout << "something went wrong!" << endl;

  auto &aut1 = *auts[0];
  auto &aut2 = *auts[1];

  // cerr << "is empty? " << pa_is_empty(aut) << endl;
  // cerr << "good set: " << seq_to_str(find_acc_pa_scc(aut)) << endl;
  // vector<state_t> pref;
  // vector<state_t> cyc;
  // tie(pref, cyc) = get_acc_pa_run(aut);
  // cerr << seq_to_str(pref) << " | " << seq_to_str(cyc) << endl;

  auto aut_prod = pa_prod(aut1, aut2, intersect);
  minimize_priorities(*aut_prod);
  minimize_pa(*aut_prod);
  print_hoa(*aut_prod);

  auto eqv = pa_equiv_states(*aut_prod);
  for (auto eq : eqv)
    cerr << seq_to_str(eq) << endl;

  eqv = pa_k_equiv(*aut_prod, 5);
  for (auto eq : eqv)
    cerr << seq_to_str(eq) << endl;
}
