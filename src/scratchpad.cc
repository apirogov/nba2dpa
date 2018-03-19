#include <iostream>
#include "io.hh"
#include "common/algo.hh"
#include "common/part_refinement.hh"
#include "pa.hh"

#include <range/v3/all.hpp>

using namespace std;
using namespace nbautils;

//this is for trying stuff out. do whatever you want here.
int main(int argc, char *argv[]) {
  string file="";
  if (argc>1)
    file = argv[1];
  auto bas = parse_hoa(file);
  if (bas.empty())
    cout << "something went wrong!" << endl;

  auto &aut = *bas.front();

  ba_to_pa(aut);
  for (auto s : aut.states())
    cout << s << " -> " << seq_to_str(aut.get_accs(s)) << endl;

  cout << "BA is empty? " << pa_is_empty(aut) << endl;
  cout << "good set: " << seq_to_str(find_acc_pa_scc(aut)) << endl;
  vector<state_t> pref;
  vector<state_t> cyc;
  tie(pref, cyc) = get_acc_pa_run(aut);
  cout << seq_to_str(pref) << " | " << seq_to_str(cyc) << endl;


}
