#include <iostream>
#include <random>
#include <bitset>
#include <set>
#include <vector>

#include "aut.hh"
#include "io.hh"
#include "graph.hh"
#include "pa.hh"
#include "preproc.hh"

#include "metrics/bench.hh"
#include "common/util.hh"

#include "range/v3/all.hpp"

using namespace std;
using namespace nbautils;

//this is for trying stuff out. do whatever you want here.
int main(int argc, char *argv[]) {
  (void)argv[0];
  (void)argc;

  auto auts = nbautils::AutStream<Aut<string>>("");
  while (auts.has_next()) {
    auto aut = auts.parse_next();

    map<unsigned, set<unsigned>> po;
    auto const simret = ba_direct_sim(aut);
    aut = simret.first;
    po = simret.second;
    for (auto const it : po) {
      cerr << it.first << " <=" << seq_to_str(it.second) << endl;
    }
    print_aut(aut, cerr);


    //ensure its a TDPA
    // aut.make_colored();
    // aut.to_tba();

    // auto aut2(aut);
    // complement_pa(aut2);

    // print_aut(aut);
    // print_aut(aut2);

    // auto paut = pa_union(aut, aut2);
    // complement_pa(paut);
    // print_aut(paut);
    // cout << pa_is_empty(paut) << endl;

    // print_aut(aut);
    // cout << pa_is_empty(aut) << endl;
    // auto run = get_acc_pa_run(aut);
    // cout << seq_to_str(run.first) << " | " << seq_to_str(run.second) << endl;
  }
}
