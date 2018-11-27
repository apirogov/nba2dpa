#include <iostream>
#include <random>
#include <bitset>
#include <set>
#include <vector>
#include <map>
#include <memory>

#include "common/types.hh"
#include "common/trie_map.hh"
#include "aut.hh"
#include "io.hh"
// #include "graph.hh"
// #include "pa.hh"
#include "preproc.hh"

// #include "metrics/bench.hh"
#include "common/util.hh"

#include "range/v3/all.hpp"

using namespace std;
using namespace nbautils;

nba_bitset vec_to_bitset(vector<int> const& v) { return to_bitset<nba_bitset>(v); }

//this is for trying stuff out. do whatever you want here.
int main(int argc, char *argv[]) {
  (void)argv[0];
  (void)argc;
  auto const log = spd::stderr_logger_mt("log");

  auto auts = nbautils::AutStream<Aut<string>>("");
  while (auts.has_next()) {
    auto aut = auts.parse_next();

    if (!aut.is_sba()) {
      cout << "Not Büchi!" << endl;
      return 1;
    }

    ba_trim(aut, log);
    aut.normalize();
    print_aut(aut);

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
