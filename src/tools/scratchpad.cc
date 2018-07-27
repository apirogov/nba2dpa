#include <iostream>
#include <random>
#include <bitset>
#include <set>
#include <vector>

#include "aut.hh"
#include "io.hh"
#include "graph.hh"
#include "pa.hh"

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
    aut.make_colored();
    minimize_priorities(aut);
    if (aut.get_patype() != PAType::MIN_EVEN)
      change_patype(aut, PAType::MIN_EVEN);
    print_aut(aut);
  }
}
