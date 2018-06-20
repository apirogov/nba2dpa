#include <iostream>

#include "aut.hh"
#include "io.hh"

using namespace std;
using namespace nbautils;

//this is for trying stuff out. do whatever you want here.
int main(int argc, char *argv[]) {
  (void)argv[0];
  (void)argc;

  auto auts = nbautils::AutStream<Aut<string>>("");

  while (auts.has_next()) {
    auto aut = auts.parse_next();

    int sigloops = 0;
    int accsigloops = 0;
    for (auto const p : aut.states()) {
      bool sigmaloop = true;
      for (auto const x : aut.syms()) {
        if (!contains(aut.succ(p,x), p)) {
          sigmaloop = false;
          break;
        }
      }
      if (sigmaloop) {
        sigloops++;
        if (aut.state_buchi_accepting(p))
          accsigloops++;
      }
    }

    cout << aut.num_states() << "\t" << sigloops << "\t" << accsigloops << endl;
  }

}
