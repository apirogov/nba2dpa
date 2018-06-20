#include <iostream>
#include <random>
#include <bitset>
#include <set>
#include <vector>

#include "aut.hh"
#include "metrics/bench.hh"
#include "common/util.hh"

using namespace std;
using namespace nbautils;

//this is for trying stuff out. do whatever you want here.
int main(int argc, char *argv[]) {
  (void)argv[0];
  (void)argc;

  set<unsigned char> inp{2,3,5,7,11,13,255};
  nba_bitset b = to_bitset<nba_bitset>(inp);
  vector<unsigned char> ret;
  from_bitset(b, back_inserter(ret));
  cerr << seq_to_str(ret) << endl;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::bernoulli_distribution d(0.5);

  int const num=100000;
  int const st=256;

  using settype = bitset<256>;
  // using settype = vector<char>;
  // using settype = set<unsigned char>;
  vector<settype> bs;
  for (int i=0; i<num; i++) {
    settype bits;
    for(int n=0; n < st; ++n) {
      if (d(gen))
        bits[n] = 1;
        // bits.push_back(n);
        // bits.emplace(n);
    }
    // cerr << seq_to_str(bits) << endl;
    bs.push_back(bits);
  }

  auto starttime = get_time();
  cout << "start" << endl;

  settype bits;
  for (int i=0; i<num; i++) {
    bits |= bs[i];
    // settype tmp;
    // set_union(cbegin(bits), cend(bits), cbegin(bs[i]), cend(bs[i]),
    //           // std::inserter(tmp, end(tmp)));
    //           std::back_inserter(tmp));
    // bits = tmp;
  }

  cerr << "total time (s): " << get_secs_since(starttime) << endl;
}
