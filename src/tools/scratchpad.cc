#include <iostream>
#include <random>
#include <bitset>
#include <set>
#include <vector>
#include <map>

#include "common/types.hh"
#include "common/trie_map.hh"
// #include "aut.hh"
// #include "io.hh"
// #include "graph.hh"
// #include "pa.hh"
// #include "preproc.hh"

// #include "metrics/bench.hh"
#include "common/util.hh"

#include "range/v3/all.hpp"

using namespace std;
using namespace nbautils;




//TODO: implement nba_bitset-labeled trie where each node can be marked.
//a Q-slice-trie stores a set of ranked slices that have exactly Q as comlete set
//a trie branch up to a marked node uniquely determines presence of a ranked slice
//given trie node n on level k, all ranked slices in the subtree rooted at n are k-equivalent

nba_bitset vec_to_bitset(vector<int> const& v) { return to_bitset<nba_bitset>(v); }

//this is for trying stuff out. do whatever you want here.
int main(int argc, char *argv[]) {
  (void)argv[0];
  (void)argc;

  auto rslice = ranked_slice{ make_pair(vec_to_bitset({3}), 5)
                            , make_pair(vec_to_bitset({2}), 2), make_pair(vec_to_bitset({5}), 4)
                            , make_pair(vec_to_bitset({4}), 3), make_pair(vec_to_bitset({1}), 1) };
  auto rslice2= ranked_slice{ make_pair(vec_to_bitset({2,3}), 2)
                            , make_pair(vec_to_bitset({4,5}), 3), make_pair(vec_to_bitset({1}), 1) };

  auto const irslice = slice_to_history(rslice);

  auto const irslice2 = slice_to_history(rslice2);

  // do this for every slice in state... concat them... then create unique different rep :
  for (auto b : irslice)
    cout << pretty_bitset(b);
  cout << endl;
  for (auto b : irslice2)
    cout << pretty_bitset(b);
  cout << endl;

  cout << "----" << endl;

  auto const rslice3 = history_to_slice(irslice);
  cout << rslice3 << endl;

  cout << "----" << endl;

  trie_map<nba_bitset, ranked_slice> trie;
  cout << trie.size() << " ";

  trie.put(irslice, rslice);
  cout << "after put: " << trie.size() << " ";
  cout << "has exact: " << trie.has(irslice2) << " ";
  cout << "has any: " << trie.has(irslice2, true) << " ";
  trie.put(irslice, rslice);
  cout << "after same put: " << trie.size() << " ";
  cout << endl;

  cout << "get exact:" << endl;
  ranked_slice fromtrie1 = trie.get(irslice);
  cout << fromtrie1 << endl;

  cout << "get any:" << endl;
  ranked_slice fromtrie2 = trie.get(irslice2, true);
  cout << fromtrie2 << endl;

  trie.put(irslice2, rslice2);
  cout << "after other put: " << trie.size() << " ";
  cout << "has exact: " << trie.has(irslice2) << " ";
  cout << endl;

  cout << "get any:" << endl;
  ranked_slice fromtrie3 = trie.get(irslice2, true);
  cout << fromtrie3 << endl;

  /*
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
    */


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
  // }

}
