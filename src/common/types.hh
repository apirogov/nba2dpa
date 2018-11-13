#pragma once

#include <iostream>
#include <vector>
#include <bitset>
#include "common/util.hh"

namespace nbautils {
using namespace std;

// type for "small" automata (as the input should be)
constexpr size_t max_nba_states = 256;
using nba_bitset = bitset<max_nba_states>;

// struct nba_bitset_cmp {
//     bool operator() (const nba_bitset &b1, const nba_bitset &b2) const {
//         return b1.to_ullong() < b2.to_ullong();
//     }
// };

// type for alphabet (<=8 aps with all <=2^8=256 combinations)
using sym_t = uint8_t;
constexpr int max_nba_syms = 1 << (8 * sizeof(sym_t));

// type for automata states
using state_t = uint32_t;

// type of priorities
using pri_t = int;

//convert a rank + event to corresponding fired prio
inline pri_t rank_to_prio(pri_t r, bool good) {
  return 2*(r+1)-(good ? 0 : 1);
}

inline pair<pri_t, bool> prio_to_event(pri_t pri) {
  return make_pair((pri+1)/2 - 1, pri%2==0);

}

// type of ranked slices
// all n nba_bitsets should be pw. disjoint (and non-empty if it is not a pre-slice)
// all pri_t values should be distinct numbers
// the pri_t values should be 0 <= p < n (unless it is just a component of a set of such slices)
using ranked_slice = vector<pair<nba_bitset, pri_t>>;

// dual rep. of ranked slices
// S(0) is any set. for i>1 we have:
// forall j<i either S(i) subset S(j) or S(i),S(j) have empty intersection
using tree_history = vector<nba_bitset>;

std::ostream& operator<<(std::ostream& os, tree_history const& th);

std::ostream& operator<<(std::ostream& os, ranked_slice const& rs);
// print a vector of ranked slices, i.e., different components
std::ostream& operator<<(std::ostream& os, vector<ranked_slice> const& rs);

// switch between redundant and non-redundant label sets
ranked_slice unprune(ranked_slice const& rslice);
ranked_slice prune(ranked_slice const& rslice);

// takes unpruned ranked slices and sorts it by rank.
tree_history in_rank_order(ranked_slice const& urs);

// calculate parent and left sibling pos in tuple
pair<vector<int>,vector<int>> unflatten(ranked_slice const& rank);

// do not use this in DetState! only works correctly on "pure" structs
tree_history slice_to_history(ranked_slice const& rs);
ranked_slice history_to_slice(tree_history const& th);

bool finer_or_equal(ranked_slice const& rs1, ranked_slice const& rs2);
pair<nba_bitset,vector<nba_bitset>> kcut_mask(tree_history const& th, int k=0);

// these just demonstrate the idea. not actually used
bool kequiv(tree_history const& th1, tree_history const& th2, int k=0);
bool not_worse(tree_history const& th1, tree_history const& th2, int k=0);

}
