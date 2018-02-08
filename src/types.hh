#pragma once

#include "triebimap.hh"
#include <iostream>
#include <cinttypes>
#include <map>
#include <memory>
#include <set>
#include <algorithm>
#include <unordered_set>
#include <string>
#include <utility>
#include <vector>

namespace nbautils {

struct ParsedMeta {
  // automaton parsed successfully
  bool success = false;
  // atomic props
  std::vector<std::string> aps;
  // automaton name
  std::string name;
};

// type for "small" automata (as the input should be)
using small_state_t = uint8_t;
constexpr size_t max_nba_states = 1<<(8*sizeof(small_state_t));

// type for alphabet (<=8 aps with all <=2^8=256 combinations)
using sym_t = uint8_t;
constexpr size_t max_nba_syms = 1<<(8*sizeof(sym_t));

// type for "big" automata
using state_t = uint32_t;

using scc_t = state_t;

template <typename T> using state_label = std::map<state_t, T>;
using state_flag = std::unordered_set<state_t>;


// acceptance-labelled KS with tagged nodes
template <typename L,typename T>
struct SWA {
  typedef std::shared_ptr<SWA<L,T>> ptr;

  size_t num_syms = 0; // number of actual symbols 0..(2^AP - 1)
  state_t init;        // initial state

  // list of successors per symbol
  std::map<state_t, std::map<sym_t, std::vector<state_t>>> adj;
  // state acceptance label
  state_label<L> acc;
  std::shared_ptr<T> tag;
};

template<typename L,typename T>
std::vector<state_t> succ(SWA<L,T> const& ks, state_t const& p, sym_t const& x) {
	if (ks.adj.find(p) == end(ks.adj))
      return {};
	auto const& edges = ks.adj.at(p);
	if (edges.find(x) == end(edges))
      return {};
	return edges.at(x);
}

//given set in sorted(!!) vector, kill states
template<typename L,typename T>
void kill_states(SWA<L,T> &aut, std::vector<state_t> const& tokill) {
  for (auto const& it : tokill) {
	  auto el =aut.adj.find(it);
	  if (el != end(aut.adj))
		aut.adj.erase(el); //kill state + edges
  }
  for (auto const& it : tokill) {
	  auto el =aut.acc.find(it);
	  if (el != end(aut.acc))
		aut.acc.erase(el); //kill acc. mark
  }

  //kill states from edge targets
  for (auto& it : aut.adj) {
    for (auto &jt : it.second) {
      std::vector<state_t> res;
      std::set_difference(begin(jt.second), end(jt.second),
                          begin(tokill), end(tokill), std::back_inserter(res));
      jt.second = res;
    }
  }
}

//given set in sorted(!!) vector, merge states into one
template<typename L,typename T>
void merge_states(SWA<L,T> &aut, std::vector<state_t>& eqclass) {
	if (eqclass.size()<2)
		return;

	//get representative
	state_t rep = eqclass.back();

	//add outgoing edges from all class members to representative
	for (auto const& st : eqclass) {
		for (auto& symmap : aut.adj.at(st)) {
			auto& sym = symmap.first;
			auto& sucs = symmap.second;

			std::vector<state_t> res;
			std::set_union(begin(sucs),end(sucs),
					begin(aut.adj.at(rep)[sym]),end(aut.adj.at(rep)[sym]),
					back_inserter(res));
			aut.adj.at(rep).at(sym) = res;
		}
	}
	//add ingoing edges into all class members to representative
	for (auto &it : aut.adj) {
		auto &st = it.first;
		for (auto& symmap : aut.adj.at(st)) {
			auto& sym = symmap.first;
			auto& sucs = symmap.second;
			std::vector<state_t> res;
			set_intersection(begin(sucs),end(sucs),
					begin(eqclass),end(eqclass),
					back_inserter(res));
			if (res.size()>0) {
				res = std::vector<state_t>{};
				auto singl = std::vector<state_t>{rep};
				std::set_union(begin(sucs),end(sucs),
					begin(singl),end(singl),
					back_inserter(res));
			sucs = res;

			}
		}
	}

	//remove representative from kill list
	eqclass.pop_back();
	//kill the others!!
	kill_states(aut, eqclass);
}

template<typename L,typename T,typename S>
std::vector<small_state_t>
powersucc(SWA<L,T> const& ks, std::vector<S> const &ps, sym_t const &x) {
  std::set<state_t> suc;
  for (auto const& p : ps) {
    auto qs = succ(ks, p, x);
    std::copy(qs.begin(), qs.end(), std::inserter(suc, suc.end()));
  }
  return std::vector<S>(suc.begin(), suc.end());
}

using BA = SWA<bool,void>;

using ps_tag = trie_bimap<small_state_t, state_t>;
template<typename L>
using PS = SWA<L, ps_tag>;
using BAPS = PS<bool>;

using priority_t = small_state_t;
using PA = SWA<priority_t,void>;

using scc_tag = state_label<scc_t>;
struct SCCInfo {
  typedef std::shared_ptr<SCCInfo> ptr;
  scc_tag scc;         // tags each state with scc num
  std::map<scc_t,state_t> sccrep;      // representative state
  state_flag unreachable; // tags a state as unreachable from initial state

  state_flag accepting; // tags an scc as fully accepting
  state_flag rejecting; // tags an scc as fully rejecting
  state_flag trivial;   // tags an scc as trivial (single state, no loop)

  std::map<scc_t,state_t> dead;      // tags an scc as dead (means that all reachable states are
                     // rejecting)
};

}
