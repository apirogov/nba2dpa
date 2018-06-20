#pragma once
#include <iostream>

#include "aut.hh"
#include "common/scc.hh"
#include "ps.hh"

#include <spdlog/spdlog.h>
namespace spd = spdlog;

namespace nbautils {
using namespace std;

//returns sorted list of accepting sinks (acc. states with self-loop for each sym)
template <typename T>
set<state_t> ba_get_acc_sinks(Aut<T> const& aut, shared_ptr<spdlog::logger> log=nullptr) {
  assert(aut.is_buchi());

  set<state_t> ret;
  for (auto const v : aut.states()) {
    auto const outsyms = aut.state_outsyms(v);
    //must be accepting and have successors for each symbol
    bool accsink = aut.state_buchi_accepting(v) && outsyms.size() == aut.num_syms();
    for (auto const i : aut.syms()) {
        // must have self-loop for each symbol
        if (!contains(aut.succ(v,i), v))
          accsink = false;
    }
    if (accsink)
      ret.emplace(v);
  }

  if (log)
    log->info("found {} accepting sinks", ret.size());
  return ret;
}

// SCC acceptance class: -1 = rejecting, 0 = mixed, 1 = accepting
using BASccAClass = map<unsigned, int>;

// classify sccs as accepting or rejecting
template <typename T>
BASccAClass ba_scc_classify_acc(Aut<T> const& aut, SCCDat const& scci) {
  BASccAClass ret;
  for (auto& scc : scci.sccs) {
    bool rej = true;
    bool acc = true;
    for (auto const s : scc.second) {
      if (aut.state_buchi_accepting(s)) //contains acc state
        rej = false; //-> not NSCC
      else //contains nonacc state
        acc = false; //-> not ASCC

      if (!rej && !acc) //we're done, have witnesses for both
        break;
    }
    if (acc)
      ret[scc.first] = 1;  //must be ASCC
    else if (rej)
      ret[scc.first] = -1; //must be NSCC
    else
      ret[scc.first] = 0;  //must be MSCC
  }

  return ret;
}

template <typename F>
void mark_dead_sccs(BASccAClass const& sccacl, F const& get_suc_sccs,
    map<unsigned,bool>& dead, unsigned num) {
  if (map_has_key(dead, num))  // done already
    return;

  // std::cout << "scc " << num << std::endl;
  auto sucsccs = get_suc_sccs(num);
  // mark children first
  for (auto sucnum : sucsccs)
    mark_dead_sccs(sccacl, get_suc_sccs, dead, sucnum);

  // if we are rejecting, assume we're dead
  bool isdead = sccacl.at(num) == -1;
  // check children and try to falsify
  for (auto sucscc : sucsccs)
    isdead = isdead && dead[sucscc];
  // if still dead, we're really dead.
  dead[num] = isdead;
}

// run dfs that marks dead sccs (assuming trivial SCCs never accepting)
template <typename T>
set<unsigned> ba_get_dead_sccs(Aut<T> const& ba, SCCDat const& scci, BASccAClass const& sccacl) {
  map<unsigned, bool> dead;
  auto const scc_succ = [&](unsigned const scc){ return succ_sccs(aut_succ(ba), scci, scc); };

  for (auto const i : ranges::view::keys(sccacl))
    mark_dead_sccs(sccacl, scc_succ, dead, i);

  return mapbool_to_set(dead);
}

// return deterministic SCCs (only det. transitions inside, can have nondet to outside)
template <typename T>
set<unsigned> ba_scc_classify_det(Aut<T> const& aut, SCCDat const& scci) {
  set<unsigned> ret;

  for (auto& scc : scci.sccs) {
    bool det = true;

    for (auto const s : scc.second) {
      for (auto const x : aut.state_outsyms(s)) {
        auto const sxsucs = aut.succ(s,x);
        if (sxsucs.size() < 2)
          continue; //at most one suc. anyway -> nice

        bool has = false;
        for (auto const suc : sxsucs) {
          auto const sucscc = scci.scc_of.at(suc);
          if (sucscc != scc.first)
            continue; //other target SCCs don't matter, we care about inside determinism

          if (has) { // s already has x-successor in current SCC
            det = false;
            goto detdone;
          }

          has = true; //mark in-SCC successor as found
        }
      }
    }

detdone:
    if (det)
      ret.emplace(scc.first);
  }

  return ret;
}

// find and purge unreachable and dead states from automaton
// (but keeps initial states, even when they are dead)
template <typename T>
void ba_trim(Aut<T>& ba, shared_ptr<spdlog::logger> log=nullptr) {
  assert(ba.is_buchi());
  size_t trimmed = 0;

  auto const unreach = unreachable_states(ba, ba.get_init());
  trimmed += unreach.size();
  // cerr << "unreach: " << unreach.size() << endl;
  ba.remove_states(unreach);

  // ----

  auto const ba_suc = aut_succ(ba);
  auto const scci = get_sccs(ba.states() | ranges::to_vector, ba_suc);

  // cerr << "#sccs: " << scci.sccs.size() << endl;

  //unmark trivial SCCs, if marked as acc. (which is useless)
  auto const sccTrv = trivial_sccs(ba_suc, scci);
  for (auto const scc : sccTrv) {
    for (auto const s : scci.sccs.at(scc))
      if (ba.has_pri(s)) {
        ba.set_pri(s, -1);
        if (log)
          log->info("trivial state {} made rejecting", s);
      }
  }

  //now classify (to get rejecting SCCs)
  auto const sccAcl = ba_scc_classify_acc(ba, scci);

  //now detect dead SCCs (rejecting that do not reach non-rej SCCs)
  auto const deadscc = ba_get_dead_sccs(ba, scci, sccAcl);

  // collect dead states
  set<state_t> dead;
  for (auto const s : ba.states()) {
      auto const scit = scci.scc_of.at(s);
      if (contains(deadscc, scit) && s != ba.get_init()) {  // is a dead scc
        dead.emplace(s);
      }
  }
  //remove them
  trimmed += dead.size();
  ba.remove_states(vector<state_t>(cbegin(dead),cend(dead)));

  if (log)
    log->info("removed {} useless states", trimmed);
}

//context state set -> relatively acc, rej subsets
using Context = unordered_map<nba_bitset, pair<nba_bitset, nba_bitset>>;

Context get_context(auto const& aut, adj_mat const& mat, nba_bitset const asinks,
                    shared_ptr<spdlog::logger> log = nullptr) {

  auto const psp = bench(log,"powerset_product",
                         WRAP(powerset_product(aut, mat, asinks)));
  auto const psp_scci = get_sccs(psp.states() | ranges::to_vector, aut_succ(psp));
  auto const psp_sccAcc = ba_scc_classify_acc(psp, psp_scci);
  // print_aut(psp);

  if (log)
    log->info("#states in 2^AxA: {}, #SCCs in 2^AxA: {}",
              psp.num_states(), psp_scci.sccs.size());

  Context ret;
  for (auto const s : psp.states()) {
    auto const& t = psp.tag.geti(s);
    int val = psp_sccAcc.at(psp_scci.scc_of.at(s));
    if (val==1)
      ret[t.first].first[t.second] = 1;
    else if (val==-1)
      ret[t.first].second[t.second] = 1;
  }

  return ret;
}

inline void print_context(Context const& ctx) {
  for (auto const it : ctx) {
    vector<state_t> tmp;
    from_bitset(it.first, back_inserter(tmp));
    cout << seq_to_str(tmp) << " -> "; //context (curr. reach. set)
    tmp.clear();
    from_bitset(it.second.first, back_inserter(tmp));
    cout << seq_to_str(tmp) << " | ";  //relative accepting st.
    tmp.clear();
    from_bitset(it.second.second, back_inserter(tmp));
    cout << seq_to_str(tmp) << endl;   //relative rejecting st.
  }
}

}

