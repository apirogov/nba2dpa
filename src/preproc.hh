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
// unmark accepting trivial states, mark nonaccepting states in inherently weak SCCs
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
  auto const scci = get_sccs(ba.states(), ba_suc);

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

  //mark all states in inherently weak SCCs accepting
  for (auto const scc : scci.sccs) {
    if (sccAcl.at(scc.first) != 0)
      continue; //rejecting component anyway

    // forbid accepting states
    auto const nacc_suc = [&](state_t v) {
      auto const filt = [&](state_t const p){ return !ba.has_pri(v) && !ba.has_pri(p); };
      auto const sucs = ba_suc(v);
      // non-accepting can go to non-accepting
      return sucs | ranges::view::filter(filt) | ranges::to_vector;
    };
    // check for non-trivial SCCs
    auto const innerscci = get_sccs(scc.second, nacc_suc);
    auto const innerTrv = trivial_sccs(nacc_suc, innerscci);

    if (innerTrv.size() == innerscci.sccs.size()) {
      //all trivial -> inherently weak -> mark all as accepting
      if (log)
      log->info("inherently weak component {} made weak", seq_to_str(scc.second));
      for (auto const s : scc.second)
        ba.set_pri(s, 0); //mark accepting
    }
  }

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

// direct simulation Algo. Paper: "Optimizing Buchi automata" (Etessami, Holzmann)
// TODO: review code, maybe make more efficient (not critical, as inputs are small)
template <typename T>
auto ba_direct_sim(Aut<T> const& ba) {
  using po_type = map<unsigned,set<unsigned>>;
  using ni_type = map<int,set<unsigned>>; //maximal neighbor i-type set N_i
  using interm_clr = pair<unsigned, ni_type>;  // <c_i-1, N_i-1>


  //color a <=_po color b
  auto const po_leq = [](po_type const& po, unsigned a, unsigned b) {
    if (a==b) //refl.
      return true;
    if (!map_has_key(po, a))
      return false;
    if (!contains(po.at(a), b)) // a <= b and a != b ==> a < b
      return false;
    return true;
  };

  //for each sym of b there is this sym of a with >= color <=> a dominates b
  auto const dominates = [&po_leq](po_type const& po, ni_type const& a, ni_type const& b) {
    for (auto const& it : b) {
      if (!map_has_key(a, it.first))
        return false; // a has not pairs for this symbol -> can not dominate b

      auto const domcand = a.at(it.first);
      for (auto const bc : it.second) { //for each maximal color for sym in b...
        //look for a maximal color in a that dominates that
        bool has_dom = false;
        for (auto const ac : domcand) {
          if (po_leq(po, bc, ac)) {
            has_dom = true;
            break;
          }
        }
        //failed to find dominating color
        if (!has_dom)
          return false;
      }
    }
    return true;
  };

  po_type old_po; //init po_-1
  po_type po{{0,{0}},{1,{0,1}}}; //init po_0
  map<state_t, unsigned> old_clr; //init C_-1
  map<state_t, unsigned> clr; //init C_0
  for (auto const s : ba.states())
    clr[s] = !ba.state_buchi_accepting(s);
  map<state_t, interm_clr> expanded_clrs; //states -> tuple-colors

  //iterate until fixpoint
  while (clr != old_clr || po != old_po) {
    //shift previous new to old
    swap(old_po, po);
    swap(old_clr, clr);

    //new colors = old color + succ sym old colors
    expanded_clrs.clear();
    // cerr << "--------" << endl;
    for (auto const s : ba.states()) {
      expanded_clrs[s] = {old_clr.at(s), {}}; //<C^{i-1}(s),...
      // cerr << s << " -> " << expanded_clrs.at(s).first << " | ";
      //... N^{i-1}(s)>
      for (auto const x : ba.state_outsyms(s)) {
        // = for each symbol, calculate the maximal (i-1)-types
        auto const sucs = ba.succ(s,x);
        set<unsigned> const loclrs = sucs | ranges::view::transform([&](state_t q){ return old_clr[q]; });
        for (auto const t : sucs) {
          auto const ct = old_clr.at(t);

          //check that ct is maximal succ color wrt po
          //i.e. no successor for same letter has a color that dominates ct
          bool maximal = false;
          if (!map_has_key(old_po, ct)) {
            maximal = true;
          } else {
            vector<unsigned> tmp;
            ranges::set_intersection(old_po.at(ct), loclrs, ranges::back_inserter(tmp));
            if (tmp.empty() || tmp==vector<unsigned>{ct})
              maximal = true;
          }

          if (maximal && !contains(expanded_clrs.at(s).second[x], ct)) {
            expanded_clrs.at(s).second[x].emplace(ct);
            // cerr << "(" << (int)x << "," << ct << ") ";
          }
        }
      }
      // cerr << endl;
    }

    //normalize colorset with lex. ordering
    vector<interm_clr> tmp = ranges::view::values(expanded_clrs);
    ranges::action::unique(ranges::action::sort(tmp));
    map<interm_clr, unsigned> norm_clr;
    int i=0;
    for (auto const& c : tmp) {
      norm_clr[c] = i++;
    }
    // map states to new number-colors
    clr.clear();
    for (auto const s : ba.states()) {
      clr[s] = norm_clr.at(expanded_clrs.at(s));
    }

    //new partial order
    po.clear();
    for (auto const s1 : ba.states()) {
      auto const c1 = expanded_clrs.at(s1);
      for (auto const s2 : ba.states()) {
        auto const c2 = expanded_clrs.at(s2);

        if (po_leq(old_po, c2.first, c1.first) && dominates(old_po, c1.second, c2.second)) {
          po[clr.at(s2)].emplace(clr.at(s1));
        }
      }
    }

    //debug print test
    /*
    for (auto const s : ba.states()) {
      cerr << s << " -> " << clr[s] << endl;
    }
    for (auto const it : po) {
      cerr << it.first << " <= " << seq_to_str(it.second) << endl;
    }
    */
  }

  //construct quotient automaton: colors as states
  Aut<string> ret(true, ba.get_name(), ba.get_aps(), clr.at(ba.get_init()));
  if (ba.state_buchi_accepting(ba.get_init()))
      ret.set_pri(ret.get_init(), 0);
  ret.tag_to_str = default_printer<string>();
  for (auto const s : ba.states()) {
    auto const c = clr.at(s);
    if (!ret.has_state(c)) {
      ret.add_state(c);
      if (ba.state_buchi_accepting(s))
        ret.set_pri(c, 0);
      // else
      //   ret.set_pri(c, 1);
    }
  }
  //edges between colors: (C(p), sym, C(q)) <=> (sym,C(q)) in N^i(p)
  for (auto const s : ba.states()) {
    for (auto const& it : expanded_clrs.at(s).second) {
      for (auto const trgclr : it.second) {
        // cerr << "edge " << clr.at(s) << " - " << (int) it.first << " > " << trgclr << endl;
        if (!ret.has_edge(clr.at(s), it.first, trgclr))
          ret.add_edge(clr.at(s), it.first, trgclr);
      }
    }
  }

  //decorate with info about original states in tag
  map<unsigned, set<state_t>> tag;
  for (auto const& it : clr) {
    tag[it.second].emplace(it.first);
  }
  for (auto const c : ret.states()) {
    ret.tag.put(seq_to_str(tag.at(c)), c);
  }
  //normalize and remap the PO
  ba_trim(ret); //remove now useless states
  auto const m = ret.normalize(); //get mapping of remaining state names (colors)
  //transform the partial order to match new ids -> gives detected language inclusions
  po_type norm_po;
  for (auto const it : po) {
    if (!map_has_key(m, it.first))
      continue; //a removed state
    norm_po[m.at(it.first)] =  it.second | ranges::view::filter([&m](auto v) { return map_has_key(m, v); })
                                         | ranges::view::transform([&m](auto v){ return m.at(v); });
  }

  return make_pair(move(ret), move(norm_po));
}

//context state set -> relatively acc, rej subsets
using Context = unordered_map<nba_bitset, pair<nba_bitset, nba_bitset>>;

Context get_context(auto const& aut, adj_mat const& mat, nba_bitset const asinks, map<unsigned, nba_bitset> const& impl,
                    shared_ptr<spdlog::logger> log = nullptr) {

  auto const psp = bench(log,"powerset_product",
                         WRAP(powerset_product(aut, mat, asinks, impl)));
  auto const psp_scci = get_sccs(psp.states(), aut_succ(psp));
  auto const psp_sccAcc = ba_scc_classify_acc(psp, psp_scci);
  // print_aut(psp);

  if (log)
    log->info("#states in 2^AxA: {}, #SCCs in 2^AxA: {}",
              psp.num_states(), psp_scci.sccs.size());

  Context ret;
  for (auto const s : psp.states()) {
    auto const& t = psp.tag.geti(s);
    int val = psp_sccAcc.at(psp_scci.scc_of.at(s));
    if (!map_has_key(ret, t.first))
      ret[t.first] = {};
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

