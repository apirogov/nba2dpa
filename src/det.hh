#pragma once

#include <functional>
#include <queue>
#include <set>
#include <unordered_set>
#include <vector>
#include "graph.hh"
#include "detstate.hh"
#include "common/scc.hh"
#include "common/types.hh"
#include "common/trie_map.hh"
#include "aut.hh"

namespace nbautils {

using PA = Aut<DetState>;

// takes: reference successor, mask for restricting candidates, valid pointer to sub-trie node,
// prefix up to sub-trie node, the prefix length and current depth
// returns: suitable candidate
DetState* trie_dfs(DetState const& ref, auto const &msk,
                   auto const nodptr, nba_bitset const pref, int const i) {
  if (!nodptr)
    return nullptr;

  if ((nodptr->key & msk.first) != 0)
    return nullptr;
  if (i<(int)msk.second.size() && (msk.second[i] & ~pref) != 0)
    return nullptr;

  //try children in trie
  for (auto const &sucnod : nodptr->suc) {
    DetState *ret = trie_dfs(ref, msk, sucnod.second.get(), pref|sucnod.first, i+1);
    if (ret)
      return ret;
  }

  //here is a possible candidate state. need to check that all states that should
  //move down are actually moved down and that tuple order is weakly preserved
  if (nodptr->value && (msk.second.back() & ~pref) == 0) {
    DetState* cand = nodptr->value.get();
    if (!cand)
      return nullptr;

    if (ref.tuples_finer_or_equal(*cand))
      return cand;
  }

  //nothing suitable found
  return nullptr;
}

// BFS-based determinization with supplied level update config
PA determinize(auto const& nba, DetConf const& dc, nba_bitset const& startset,
    auto const& pred, map<state_t, nba_bitset>* backmap = nullptr) {
  assert(nba.is_buchi());
  // create automaton with same letters etc
  state_t const myinit = 0;
  auto pa = PA(false, nba.get_name(), nba.get_aps(), myinit);
  pa.set_patype(PAType::MIN_EVEN);
  pa.tag_to_str = default_printer<DetState>();
  pa.tag.put(DetState(dc, startset), myinit); // initial state tag

  trie_map<nba_bitset, DetState> existing; //existing states organized in trie
  // use Mueller/Schupp update for reference successor in trie query
  auto dc2 = dc;
  dc2.update = UpdateMode::MUELLERSCHUPP;
  // dc2.puretrees = false;

  int numvis=0;
  //always track normal successor powerset and det state in parallel
  if (backmap)
    (*backmap)[myinit] = startset;
  unordered_set<state_t> vis2nd;
  bfs(make_pair(startset, myinit), [&](auto const& stp, auto const& visit, auto const&) {
    // get inner states of current macro state
    auto const cur = pa.tag.geti(stp.second);

    if (contains(vis2nd,stp.second)) //only explore if the PA state unexplored
      return;
    vis2nd.emplace(stp.second);

    // cout << "visit " << curlevel.to_string() << endl;
    ++numvis;
    if (numvis % 5000 == 0) //progress indicator
      cerr << numvis << endl;

    for (auto const i : pa.syms()) {
      // calculate successor level
      DetState suclevel;
      pri_t sucpri;
      tie(suclevel, sucpri) = cur.succ(dc, i);
      // cout << "suc " << suclevel.to_string() << endl;

      if (suclevel.powerset == 0) //is an empty set -> invalid successor
        continue;

      //calculate powerset successor to track scc
      // cerr << "visiting " << pretty_bitset(stp.first) << endl;
      nba_bitset sucset = powersucc(dc.aut_mat, stp.first, i, dc.aut_asinks, dc.impl_mask);
      if (!pred(sucset)) //predicate not satisfied -> don't explore this node
        continue;

      //check whether there is a (k-1) equivalent successor already
      //and replace if some suitable is found
      if (dc.opt_suc) {
        auto const ev = prio_to_event(sucpri); //get dominant rank event

        // Get MullerSchupp successor to span largest trie subtree possible
        DetState refsuc;
        pri_t refpri;
        tie(refsuc, refpri) = cur.succ(dc2, i);
        auto const th = refsuc.to_tree_history(); //get dual structure

        //calculate k-equivalence level
        int k = ev.first + 2; //+1 for prepended powerset, +1 because first rank is 0
        if (k > (int)th.size()) //may happen due to breakpoint component becoming empty
          k = th.size();

        auto tht = th;
        tht.resize(k); //keep ranks 0 to k -> path to maximal collapsed k-equiv state in trie
        auto const msk = kcut_mask(th, k-1);
        auto const ini = existing.traverse(tht);
        DetState* cand = nullptr;
        if (ini) //if the corresponding trie subtree exists, search for successors
          cand = trie_dfs(refsuc, msk, ini, tht.back(), 0);
        if (cand) {
          // if we found a suitable successor in trie, use that
          if (suclevel != *cand && dc.debug) {
            cerr << "suc of:\t" << cur << " : " << ev.first << " " << (ev.second ? "+" : "-") << endl
                 << "replcd:\t" << suclevel << endl
                 << "via:\t" << refsuc << endl
                 << "with:\t" << *cand << endl;
            cerr << "node:\t" << tht << endl;
            cerr << "msk:\t" << pretty_bitset(msk.first) << " | " << msk.second << endl;
            cerr << "addr:\t" << suclevel.to_tree_history() << endl;
            cerr << "raddr:\t" << refsuc.to_tree_history() << endl;
            cerr << "naddr:\t" << cand->to_tree_history() << endl;
          }
          suclevel = *cand;
        } else {
          // otherwise, insert and use the one we calculated using user config
          existing.put(suclevel.to_tree_history(), suclevel);
        }
      }

      // now suclevel definitely has some suitable successor we decided on
      // ----

      //check whether there is already a state in the graph with this label
      auto const sucst = pa.tag.put_or_get(suclevel, pa.num_states());

      //if this is a new successor, add it to graph and enqueue it:
      if (!pa.has_state(sucst)) {
        pa.add_state(sucst);
        if (backmap) //assign a language equiv. original powerset from SCC to the state
          (*backmap)[sucst] = sucset;
      }
      // create edge
      pa.add_edge(stp.second, i, sucst, sucpri);
      // schedule for bfs
      visit(make_pair(sucset, sucst));
    }
  });

  // cerr << "determinized to " << pa->num_states() << " states" << endl;
  return pa;
}

// start with initial state of NBA, explore by DFS completely
PA determinize(auto const& nba, DetConf const& dc) {
  nba_bitset initset = 0;
  initset[nba.get_init()] = 1; //1<<x does not work as expected
  return determinize(nba, dc, initset, const_true);
}

//find smallest bottom SCC (bottom ensures that all powersets in PS SCC are reachable)
int get_min_term_scc(PA const& pa, SCCDat const& pai) {
    int mintermscc = -1;
    size_t mintermsz = pa.num_states()+1;
    for (auto const it : pai.sccs) {
      auto const& sccnum = it.first;
      auto const ssccsz = it.second.size();
      auto const sucsccs = succ_sccs(aut_succ(pa), pai, sccnum);

      if (sucsccs.empty() && ssccsz < mintermsz) {
        mintermscc = sccnum;
        mintermsz  = ssccsz;
      }
    }
    assert(mintermscc != -1);
    return mintermscc;
}

// determinization of each powerset component separately, then fusing
PA determinize(auto const& nba, DetConf const& dc, PS const& psa, SCCDat const& psai) {
  map<state_t, state_t> ps2pa;
  map<state_t, nba_bitset> origps;
  PA ret(false, nba.get_name(), nba.get_aps(), 0);
  ret.remove_states({0}); //we want a blank graph without states
  ret.set_patype(PAType::MIN_EVEN);
  ret.tag_to_str = default_printer<DetState>();

  for (auto it : ranges::view::all(psai.sccs) | ranges::view::reverse) {
    auto const& scc = it.first;
    auto const& rep = it.second.front();

    nba_bitset const repps = psa.tag.geti(rep); //powerset of scc representative
    if (repps == 0) //empty powerset
      continue;

    // cerr << "repps: " << pretty_bitset(repps) << endl;

    //this map will map tuples with weird optimizations to the powerset states they represent
    auto backmap = map<state_t, nba_bitset>();
    auto sccpa = determinize(nba, dc, repps, [&psa,&psai,&scc](nba_bitset const& ds){
        if (!psa.tag.has(ds)) {
          // cerr << "reached " << pretty_bitset(ds) << endl;
          throw runtime_error("we reached a weird successor!");
        }
        auto const s = psa.tag.get(ds);
        //don't explore levels with powerset in other scc
        if (psai.scc_of.at(s) != scc)
          return false;
        return true;
      }, &backmap);

    // for (auto const it : backmap) {
    //   cerr << pretty_bitset(sccpa.tag.geti(it.first).powerset) << " -> "
    //        << pretty_bitset(it.second) << endl;
    // }

    auto const sccpai = get_sccs(sccpa.states(), aut_succ(sccpa));

    //get states that belong to bottom SCC containing current powerset SCC rep
    auto const mintermscc = get_min_term_scc(sccpa, sccpai);
    auto sccstates = sccpai.sccs.at(mintermscc);
    vec_to_set(sccstates);

    //trim to that bottom SCC, normalize and insert into result automaton
    sccpa.remove_states(set_diff(sccpa.states() | ranges::to_vector, sccstates));
    auto const normmap = sccpa.normalize(ret.num_states());
    for (auto const st : sccstates) {
      //map detstate to the semantic/original powerset
      origps[normmap.at(st)] = backmap.at(st);
    }
    ret.insert(sccpa);

    // cerr << ret.num_states() << " " << origps.size() << endl;
    // cerr << "mintermscc: " << mintermscc << " , " << seq_to_str(sccstates) << " -- " << seq_to_str(sccpa.states()) << endl;
    // for (auto const st : sccstates) {
    //   cerr << st << " -> " << normmap.at(st) << " -> " << pretty_bitset(origps.at(normmap.at(st))) << endl;
    // }

    //find representative in trimmed SCC PA graph (which is start for exploration)
    // cerr <<  pretty_bitset(psa.tag.get(origps.at(sccpa.get_init()))) << "->" << pretty_bitset(repps) << endl;
    // TODO: do we even care where to start exploration?
    state_t repst = sccpa.get_init();
    state_t const entry = psa.tag.get(origps.at(sccpa.get_init()));
    if (entry != rep) {
      vector<state_t> const path = find_path_from_to(psa, entry, rep);
      vector<sym_t> const word = get_word_from_path(psa, path);
      for (auto const x : word) {
        repst = sccpa.succ(repst, x).front();
      }
    }

    // int repst=-1;
    //     for (auto const st : sccpa.states()) {
    //       if (sccpa.tag.hasi(st) && origps.at(st) == repps) {
    //         repst = st;
    //         break;
    //       }
    // }
    assert(repst>=0);
    ps2pa[rep] = repst;

    //update PS State -> PA State inter-SCC map
    //by exploring the scc of PS and simulating it in PA
    bfs(rep, [&](auto const& st, auto const& visit, auto const&) {
        auto const pst = ps2pa.at(st);
        for (auto sym : psa.state_outsyms(st)) {
          for (auto sucst : psa.succ(st, sym)) {
            //add successor powerset states in same scc
            if (map_has_key(ps2pa, sucst) || psai.scc_of.at(sucst) != psai.scc_of.at(st))
              continue;

            auto const psucst = sccpa.succ(pst, sym);
            assert(psucst.size() == 1); //PA SCC subautomaton must be deterministic!
            ps2pa[sucst] = psucst.front();

            visit(sucst);
          }
        }
    });
  }

  //set initial state to mapped initial state
  ret.set_init(ps2pa.at(psa.get_init()));

  //traverse resulting DPA and add missing edges (which are present in PSA)
  bfs(ret.get_init(), [&](auto const& st, auto const& visit, auto const&) {
      //get corresponding state in PSA
      // auto const pst = psa.tag.get(ret.tag.geti(st).powerset);
      auto const pst = psa.tag.get(origps.at(st));

      for (auto i : ret.syms()) {
        if (!ret.state_has_outsym(st, i)) { //missing successor candidate
          //check out its successors
          auto const psucs = psa.succ(pst, i);

          if (!psucs.empty()) { //indeed missing successor!
            assert(psucs.size()==1); //PS subautomaton is deterministic
            // set PA rep of PS successor as PA successor
            ret.add_edge(st, i, ps2pa.at(psucs.front()), 0);
          }
        }
        // if now a successor is present, will also be visited
        for (auto const sucst : ret.succ(st,i))
          visit(sucst);
      }
    });

  return ret;
}

}  // namespace nbautils
