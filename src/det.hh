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
#include "common/hitset.hh"
#include "common/maxsat.hh"
#include "aut.hh"

namespace nbautils {

using PA = Aut<DetState>;

// takes: reference successor, mask for restricting candidates, valid pointer to sub-trie node,
// prefix up to sub-trie node, the prefix length and current depth
// returns: suitable candidate(s)
vector<DetState*> trie_dfs(DetState const& ref, auto const &msk,
                   auto const nodptr, nba_bitset const pref, int const i,
                   bool getAll=false) {
  if (!nodptr)
    return {};

  if ((nodptr->key & msk.first) != 0)
    return {};
  if (i<(int)msk.second.size() && (msk.second[i] & ~pref) != 0)
    return {};

  //try children in trie
  vector<DetState*> ret;
  for (auto const &sucnod : nodptr->suc) {
    auto cret = trie_dfs(ref, msk, sucnod.second.get(), pref|sucnod.first, i+1, getAll);
    ret.insert(end(ret),begin(cret),end(cret));
    if (!getAll && !ret.empty())
      return ret;
  }

  //here is a possible candidate state. need to check that all states that should
  //move down are actually moved down and that tuple order is weakly preserved
  if (nodptr->value && (msk.second.back() & ~pref) == 0) {
    DetState* cand = nodptr->value.get();

    if (cand && ref.tuples_finer_or_equal(*cand))
      ret.push_back(cand);
  }

  return ret;
}

vector<DetState*> existing_succ(DetConf const& dc, trie_map<nba_bitset, DetState>& existing,
    DetState const& cur, sym_t i, bool getAll=false) {
  // use Mueller/Schupp update for reference successor in trie query
  // TODO: make some override without replicating the whole thing
  auto dc2 = dc;
  dc2.update = UpdateMode::MUELLERSCHUPP;

  // Get MullerSchupp successor to span largest trie subtree possible
  DetState refsuc;
  pri_t refpri;
  tie(refsuc, refpri) = cur.succ(dc2, i);
  auto const ev = prio_to_event(refpri); //get dominant rank event
  auto const th = refsuc.to_tree_history(); //get dual structure

  //calculate k-equivalence level
  int k = ev.first + 2; //+1 for prepended powerset, +1 because first rank is 0
  if (k > (int)th.size()) //may happen due to breakpoint component becoming empty
    k = th.size();

  auto tht = th;
  tht.resize(k); //keep ranks 0 to k -> path to maximal collapsed k-equiv state in trie
  auto const msk = kcut_mask(th, k-1);
  auto const ini = existing.traverse(tht);
  vector<DetState*> cands;

  if (ini) //if the corresponding trie subtree exists, search for successors
    cands = trie_dfs(refsuc, msk, ini, tht.back(), 0, getAll);
  return cands;
}

// BFS-based determinization with supplied level update config
PA determinize(auto const& nba, DetConf const& dc, nba_bitset const& startset,
    auto const& pred, map<state_t, nba_bitset>* backmap = nullptr,
    map<state_t,map<sym_t, vector<state_t>>>* altmap = nullptr) {
  assert(nba.is_buchi());
  // create automaton with same letters etc
  state_t const myinit = 0;
  auto pa = PA(false, nba.get_name(), nba.get_aps(), myinit);
  pa.set_patype(PAType::MIN_EVEN);
  pa.tag_to_str = default_printer<DetState>();
  pa.tag.put(DetState(dc, startset), myinit); // initial state tag

  trie_map<nba_bitset, DetState> existing; //existing states organized in trie
  existing.put(pa.tag.geti(myinit).to_tree_history(), pa.tag.geti(myinit));
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

      if (dc.opt_suc || dc.hitset) {
        // if we try to reuse states during construction (smart successor selection)
        if (dc.opt_suc) {
          //check whether there is a (k-1) equivalent successor already
          //and replace if some suitable is found
          vector<DetState*> cands = existing_succ(dc, existing, cur, i);
          DetState* cand = cands.empty() ? nullptr : cands.front();
          if (cand) {
            // if we found a suitable successor in trie, use that
            /*
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
            */
            suclevel = *cand;
          }
        }
        //insert the corresponding successor if we use smart succ selection or postproc trie opt
        existing.put(suclevel.to_tree_history(), suclevel);

        //TODO: find bug
        /*
        auto tmp = existing.traverse(suclevel.to_tree_history());
        auto tmp2 = existing_succ(dc, existing, cur, i);
        assert(tmp != nullptr);
        assert(tmp->value != nullptr);
        assert(!tmp2.empty());
        */
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

  // cerr << "In trie: " << existing.size() << endl;
  // collect alternative edge targets from trie
  if (altmap) {
    auto& alts = *altmap;
    for (state_t const st : pa.states()) {
      DetState cur = pa.tag.geti(st);
      alts[st] = {};
      for (sym_t const i : pa.state_outsyms(st)) {
        alts[st][i] = pa.succ(st,i); //we definitely have the assigned succ.
        // and we possibly could have chosen another existing successor state
        for (DetState* const pst : existing_succ(dc, existing, cur, i, true)) {
          alts[st][i].push_back(pa.tag.get(*pst));
        }
        vec_to_set(alts[st][i]);
      }
    }
  }

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
int get_min_term_scc(auto const& succfun, SCCDat const& pai) {
    int mintermscc = -1;
    int mintermsz = -1;
    for (auto const it : pai.sccs) {
      auto const& sccnum = it.first;
      auto const ssccsz = it.second.size();
      auto const sucsccs = succ_sccs(succfun, pai, sccnum);

      if (sucsccs.empty() && (mintermsz<0 || ((int)ssccsz) < mintermsz)) {
        mintermscc = sccnum;
        mintermsz  = ssccsz;
      }
    }
    assert(mintermscc != -1);
    return mintermscc;
}

// remove useless mappings, i.e. restrict to states in keep set
void restrict_altmap(map<state_t, map<sym_t, vector<state_t>>>& altmap, set<state_t> const& keep) {
  // first remove outgoing
  auto altit = altmap.begin();
  while (altit != end(altmap)) {
    auto const tmp = altit++;
    //state removed -> remove constraints
    if (keep.find(tmp->first) == end(keep))
      altmap.erase(tmp);
  }
  // next remove incoming
  for (auto& altit : altmap) {
    for (auto& symtoes : altit.second) {
      vector<state_t> tmp;
      // assert(!symtoes.second.empty());
      set_intersection(begin(symtoes.second),end(symtoes.second),
                        begin(keep),end(keep),back_inserter(tmp));
      // assert(!tmp.empty());
      symtoes.second.swap(tmp);
    }
  }
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
    auto altmap = map<state_t, map<sym_t, vector<state_t>>>();
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
      }, &backmap, dc.hitset ? &altmap : nullptr);

    // for (auto const it : backmap) {
    //   cerr << pretty_bitset(sccpa.tag.geti(it.first).powerset) << " -> "
    //        << pretty_bitset(it.second) << endl;
    // }

    auto const sccpai = get_sccs(sccpa.states(), aut_succ(sccpa));

    //get states that belong to bottom SCC containing current powerset SCC rep
    auto const sccpa_succ = aut_succ(sccpa);
    auto mintermscc = get_min_term_scc(sccpa_succ, sccpai);
    auto sccstates = sccpai.sccs.at(mintermscc);
    vec_to_set(sccstates);

    //trim to that bottom SCC
    sccpa.remove_states(set_diff(sccpa.states() | ranges::to_vector, sccstates));
    // restrict_altmap(altmap, set<state_t>(begin(sccstates),end(sccstates)));

    if (dc.hitset) {
      //also trim constraint map of alternative edge targets
      vector<state_t> hitset; //holds hitset in greedy order
      set<state_t> sccsts(sccstates.begin(), sccstates.end()); //holds hitset as set
      size_t szbefore = sccsts.size();
      int hitsetround = 0;

      // cerr << "processing bottom SCC: " << sccsts.size();

      size_t oldsz = 0;
      while (oldsz != sccsts.size()) {
        hitsetround++;
        oldsz = sccsts.size();

        // remove useless mappings
        restrict_altmap(altmap, sccsts);

        /*
        // DEBUG OUTPUT
        print_aut(sccpa);
        cerr << "curr sccsts: " << seq_to_str(sccsts) << endl;
        cerr << "curr altmap:" << endl;
        for (auto const& xxx : altmap)
          for (auto const& yyy : xxx.second)
            cerr << xxx.first << "," << yyy.first << " -> " << seq_to_str(yyy.second) << endl;
        */

        // ----
        // GREEDY HITSET CALCULATION (can profit from iteration)
        /*
        // collect constraints
        map<state_t,set<state_t>> constraints;
        for (auto const& altmaps : altmap) {
          for (auto const& symtoes : altmaps.second) {
            constraints[constraints.size()] = set<state_t>(begin(symtoes.second), end(symtoes.second));
          }
        }

        // get a hitset
        // cerr << "before hitset: " << sccsts.size() << " - " << seq_to_str(sccsts) << endl;
        hitset = greedy_hitting_set(constraints);
        */
        // ----

        //TODO: add another flag for this variant, make this neat
        hitset = altmap_to_maxsat(altmap);
        // cerr << "calculated hitset: " << seq_to_str(hitset) << endl;

        // ----
        sccsts = set<state_t>(begin(hitset),end(hitset));

        restrict_altmap(altmap, sccsts);
        // cerr << "after hitset: " << sccsts.size() << " - " << seq_to_str(sccsts) << endl;

        // shrink hitset by computing another bottom SCC in rest
        auto const hitset_succ = [&](state_t s){
          set<state_t> sucs;
          for (auto const symtoes : altmap.at(s)) {
            set_intersection(begin(symtoes.second),end(symtoes.second),
                             begin(sccsts),end(sccsts),inserter(sucs,sucs.begin()));
          }
          return sucs;
        };
        auto const hitsetsccs = get_sccs(ranges::view::keys(altmap), hitset_succ);
        auto const minterm = get_min_term_scc(hitset_succ, hitsetsccs);
        sccsts = set<state_t>(begin(hitsetsccs.sccs.at(minterm)),end(hitsetsccs.sccs.at(minterm)));
        // cerr << "after botscc: " << sccsts.size() << " - " << seq_to_str(sccsts) << endl;
      }

      if (szbefore != sccsts.size()) {
        cerr << "performed " << hitsetround << " hitset rounds" << endl;
        cerr << "hitset state reduction: " << szbefore << " to " << sccsts.size() << endl;
      }

      // now altmap also only contains the hitset successors, we can redirect edges and remove useless
      int redirected = 0;
      for (auto const& altit : altmap) {
        //for all kept states
        state_t const st = altit.first;
        for (auto const& symtoes : altit.second) {
          sym_t const sym = symtoes.first;
          //remap edges to other alternative state, preserving priority
          auto const tmp = sccpa.succ_edges_raw(st, sym).cbegin();

          // assert(sccpa.succ(st,sym).size()==1);

          //old edge target and its priority (latter must be preserved)
          state_t const trg = tmp->first;
          pri_t const pri = tmp->second;

          // if (find(begin(symtoes.second),end(symtoes.second),trg) == end(symtoes.second)) {

            //pick some valid successor from hitset
            state_t ntrg = *(symtoes.second.cbegin());
            //improvement:
            //take earliest in hitset (as they are sorted by how "good" they are)
            //-> makes most happy, makes more paths "same", can improve Hopcroft minimization later
            for (state_t const hst : hitset) {
              if (sorted_contains(symtoes.second, hst)) {
                ntrg = hst;
                break;
              }
            }

            if (trg != ntrg) { //count edges that were effectively redirected
              // cerr << "redirecting " << st << "," << sym << " from " << trg << " to " << ntrg << endl;
              sccpa.remove_edge(st, sym, trg);
              sccpa.add_edge(st, sym, ntrg, pri);
              redirected++;
            }
            // assert(sccpa.succ(st,sym).size()==1);

          // }
        }
      }
      if (redirected > 0)
        cerr << "redirected " << redirected << " edges" << endl;

      //remove useless - again calculate a minimal bottom SCC after redirection and trim
      auto const sccpai2 = get_sccs(sccpa.states(), aut_succ(sccpa));
      mintermscc = get_min_term_scc(sccpa_succ, sccpai2);
      sccstates = sccpai2.sccs.at(mintermscc);
      vec_to_set(sccstates);

      // sccstates.clear();
      // copy(begin(sccsts),end(sccsts),back_inserter(sccstates));
      // vec_to_set(sccstates);
      // cerr << "sccstates size: " << sccstates.size() << endl;
      auto const tokill = set_diff(sccpa.states() | ranges::to_vector, sccstates);
      sccpa.remove_states(tokill);
      // cerr << "killed: " << seq_to_str(tokill) << endl;
      // cerr << "sccpa size after remove: " << sccpa.num_states() << endl;
      // cerr << "initial: " << sccpa.get_init() << endl;

      // cerr << "Keeping: " << sccsts.size() << endl;
    // cerr << "SCC states after norm: " << sccpa.num_states() << endl;
    // cerr << "altmap states after norm: " << altmap.size() << endl;
    // for (auto it : altmap) {
    //   for (auto jt : it.second)
    //     if (jt.second.size()>1)
    //       cerr << it.first << "," << jt.first << " -> alts: " << jt.second.size() << endl;
    // }
    }

    //normalize and insert into result automaton
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
    // cerr << "repst before: " << repst << endl;

    state_t const entry = psa.tag.get(origps.at(sccpa.get_init()));
    if (entry != rep) {
      vector<state_t> const path = find_path_from_to(psa, entry, rep);
      vector<sym_t> const word = get_word_from_path(psa, path);
      for (auto const x : word) {
        repst = sccpa.succ(repst, x).front();
      }
    }

    // cerr << "repst after: " << repst << endl;

    // int repst=-1;
    //     for (auto const st : sccpa.states()) {
    //       if (sccpa.tag.hasi(st) && origps.at(st) == repps) {
    //         repst = st;
    //         break;
    //       }
    // }
    assert(repst>=0);
    ps2pa[rep] = repst;

    // assert(sccpa.is_deterministic());

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

            // cerr << pst << ", " << sym  << endl; //<< " -> " << seq_to_str(psucst) << endl;

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
