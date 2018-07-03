#pragma once

#include <functional>
#include <queue>
#include <set>
#include <vector>
#include "detstate.hh"
#include "common/scc.hh"
#include "aut.hh"

namespace nbautils {

using PA = Aut<DetState>;

// BFS-based determinization with supplied level update config
PA determinize(auto const& nba, DetConf const& dc, nba_bitset const& startset, auto const& pred) {
  assert(nba.is_buchi());
  // create automaton with same letters etc
  state_t const myinit = 0;
  auto pa = PA(false, nba.get_name(), nba.get_aps(), myinit);
  pa.set_patype(PAType::MIN_EVEN);
  // pa.tag_to_str = [](ostream& out, auto const& ds){ out << ds; };
  pa.tag_to_str = default_printer<DetState>();
  pa.tag.put(DetState(dc, startset), myinit); // initial state tag

  int numvis=0;
  bfs(myinit, [&](auto const& st, auto const& visit, auto const&) {
    // get inner states of current ps state
    auto const cur = pa.tag.geti(st);

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

      if (!pred(suclevel)) //predicate not satisfied -> don't explore this node
        continue;

      //check whether there is already a state in the graph with this label
      auto const sucst = pa.tag.put_or_get(suclevel, pa.num_states());

      //if this is a new successor, add it to graph and enqueue it:
      if (!pa.has_state(sucst))
        pa.add_state(sucst);
      // create edge
      pa.add_edge(st, i, sucst, sucpri);
      // schedule for bfs
      visit(sucst);
    }
  });

  // cerr << "determinized to " << pa->num_states() << " states" << endl;
  return pa;
}

// start with initial state of NBA, explore by DFS completely
PA determinize(auto const& nba, DetConf const& dc) {
  return determinize(nba, dc, nba_bitset(1<<nba.get_init()), const_true);
}

//input: automaton with sccinfo, a powerset that exists in the automaton
//output: scc number of smallest terminal automaton containing that powerset
int get_min_term_scc_with_powerset(PA const& pa, SCCDat const& pai, nba_bitset const& s) {
    // get all levels with same powerset (take all, remove non-copies)
    auto copies = pa.states() | ranges::to_vector;
    copies = std::move(copies)
           | ranges::action::remove_if([&pa,&s](state_t st){ return pa.tag.geti(st).powerset != s; });

    assert(!copies.empty()); // provided powerset should be contained in automaton

    //keep one copy per SCC (we need to check each candidate SCC only once)
    copies = std::move(copies)
               | ranges::action::sort(  [&](auto s, auto t){ return pai.scc_of.at(s) <  pai.scc_of.at(t); })
               | ranges::action::unique([&](auto s, auto t){ return pai.scc_of.at(s) == pai.scc_of.at(t); });

    bool found = false;

    //find smallest bottom SCC (bottom ensures that all powersets in PS SCC are reachable)
    int mintermscc = 0;
    size_t mintermsz = pa.num_states()+1;
    for (auto const s : copies) {
      auto const sccnum = pai.scc_of.at(s);
      auto const ssccsz = pai.sccs.at(sccnum).size();
      // cerr << "in scc " << sccnum << " (size " << ssccsz << ")"<< " with succ sccs: ";
      auto const sucsccs = succ_sccs(aut_succ(pa), pai, sccnum);
      // for (auto v : sucsccs)
      //   cout << (int)v << ",";
      // cout << endl;

      if (sucsccs.empty() && ssccsz < mintermsz) {
        found = true;
        mintermscc = sccnum;
        mintermsz  = ssccsz;
      }
    }

    assert(found); //there must be a bottom SCC with this powerset by construction
    return mintermscc;
}

// determinization of each powerset component separately, then fusing
PA determinize(auto const& nba, DetConf const& dc, PS const& psa, SCCDat const& psai) {
  map<state_t, state_t> ps2pa;
  PA ret(false, nba.get_name(), nba.get_aps(), 0);
  ret.set_patype(PAType::MIN_EVEN);
  ret.tag_to_str = default_printer<DetState>();

  // for (auto &it : psai.sccrep) {
  for (auto it : ranges::view::all(psai.sccs) | ranges::view::reverse) {
  // for (auto it=psai.sccs.crbegin(); it!=psai.sccs.crend(); ++it) {
    auto const& scc = it.first;
    auto const& rep = it.second.front();

    nba_bitset const repps = psa.tag.geti(rep); //powerset of scc representative
    if (repps == 0) //empty powerset
      continue;

    // cerr << "determinize called on SCC with statenumber " << psai.sccs.at(i).size() << endl;

    auto sccpa = determinize(nba, dc, repps, [&psa,&psai,&scc](DetState const& ds){
        if (!psa.tag.has(ds.powerset))
          throw runtime_error("we reached a weird successor!");
        auto const s = psa.tag.get(ds.powerset);
        //don't explore levels with powerset in other scc
        if (psai.scc_of.at(s) != scc)
          return false;
        return true;
      });
    auto const sccpai = get_sccs(sccpa.states() | ranges::to_vector, aut_succ(sccpa));

    //get states that belong to bottom SCC containing current powerset SCC rep
    auto const mintermscc = get_min_term_scc_with_powerset(sccpa, sccpai, repps);
    auto sccstates = sccpai.sccs.at(mintermscc);
    vec_to_set(sccstates);

    //trim to that bottom SCC, normalize and insert into result automaton
    sccpa.remove_states(set_diff(sccpa.states() | ranges::to_vector, sccstates));
    sccpa.normalize(ret.num_states());
    ret.insert(sccpa);

    //find representative in trimmed SCC PA graph (which is start for exploration)
    int repst=-1;
        for (auto const st : sccpa.states()) {
          if (sccpa.tag.hasi(st) && sccpa.tag.geti(st).powerset == repps) {
            repst = st;
            break;
          }
    }
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
      for (auto i : ret.syms()) {
        if (!ret.state_has_outsym(st, i)) { //missing successor candidate
          //get corresponding state in PSA
          auto const pst = psa.tag.get(ret.tag.geti(st).powerset);
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
