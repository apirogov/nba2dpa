#include "aut.hh"
#include "det.hh"

namespace nbautils {
  using namespace std;


/*
//input: automaton with sccinfo, a powerset that exists in the automaton
//output: scc number of smallest terminal automaton containing that powerset
int get_min_term_scc_with_powerset(PA const& pa, SCCDat<state_t> const& pai, vector<small_state_t> const& s) {
    // get all levels with same powerset (take all, remove non-copies)
    auto copies = pa.states();
    copies.erase(remove_if(begin(copies), end(copies),
                 [&pa,&s](state_t st){ return pa.tag->geti(st).states() != s; }), end(copies));

    assert(!copies.empty()); // provided powerset contained in automaton

    //keep one copy per SCC (we need to check each candidate SCC only once)
    sort(begin(copies), end(copies), [&](state_t s, state_t t){
        return pai.scc_of.at(s) < pai.scc_of.at(t); });
    auto it = unique(begin(copies), end(copies), [&](state_t s, state_t t){
        return pai.scc_of.at(s) == pai.scc_of.at(t); });
    copies.erase(it, end(copies));

#ifndef NDEBUG
    bool found = false;
#endif

    //find smallest bottom SCC (bottom ensures that all powersets in PS SCC are reachable)
    int mintermscc = 0;
    size_t mintermsz = pa.num_states()+1;

    for (auto const s : copies) {
      auto const sccnum = pai.scc_of.at(s);
      auto const ssccsz = pai.sccs.at(sccnum).size();

      // cerr << "in scc " << sccnum << " (size " << ssccsz << ")"<< " with succ sccs: ";

      succ_fun<state_t> pa_sucs = [&pa](state_t s){ return pa.succ(s); };
      auto const sucsccs = succ_sccs(pai, sccnum, pa_sucs);

      // for (auto v : sucsccs)
      //   cout << (int)v << ",";
      // cout << endl;

      if (sucsccs.empty() && ssccsz < mintermsz) {
#ifndef NDEBUG
        found = true;
#endif
        mintermscc = sccnum;
        mintermsz  = ssccsz;
      }
    }

#ifndef NDEBUG
    assert(found); //there must be a bottom SCC with this powerset by construction
#endif
    return mintermscc;
}

// determinization of each powerset component separately, then fusing
PA::uptr determinize(LevelConfig const& lc, PS const& psa, SCCDat<state_t> const& psai) {
  map<state_t, state_t> ps2pa;
  auto ret = std::make_unique<PA>(Acceptance::PARITY, lc.aut->get_name(), lc.aut->get_aps());
  ret->set_patype(PAType::MIN_EVEN);
  ret->tag_to_str = [](Level const& l){ return l.to_string(); };

  // for (auto &it : psai.sccrep) {
  unsigned i=psai.sccs.size()-1;
  for (auto it=psai.sccs.crbegin(); it!=psai.sccs.crend(); ++it) {
    auto const& scc = i;
    auto const& rep = it->front();

    auto const repps = psa.tag->geti(rep); //powerset of scc representative
    if (repps.empty())
      continue;

    // cerr << "determinize called on SCC with statenumber " << psai.sccs.at(i).size() << endl;

    auto sccpa = determinize(lc, repps, [&psa,&psai,&scc](Level const& l){
        auto const pset = l.states();
        if (!psa.tag->has(pset))
          throw runtime_error("we reached a weird successor!");

        auto s = psa.tag->get(pset);
        //don't explore levels with powerset in other scc
        if (psai.scc_of.at(s) != scc)
          return false;
        return true;
        });

    auto const sccpai = get_sccs(sccpa->states(), swa_succ(*sccpa));

    auto const mintermscc = get_min_term_scc_with_powerset(*sccpa, *sccpai, repps);
    auto sccstates = sccpai->sccs.at(mintermscc);
    vec_to_set(sccstates);
    sccpa->remove_states(set_diff(sccpa->states(), sccstates));
    sccpa->normalize(ret->num_states());

    //find representative in trimmed SCC PA graph
    int repst=-1;
    for (auto const st : sccpa->states()) {
      if (sccpa->tag->hasi(st) && sccpa->tag->geti(st).states() == repps) {
        repst = st;
        break;
      }
    }
    assert(repst >= (int)ret->num_states());

    //copy the SCC
    ret->insert(*sccpa);

    //update PS State -> PA State inter-SCC map
    //by exploring the scc of PS and simulating it in PA
    ps2pa[rep] = repst;
    bfs(rep, [&](auto const& st, auto const& visit, auto const&) {
        auto const pst = ps2pa.at(st);
        //add successor powerset states in same scc
        for (auto sym : psa.outsyms(st)) {
          for (auto sucst : psa.succ(st, sym)) {
            if (map_has_key(ps2pa,sucst) || psai.scc_of.at(sucst) != psai.scc_of.at(st))
              continue;

            auto psucst = sccpa->succ(pst, sym);
            assert(psucst.size()==1); //PA SCC subautomaton must be deterministic!
            ps2pa[sucst] = psucst.front();

            visit(sucst);
          }
        }
    });
    // cerr << "after norm " << sccpa->num_states() << endl;
    // cout << mintermscc << " with " << sccpa->num_states() << endl;
    --i;
  }

#ifndef NDEBUG
  set<state_t> mapvals;
  for (auto it : ps2pa)
    mapvals.emplace(it.second);
#endif
  assert(mapvals.size() == map_get_keys(ps2pa).size()); //i, "not injective!");

  //set initial state to mapped initial state
  ret->set_init({ps2pa.at(psa.get_init().front())});

  //traverse resulting DPA and add missing edges (which are present in PSA)
  bfs(ret->get_init().front(), [&](auto const& st, auto const& visit, auto const&) {
      for (auto i=0; i<ret->num_syms(); i++) {
        if (!ret->state_has_outsym(st,i)) { //missing successor candidate
          //get corresponding state in PSA
          auto const pst = psa.tag->get(ret->tag->geti(st).states());
          //check out its successors
          auto const psucs = psa.succ(pst, i);

          if (!psucs.empty()) { //indeed missing successor!
            assert(psucs.size()==1); //PS subautomaton is deterministic
            // set PA rep of PS successor as PA successor
            ret->set_succs(st, i, {ps2pa.at(psucs.front())});
          }
        }
        // if now a successor is present, visit it too
        for (auto const sucst : ret->succ(st,i))
          visit(sucst);
      }
    });

  return move(ret);
}
*/

}  // namespace nbautils
