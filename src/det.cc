#include "det.hh"
#include "types.hh"

namespace nbautils {
  using namespace std;

// BFS-based determinization with supplied level update config
PA::uptr determinize(LevelConfig const& lc, vector<small_state_t> const& startset, auto const& pred) {
  // same letters etc
  auto pa = std::make_unique<PA>(lc.aut->get_name(), lc.aut->get_aps());

  state_t myinit = 0;
  pa->set_init({myinit});
  pa->set_accs(myinit, {0}); // initial state priority does not matter
  pa->tag->put(Level(lc, startset), myinit); // initial state tag

  bfs(myinit, [&](auto const& st, auto const& visit) {
    // get inner states of current ps state
    auto const curlevel = pa->tag->geti(st);

    for (auto i = 0; i < pa->num_syms(); i++) {
      // calculate successor level
      auto const suclevel = curlevel.succ(lc, i);

      if (suclevel.powerset == 0) //is an empty set -> invalid successor
        continue;

      if (!pred(suclevel)) //predicate not satisfied -> don't explore this node
        continue;

      //check whether there is already a state in the graph with this label
      auto const sucst = pa->tag->put_or_get(suclevel, pa->num_states());

      //if this is a new successor, add it to graph and enqueue it:
      if (!pa->has_state(sucst))
        pa->add_state(sucst);
      // assign priority according to resulting level
      if (!pa->has_accs(sucst))
        pa->set_accs(sucst,{suclevel.prio});
      // create edge
      pa->adj[st][i].push_back(sucst);
      // schedule for bfs
      visit(sucst);
    }
  });

  return move(pa);
}

// start with initial state of NBA, explore by DFS completely
PA::uptr determinize(LevelConfig const& lu) {
  auto const explore_all = [](auto const&){ return true; }; //trivial predicate
  //start with the set of initial states of provided automaton
  auto const startset = vector<small_state_t>(cbegin(lu.aut->get_init()), cend(lu.aut->get_init()));

  return move(determinize(lu,startset,explore_all));
}

//input: automaton with sccinfo, a powerset that exists in the automaton
//output: scc number of smallest terminal automaton containing that powerset
scc_t get_min_term_scc_with_powerset(PA const& pa, SCCInfo const& pai, vector<small_state_t> const& s) {
    // get all levels with same powerset (take all, remove non-copies)
    auto copies = pa.states();
    copies.erase(remove_if(begin(copies), end(copies),
                 [&pa,&s](state_t st){ return pa.tag->geti(st).states() != s; }), end(copies));

    assert(!copies.empty()); // provided powerset contained in automaton

    // cout << "num copies: " << copies.size() << endl;

    bool found = false;

    scc_t mintermscc = 0;
    size_t mintermsz = pa.num_states()+1;

    for (auto const s : copies) {
      auto const sccnum = pai.scc.at(s);
      auto const ssccsz = pai.sccsz.at(sccnum);

      // cout << "in scc " << sccnum << " (size " << ssccsz << ")"<< " with succ sccs: ";

      auto const sucsccs = succ_sccs(pa, pai, sccnum);

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
PA::uptr determinize(LevelConfig const& lc, PS<Acceptance::BUCHI> const& psa, SCCInfo const& psai) {
  map<state_t, state_t> ps2pa;
  int curnumstates = 0;
  auto ret = std::make_unique<PA>(lc.aut->get_name(), lc.aut->get_aps());
  for (auto &it : psai.sccrep) {
    auto const& scc = it.first;
    auto const& rep = it.second;
    auto const repps = psa.tag->geti(rep); //powerset of scc representative

    // cout << "scc " << scc <<" (" << psai.sccsz.at(scc) << " states) -> ";

    auto sccpa = determinize(lc, repps, [&psa,&psai,&scc](Level const& l){
        auto const pset = l.states();
        if (!psa.tag->has(pset))
          throw runtime_error("we reached a weird successor!");

        auto s = psa.tag->get(pset);
        //don't explore levels with powerset in other scc
        if (psai.scc.at(s) != scc)
          return false;
        return true;
        });

    auto const sccpai = get_scc_info(*sccpa);

    auto const mintermscc = get_min_term_scc_with_powerset(*sccpa, *sccpai, repps);
    auto const sccstates = scc_states(*sccpa, *sccpai, mintermscc);
    sccpa->remove_states(set_diff(sccpa->states(), sccstates));
    cout << "before norm " << sccpa->num_states() << endl;
    sccpa->normalize(curnumstates);
    cout << "after norm " << sccpa->num_states() << endl;
    cout << seq_to_str(sccpa->states(),",") << endl;
    ret->insert(*sccpa);
    // cout << "after insert " << ret->num_states() << endl;

    //update PS State -> PA State inter-SCC map
    //by exploring the scc of PS and simulating it in PA
    ps2pa[rep] = curnumstates;
    bfs(rep, [&](auto const& st, auto const& visit) {
        auto const pst = ps2pa.at(st);
        //add successor powerset states in same scc
        for (auto sym : psa.outsyms(st)) {
          for (auto sucst : psa.succ(st, sym)) {
            if (map_has_key(ps2pa,sucst) || psai.scc.at(sucst) != psai.scc.at(st))
              continue;

            auto psucst = sccpa->succ(pst, sym);
            assert(psucst.size()==1); //PA SCC subautomaton must be deterministic!
            ps2pa[sucst] = psucst.front();
            cout << (int)sucst << " -> " << psucst.front() << endl;
            if (psucst.front()==84) {
            for (auto s : psa.tag->geti(sucst))
            cout << (int)s << ",";
            cout << endl;
            cout << ret->tag->geti(psucst.front()) << endl;
            }
            visit(sucst);
          }
        }
    });
    // cout << "after norm " << sccpa->num_states() << endl;
    // cout << mintermscc << " with " << sccpa->num_states() << endl;
    // curnumstates += sccpa->num_states();
    curnumstates = ret->num_states();
  }

#ifndef NDEBUG
  set<state_t> mapvals;
  for (auto it : ps2pa)
    mapvals.emplace(it.second);
#endif
  assert(mapvals.size() == map_get_keys(ps2pa).size()); //i, "not injective!");
  // cout << "total: " << curnumstates << endl;

  //set initial state to mapped initial state
  ret->set_init({ps2pa.at(psa.get_init().front())});

  //traverse PSA again to construct interSCC edges
  bfs(psa.get_init().front(), [&](auto const& st, auto const& visit) {
      cout << "visit " << st << endl;
      for (auto sym : psa.outsyms(st)) {
        for (auto sucst : psa.succ(st, sym)) {
          if (psai.scc.at(sucst) != psai.scc.at(st)) {
            //we can't have this edge yet in the PA so add it
            auto const past = ps2pa.at(st);
            auto const pasuc = ps2pa.at(sucst);
            cout << (int)past << " - " << (int)sym << " > " << (int)pasuc << endl;

            auto const old = ret->succ(past, sym);
            if (old.size() > 0)
              cout << "weird!" << endl;
            auto const tmp = set_merge(old, {pasuc});
            ret->adj.at(past)[sym] = tmp;
          }

          visit(sucst);
        }
      }
    });

  return move(ret);
}

}  // namespace nbautils
