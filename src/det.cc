#include "det.hh"
#include "types.hh"

namespace nbautils {

// BFS-based determinization with supplied level update config
PA::uptr determinize(LevelConfig const& lc) {
  auto& aut = *lc.aut;

  auto pap = std::make_unique<PA>(PA());
  auto& pa = *pap;

  pa.tag = std::make_unique<pa_tag_store>(pa_tag_store());
  auto& tag = pa.tag;

  // same letters etc
  pa.meta = aut.meta;
  // initial state is 0
  pa.init = 0;
  pa.add_state(0);
  pa.set_accs(pa.init, {0});  // priority does not matter

  // q_0 tag
  auto pinit = Level(lc, std::vector<small_state_t>{(small_state_t)aut.init});
  tag->put(pinit, pa.init);

  bfs(pa.init, [&](auto const& st, auto const& visit) {
    // get inner states of current ps state
    auto curlevel = tag->get(st);

    for (auto i = 0; i < (int)pa.num_syms(); i++) {
      auto suclevel = curlevel.succ(lc, i);
      if (suclevel.powerset == 0) //emptyset
        continue;
      state_t sucst = tag->put_or_get(suclevel, tag->size());
      if (!pa.has_state(sucst))
        pa.add_state(sucst);
      if (!pa.has_accs(sucst))  // assign priority according to resulting level
        pa.set_accs(sucst,{suclevel.prio});
      pa.adj[st][i].push_back(sucst);

      visit(sucst);
    }
  });

  return move(pap);
}

// determinization of each powerset component separately, then fusing
PA::uptr determinize(LevelConfig const& lc, PS<Acceptance::BUCHI> const& psa, SCCInfo const& psai) {
  auto& aut = *lc.aut;

  int total = 0;
  for (auto &it : psai.sccrep) {
    auto const& scc = it.first;
    auto const& rep = it.second;
    auto const& repps = psa.tag->get(rep);

    cout << "scc " << scc <<" (" << psai.sccsz.at(scc) << " states) -> ";

    auto pap = std::make_unique<PA>(PA());
    auto& pa = *pap;
    pa.tag = std::make_unique<pa_tag_store>(pa_tag_store());
    auto& tag = pa.tag;

    // same letters etc
    pa.meta = aut.meta;

    // initial state is 0
    pa.init = 0;
    pa.add_state(0);
    pa.set_accs(pa.init, {0});  // priority does not matter

    // initial tag - powerset of representative
    auto pinit = Level(lc, repps);
    tag->put(pinit, pa.init);

    bfs(pa.init, [&](auto const& st, auto const& visit) {
      // get inner states of current ps state
      auto curlevel = tag->get(st);

      for (auto i = 0; i < (int)pa.num_syms(); i++) {
        auto suclevel = curlevel.succ(lc, i);
        if (suclevel.powerset == 0) //emptyset
          continue;

        auto pset = suclevel.states();
        if (psa.tag->has(pset)) {
          auto s = psa.tag->get(pset);
          if (psai.scc.at(s) != scc)
            continue; //we left scc
        }

        state_t sucst = tag->put_or_get(suclevel, tag->size());
        if (!pa.has_state(sucst))
          pa.add_state(sucst);
        if (!pa.has_accs(sucst))  // assign priority according to resulting level
          pa.set_accs(sucst,{suclevel.prio});
        pa.adj[st][i].push_back(sucst);

        visit(sucst);
      }
    });

    auto pai = get_scc_info(pa);

    vector<state_t> copies;
    for (auto s : pa.states()) {
      if (pa.tag->get(s).states() == repps)
        copies.push_back(s);
    }


    cout << "det states: " << pa.num_states()
         << ", num scc: " << pai->sccrep.size()
         << ", num rep copy: " << copies.size() << endl;

    scc_t mintermscc;
    size_t mintermsz = pa.num_states()+1;
    for (auto s : copies) {
      auto sscc = pai->scc.at(s);
      auto ssccsz = pai->sccsz.at(sscc);
      auto sucsccs = succ_sccs(pa, *pai, sscc);
      if (sucsccs.empty() && ssccsz < mintermsz) {
        mintermscc = sscc;
        mintermsz = ssccsz;
      }
      // cout << sscc << "(sz:" << pai->sccsz.at(sscc) << ",sucs: ";
      // for (auto tmp : sucsccs)
      //   cout << tmp << ",";
      // cout << ") " << endl;
    }
    cout << mintermscc << " with " << mintermsz << endl;
    total += mintermsz;
  }
  cout << "total: " << total << endl;

  auto pap = std::make_unique<PA>(PA());
  auto& pa = *pap;

  pa.tag = std::make_unique<pa_tag_store>(pa_tag_store());
  auto& tag = pa.tag;

  return move(pap);
}

}  // namespace nbautils
