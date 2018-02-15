#include "det.hh"
#include "types.hh"

namespace nbautils {

// BFS-based determinization with supplied level update config
PA::uptr determinize(LevelConfig const& lc) {
  spd::get("log")->debug("determinize(aut)");
  auto starttime = get_time();

  auto& aut = *lc.aut;

  auto pap = std::make_unique<PA>(PA());
  auto& pa = *pap;

  pa.tag = std::make_unique<pa_tag_store>(pa_tag_store());
  auto& tag = pa.tag;

  // same letters etc
  pa.meta = aut.meta;
  // initial state is 0
  pa.init = 0;
  pa.acc[pa.init] = 0;  // priority does not matter

  // q_0 tag
  auto pinit = Level(lc, std::vector<small_state_t>{(small_state_t)aut.init});
  tag->put(pinit, pa.init);

  // explore powerset by bfs
  std::queue<state_t> bfsq;
  bfsq.push(pa.init);
  while (!bfsq.empty()) {
    auto const st = bfsq.front();
    bfsq.pop();
    if (pa.has_state(st)) continue;  // have visited this one

    // get inner states of current ps state
    auto curlevel = tag->get(st);

    for (auto i = 0; i < (int)pa.num_syms(); i++) {
      auto suclevel = curlevel.succ(lc, i);
      if (suclevel.powerset == 0) //emptyset
        continue;
      state_t sucst = tag->put_or_get(suclevel, tag->size());
      if (!pa.has_acc(sucst))  // assign priority according to resulting level
        pa.acc[sucst] = suclevel.prio;
      pa.adj[st][i].push_back(sucst);
      bfsq.push(sucst);
    }
  }

  spd::get("log")->debug("determinize completed ({:.4f} s)", get_secs_since(starttime));
  return move(pap);
}

}  // namespace nbautils
