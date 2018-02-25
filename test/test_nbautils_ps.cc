#include <vector>

#include <catch.hpp>

// #include "debug.hh"
#include "swa.hh"
#include "io.hh"
#include "ps.hh"
#include "common/util.hh"

using namespace nbautils;

string const filedir = "test/";

TEST_CASE("Powerset construction") {
  auto const bas = parse_hoa_ba(filedir+"ps_test.hoa");
  auto const& aut = bas.front();
  // auto const auti = get_scc_info(*aut, true);
  auto const accsinks = get_accepting_sinks(*aut);

  SECTION("powerset without accepting sink collaps") {
    auto psa = powerset_construction(*aut, {});
    REQUIRE(psa->get_name() == aut->get_name());
    REQUIRE(psa->get_aps() == aut->get_aps());
    REQUIRE(psa->get_init().size() == 1);
    REQUIRE(is_deterministic(*psa));

    auto st = [&](vector<small_state_t> s){ return psa->tag->get(s); };
    auto v = [](vector<small_state_t> v){ return v; };
    auto sucst = [&](vector<small_state_t> s, sym_t x){
      return psa->tag->geti(psa->succ(st(s), x).front());
    };

    REQUIRE(psa->tag->geti(psa->get_init().front()) == to_small_state_t(aut->get_init()));
    REQUIRE(psa->num_states() == 6);

    REQUIRE(sucst(v({0}), 0) == v({0,1}));
    REQUIRE(sucst(v({0}), 1) == v({2}));
    REQUIRE(sucst(v({2}), 0) == v({0,2}));
    REQUIRE(sucst(v({2}), 1) == v({1,2}));
    REQUIRE(sucst(v({0,1}), 0) == v({0,1,2}));
    REQUIRE(sucst(v({0,1}), 1) == v({2}));
    REQUIRE(sucst(v({1,2}), 0) == v({0,2}));
    REQUIRE(sucst(v({1,2}), 1) == v({1,2}));
    REQUIRE(sucst(v({0,2}), 0) == v({0,1,2}));
    REQUIRE(sucst(v({0,2}), 1) == v({1,2}));
    REQUIRE(sucst(v({0,1,2}), 0) == v({0,1,2}));
    REQUIRE(sucst(v({0,1,2}), 1) == v({1,2}));
  }

  SECTION("powerset with accepting sink collaps") {
    REQUIRE(accsinks == vector<small_state_t>{2});
    auto psa = powerset_construction(*aut,accsinks);
    REQUIRE(psa->get_name() == aut->get_name());
    REQUIRE(psa->get_aps() == aut->get_aps());
    REQUIRE(psa->get_init().size() == 1);
    REQUIRE(is_deterministic(*psa));

    auto st = [&](vector<small_state_t> s){ return psa->tag->get(s); };
    auto v = [](vector<small_state_t> v){ return v; };
    auto sucst = [&](vector<small_state_t> s, sym_t x){
      return psa->tag->geti(psa->succ(st(s), x).front());
    };

    REQUIRE(psa->tag->geti(psa->get_init().front()) == to_small_state_t(aut->get_init()));
    REQUIRE(psa->num_states() == 3);

    REQUIRE(sucst(v({0}), 0) == v({0,1}));
    REQUIRE(sucst(v({0}), 1) == v({2}));
    REQUIRE(sucst(v({0,1}), 0) == v({2}));
    REQUIRE(sucst(v({0,1}), 1) == v({2}));
    REQUIRE(sucst(v({2}), 0) == v({2}));
    REQUIRE(sucst(v({2}), 1) == v({2}));
  }

  SECTION("edge cases") {
    BA aut("edgecase",{"x"});
    //empty -> emptyset
    auto psa = powerset_construction(aut, {});
    REQUIRE(psa->num_states() == 1);
    REQUIRE(psa->tag->geti(0) == vector<small_state_t>{});
    //singleton, no initial -> emptyset
    aut.add_state(0);
    psa = powerset_construction(aut, {});
    REQUIRE(psa->num_states() == 1);
    REQUIRE(psa->tag->geti(0) == vector<small_state_t>{});
    //singleton state
    aut.set_init({0});
    psa = powerset_construction(aut, {});
    REQUIRE(psa->tag->geti(0) == vector<small_state_t>{0});
    REQUIRE(psa->tag->geti(1) == vector<small_state_t>{});
    REQUIRE(psa->num_states() == 2);
  }

  SECTION("Powerset of deterministic automaton is same automaton") {
    auto psa = powerset_construction(*aut, {});
    auto psa2 = powerset_construction(*psa, {});
    REQUIRE(psa->num_states() == psa2->num_states());
    REQUIRE(is_deterministic(*psa));
    REQUIRE(is_deterministic(*psa2));

    // simultaneous traversal of det. automaton to witness isomorphy
    map<state_t, state_t> m;
    state_t start = psa->get_init().front();
    m[start] = psa2->get_init().front();
    set<state_t> discovered;
    bfs(start, [&](auto const& st, auto const& visit, auto const&) {
        auto st2 = m.at(st);
        REQUIRE(psa->outsyms(st) == psa2->outsyms(st2));
        for (auto sym : psa->outsyms(st)) {
          state_t suc = psa->succ(st, sym).front();
          state_t suc2 = psa2->succ(st2, sym).front();

          if (!contains(discovered, suc)) {
            REQUIRE(!map_has_key(m, suc));
            discovered.emplace(suc);
            m[suc] = suc2;
          } else {
            REQUIRE(map_has_key(m, suc));
            REQUIRE(m.at(suc)==suc2);
          }

          visit(suc);
        }
    });
  }

  //TODO powerset product test

  //TODO randomized tests on generated NBAs
}
