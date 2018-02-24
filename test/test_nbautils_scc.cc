#include <vector>

#include <catch.hpp>

// #include "debug.hh"
#include "swa.hh"
#include "io.hh"
#include "scc.hh"

using namespace nbautils;

string const filedir = "test/";

TEST_CASE("SCC calculation and analysis") {
  auto const bas = parse_hoa_ba(filedir+"scc_test.hoa");
  auto const& aut = bas.front();
  auto const without = get_scc_info(*aut);
  auto const auti = get_scc_info(*aut, true);

  SECTION("get_scc_info (with and without analysis)") {
    // same in common part
    REQUIRE(without->scc == auti->scc);
    REQUIRE(without->sccrep == auti->sccrep);
    REQUIRE(without->sccsz == auti->sccsz);
    REQUIRE(without->unreachable == auti->unreachable);

    //check expected results
    REQUIRE(auti->unreachable == set<state_t>{11,12});

    REQUIRE(auti->scc.at(0) != auti->scc.at(1));
    REQUIRE(auti->scc.at(0) != auti->scc.at(3));
    REQUIRE(auti->scc.at(1) == auti->scc.at(2));
    REQUIRE(auti->scc.at(3) != auti->scc.at(4));
    REQUIRE(auti->scc.at(3) != auti->scc.at(7));
    REQUIRE(auti->scc.at(3) != auti->scc.at(10));
    REQUIRE(auti->scc.at(4) == auti->scc.at(5));
    REQUIRE(auti->scc.at(5) != auti->scc.at(6));
    REQUIRE(auti->scc.at(7) == auti->scc.at(8));
    REQUIRE(auti->scc.at(8) != auti->scc.at(9));
    REQUIRE(auti->num_sccs() == 8);

    REQUIRE(auti->sccsz.at(auti->scc.at(0)) == 1);
    REQUIRE(auti->sccsz.at(auti->scc.at(1)) == 2);
    REQUIRE(auti->sccsz.at(auti->scc.at(3)) == 1);
    REQUIRE(auti->sccsz.at(auti->scc.at(4)) == 2);
    REQUIRE(auti->sccsz.at(auti->scc.at(6)) == 1);
    REQUIRE(auti->sccsz.at(auti->scc.at(7)) == 2);
    REQUIRE(auti->sccsz.at(auti->scc.at(9)) == 1);
    REQUIRE(auti->sccsz.at(auti->scc.at(10)) == 1);

    REQUIRE(contains(auti->accepting, auti->scc.at(9)));
    REQUIRE(contains(auti->accepting, auti->scc.at(9)));
    REQUIRE(contains(auti->accepting, auti->scc.at(10)));
    REQUIRE(contains(auti->accepting, auti->scc.at(5)));
    REQUIRE(contains(auti->rejecting, auti->scc.at(0)));
    REQUIRE(contains(auti->rejecting, auti->scc.at(2)));
    REQUIRE(contains(auti->rejecting, auti->scc.at(3)));
    REQUIRE(contains(auti->rejecting, auti->scc.at(6)));
    REQUIRE(!contains(auti->rejecting, auti->scc.at(8)));
    REQUIRE(!contains(auti->accepting, auti->scc.at(8)));
    REQUIRE(contains(auti->trivial, auti->scc.at(3)));
    REQUIRE(contains(auti->trivial, auti->scc.at(6)));
    REQUIRE(contains(auti->trivial, auti->scc.at(9)));
  }

  SECTION("testing scc_states and succ_sccs") {
    auto sccof_equals = [&](state_t s, vector<state_t> const& sts){
      REQUIRE(scc_states(*aut, *auti, auti->scc.at(s)) == sts);
    };
    for (auto v : vector<state_t>{0,3,6,9,10})
      sccof_equals(v,{v});
    sccof_equals(1,{1,2});
    sccof_equals(4,{4,5});
    sccof_equals(7,{7,8});

    auto ssucc = [&](state_t p) { return succ_sccs(*aut, *auti, auti->scc.at(p)); };
    auto sccof_succ_sccof = [&](state_t p, state_t q){
      return contains(ssucc(p), auti->scc.at(q));
    };
    REQUIRE(ssucc(1).empty());
    REQUIRE(ssucc(6).empty());
    REQUIRE(ssucc(9).empty());
    REQUIRE(ssucc(10).empty());
    REQUIRE(sccof_succ_sccof(0,1));
    REQUIRE(sccof_succ_sccof(0,3));
    REQUIRE(sccof_succ_sccof(3,4));
    REQUIRE(sccof_succ_sccof(3,7));
    REQUIRE(sccof_succ_sccof(3,10));
    REQUIRE(sccof_succ_sccof(4,6));
    REQUIRE(sccof_succ_sccof(8,9));
  }

  SECTION("dead SCCs + trim BA") {
    auto const dead = get_dead_sccs(*aut, *auti);
    for (auto v : vector<state_t>{1,2,6,9})
      REQUIRE(contains(dead, auti->scc.at(v)));

    trim_ba(*aut, *auti, dead);
    REQUIRE(auti->unreachable.empty());
    for (auto v : vector<state_t>{1,2,6,9,11,12}) {
      REQUIRE(!aut->has_state(v));
      REQUIRE(!map_has_key(auti->scc, v));
    }
    for (auto s : dead) {
      REQUIRE(!map_has_key(auti->sccrep, s));
      REQUIRE(!map_has_key(auti->sccsz, s));
      REQUIRE(!contains(auti->accepting, s));
      REQUIRE(!contains(auti->rejecting, s));
      REQUIRE(!contains(auti->trivial, s));
    }

    //trimming again does nothing
    REQUIRE(get_dead_sccs(*aut, *auti).empty());
    REQUIRE(trim_ba(*aut, *auti, {})==0);
  }

  //TODO randomized tests on generated NBAs
}
