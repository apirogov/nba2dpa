#include <vector>

#include <catch.hpp>

// #include "debug.hh"
#include "swa.hh"
#include "io.hh"
#include "common/scc.hh"

using namespace nbautils;

string const filedir = "test/";

TEST_CASE("SCC calculation and analysis") {
  auto const bas = parse_hoa(filedir+"scc_test.hoa");
  auto const& aut = bas.front();
  succ_fun<state_t> aut_sucs = [&aut](state_t v){ return aut->succ(v); };
  function<bool(state_t)> aut_acc = [&aut](state_t v){ return aut->has_accs(v); };
  auto const auti = get_sccs(aut->states(), aut_sucs);
  auto const autcl = ba_classify_sccs(*auti, aut_acc);

  SECTION("get_scc_info (and classify)") {
    //check expected results
    auto unreachable = unreachable_states(aut->states(), aut->get_init().front(), aut_sucs);
    REQUIRE(unreachable == vector<state_t>{11,12});

    REQUIRE(auti->scc_of.at(0) != auti->scc_of.at(1));
    REQUIRE(auti->scc_of.at(0) != auti->scc_of.at(3));
    REQUIRE(auti->scc_of.at(1) == auti->scc_of.at(2));
    REQUIRE(auti->scc_of.at(3) != auti->scc_of.at(4));
    REQUIRE(auti->scc_of.at(3) != auti->scc_of.at(7));
    REQUIRE(auti->scc_of.at(3) != auti->scc_of.at(10));
    REQUIRE(auti->scc_of.at(4) == auti->scc_of.at(5));
    REQUIRE(auti->scc_of.at(5) != auti->scc_of.at(6));
    REQUIRE(auti->scc_of.at(7) == auti->scc_of.at(8));
    REQUIRE(auti->scc_of.at(8) != auti->scc_of.at(9));
    REQUIRE(auti->sccs.size() == 10);

    REQUIRE(auti->sccs.at(auti->scc_of.at(0)).size() == 1);
    REQUIRE(auti->sccs.at(auti->scc_of.at(1)).size() == 2);
    REQUIRE(auti->sccs.at(auti->scc_of.at(3)).size() == 1);
    REQUIRE(auti->sccs.at(auti->scc_of.at(4)).size() == 2);
    REQUIRE(auti->sccs.at(auti->scc_of.at(6)).size() == 1);
    REQUIRE(auti->sccs.at(auti->scc_of.at(7)).size() == 2);
    REQUIRE(auti->sccs.at(auti->scc_of.at(9)).size() == 1);
    REQUIRE(auti->sccs.at(auti->scc_of.at(10)).size() == 1);

    REQUIRE(contains(autcl->accepting, auti->scc_of.at(9)));
    REQUIRE(contains(autcl->accepting, auti->scc_of.at(9)));
    REQUIRE(contains(autcl->accepting, auti->scc_of.at(10)));
    REQUIRE(contains(autcl->accepting, auti->scc_of.at(5)));
    REQUIRE(contains(autcl->rejecting, auti->scc_of.at(0)));
    REQUIRE(contains(autcl->rejecting, auti->scc_of.at(2)));
    REQUIRE(contains(autcl->rejecting, auti->scc_of.at(3)));
    REQUIRE(contains(autcl->rejecting, auti->scc_of.at(6)));
    REQUIRE(!contains(autcl->rejecting, auti->scc_of.at(8)));
    REQUIRE(!contains(autcl->accepting, auti->scc_of.at(8)));
    auto trivial = trivial_sccs(*auti, aut_sucs);
    REQUIRE(contains(trivial, auti->scc_of.at(3)));
    REQUIRE(contains(trivial, auti->scc_of.at(6)));
    REQUIRE(contains(trivial, auti->scc_of.at(9)));
  }

  SECTION("testing scc_states and succ_sccs") {
    auto sccof_equals = [&](state_t s, vector<state_t> const& sts){
      REQUIRE(auti->sccs.at(auti->scc_of.at(s)) == sts);
    };
    for (auto v : vector<state_t>{0,3,6,9,10})
      sccof_equals(v,{v});
    sccof_equals(1,{1,2});
    sccof_equals(4,{4,5});
    sccof_equals(7,{7,8});

    auto ssucc = [&](state_t p) { return succ_sccs(*auti, auti->scc_of.at(p), aut_sucs); };
    auto sccof_succ_sccof = [&](state_t p, state_t q){
      return contains(ssucc(p), auti->scc_of.at(q));
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

  /*
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
  */

  //TODO randomized tests on generated NBAs
}
