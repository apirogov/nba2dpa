#include <vector>
#include <map>

#include <catch.hpp>
// #include <rapidcheck.h>

#include "swa.hh"
#include "io.hh"
#include "common/relorder.hh"
#include "common/part_refinement.hh"
#include "common/util.hh"
#include "common/algo.hh"
#include "common/triebimap.hh"
#include "common/bimap.hh"

using namespace nbautils;

// test generic BFS by using it to count sum of state ids in a virtual full graph
TEST_CASE("Generic bfs", "[parse_ba]") {
  int count = 0;
  bfs(1, [&](auto const& st, auto const& v, auto const&){
      count += st;
      for (auto i=1; i<=10; i++)
        v(i);
  });
  REQUIRE(count == 55); //this suceeds if every "node" was visited once
}

TEST_CASE("STL helpers", "[stl_helpers]") {
  map<int, int> m = {{5,3},{3,1},{7,2},{1,2}};
  REQUIRE(map_get_keys(m) == vector<int>{1,3,5,7});
  REQUIRE(map_get_vals(m) == vector<int>{2,1,3,2});
  REQUIRE(map_has_key(m, 5));
  REQUIRE(!map_has_key(m, 6));
  vector<int> v = {5,3,7,1,1};
  set<int> s(cbegin(v), cend(v));
  REQUIRE(contains(v, 5));
  REQUIRE(!contains(v, 6));
  for (auto it : s)
    REQUIRE(contains(v, it) == contains(s, it));

  // auto w = fmap(v, [](auto const & i) { return i+1; });
  // REQUIRE(v.size() == w.size());
  // for (auto i=0; i<(int)v.size(); ++i)
  //   REQUIRE(w[i] == v[i]+1);
  // badvec=vector<int>{2,1}
  REQUIRE(is_set_vec(vector<int>{1,2,3}));
  REQUIRE(!is_set_vec(vector<int>{2,1,3}));
  REQUIRE(!is_set_vec(vector<int>{1,1,2,3}));
  REQUIRE(!is_set_vec(vector<int>{2,1,1,3}));
  v = vector<int>{2,1,1,3};
  REQUIRE(seq_to_str(v) == "2,1,1,3");
  REQUIRE(seq_to_str(set<int>{1,3,5}," ") == "1 3 5");
  vec_to_set(v);
  REQUIRE((is_set_vec(v) && v == vector<int>{1,2,3}));

  v = {1,2,3,4};
  vector<int> x = {3,4,5,6};
  REQUIRE(set_intersect(v,x) == vector<int>{3,4});
  REQUIRE(set_merge(v,x) == vector<int>{1,2,3,4,5,6});
  REQUIRE(set_diff(v,x) == vector<int>{1,2});
  REQUIRE(set_diff(x,v) == vector<int>{5,6});

  auto bs = vec_fmap(v, [](auto const& v){ return v % 2 == 0; });
  REQUIRE(bs == vector<bool>{false, true, false, true});

  v = vector<int>{1,3,5,7,2,4,6,9};
  auto samepar = [](int a, int b){return (a%2==0)==(b%2==0);};
  REQUIRE(group_by(v, samepar) == vector<vector<int>>{{1,3,5,7},{2,4,6},{9}});
}

TEST_CASE("Testing the relative order structure", "[relorder]") {
	RelOrder order = RelOrder(5);
	vector<RelOrder::ord_t> testord{4,0,2,1,3};
	auto virtord = order.from_ranks(testord);

	SECTION("Converting to and from relorder (identity) ", "[relorder:1]") {
		auto backord = order.to_ranks(virtord);
		REQUIRE(backord == testord); //must return back unchanged

		testord[4]=5; // make an input rank invalid
		virtord = order.from_ranks(testord);
		REQUIRE(virtord.size()==0); //must return empty vector as error
	}

	SECTION("Sanity check killing in relorder", "[relorder:2]") {
		auto newelem = order.kill(virtord[2]);
		REQUIRE(3 == *virtord[4]); //must be unchanged
		REQUIRE(*newelem == 5); //must get next free
		REQUIRE(virtord[2] == newelem); //must be overwritten inplace

		vector<RelOrder::ord_t> targetord{3,0,4,1,2};
		REQUIRE(*virtord[2] == 5); //must be assigned next free

		auto backord = order.to_ranks(virtord);
		REQUIRE(backord == targetord); //must be normalized again
	}

    /* TODO: permutation generator
	SECTION("Randomized tests", "[relorder:2]") {
      rc::check("to_ranks . from_ranks == id",
        [](vector<RelOrder::ord_t> const& v) {
          RelOrder::ord_t max = 0;
          for (auto i : v)
            if (i > max)
              max = i;
          RelOrder order(max+1);
          auto tokens = order.from_ranks(v);
          // if (!tokens.empty())
          //   order.kill(tokens[0]);
          auto back = order.to_ranks(tokens);
          RC_ASSERT(v == back);
        });
    }
    */
}

TEST_CASE("testing partition refinement structure", "[part_refine]") {
  vector<int> vals{0,1,2,3,4,5,6,7,8,9};
  PartitionRefiner<int> pr({vals});

  REQUIRE(pr.num_sets()==1);
  auto ini = pr.get_set_ids().front();
  REQUIRE(pr.get_set_size(ini)==10);
  for (auto v : vals)
    REQUIRE(pr.get_set_of(v)==ini);

  auto newset = pr.separate(ini, [](int i){ return i%2==0; });
  REQUIRE(newset != nullptr);
  REQUIRE(pr.get_set_ids().front()==*newset);
  REQUIRE(pr.num_sets()==2);
  REQUIRE(pr.num_sets()==pr.get_set_ids().size());
  REQUIRE(pr.get_elements_of(pr.get_set_ids().front())==vector<int>{0,2,4,6,8});
  REQUIRE(pr.get_elements_of(pr.get_set_ids().back())==vector<int>{1,3,5,7,9});

  auto newset2 = pr.separate(pr.get_set_ids().back(), [](int i){ return i%3==0; });
  REQUIRE(newset2 != nullptr);
  REQUIRE(pr.num_sets()==3);
  auto sids = pr.get_set_ids();
  REQUIRE(pr.num_sets()==sids.size());
  REQUIRE(sids[1]==*newset2);
  REQUIRE(pr.get_elements_of(sids[0])==vector<int>{0,2,4,6,8});
  REQUIRE(pr.get_elements_of(sids[1])==vector<int>{3,9});
  REQUIRE(pr.get_elements_of(sids[2])==vector<int>{1,5,7});

  for (auto sid : pr.get_set_ids()) {
    REQUIRE(pr.separate(sid, [](int i){ return i%2==0; }) == nullptr);
  }

  for (auto sid : pr.get_set_ids()) {
    auto els = pr.get_elements_of(sid);
    REQUIRE(pr.get_set_size(sid) == els.size());
    for (auto el : els) {
      REQUIRE(pr.get_set_of(el)==sid);
    }
  }

  vector<vector<int>> vals2{{1,3},{2,4},{7}};
  PartitionRefiner<int> pr2(vals2);
  REQUIRE(pr2.num_sets()==3);
  auto sids2 = pr2.get_set_ids();
  REQUIRE(pr2.get_elements_of(sids2[0])==vector<int>{1,3});
  REQUIRE(pr2.get_elements_of(sids2[1])==vector<int>{2,4});
  REQUIRE(pr2.get_elements_of(sids2[2])==vector<int>{7});
  REQUIRE(pr2.get_refined_sets() == vals2);
  for (auto sid : sids2) {
    auto els = pr2.get_elements_of(sid);
    REQUIRE(pr2.get_set_size(sid) == els.size());
    for (auto el : els) {
      REQUIRE(pr2.get_set_of(el)==sid);
    }
  }
}

TEST_CASE("minimize DFA") {
  SWA<string> aut(Acceptance::PARITY,"dfa",{"x","y"},{0});
  aut.add_state(1); aut.add_state(2); aut.add_state(3); aut.add_state(4);
  aut.set_accs(0, {0}); aut.set_accs(1, {1}); aut.set_accs(2, {2}); aut.set_accs(3, {2}); aut.set_accs(4, {1});
  aut.set_succs(0, 0, {1}); aut.set_succs(0, 1, {3}); aut.set_succs(0, 2, {4}); aut.set_succs(0, 3, {4});
  aut.set_succs(1, 0, {1}); aut.set_succs(1, 1, {2}); aut.set_succs(1, 2, {4}); aut.set_succs(1, 3, {4});
  aut.set_succs(2, 0, {0}); aut.set_succs(2, 1, {2}); aut.set_succs(2, 2, {4}); aut.set_succs(2, 3, {4});
  aut.set_succs(3, 0, {0}); aut.set_succs(3, 1, {3}); aut.set_succs(3, 2, {4}); aut.set_succs(3, 3, {4});
  aut.set_succs(4, 0, {4}); aut.set_succs(4, 1, {4}); aut.set_succs(4, 2, {4}); aut.set_succs(4, 3, {4});
  REQUIRE(is_deterministic(aut));
  REQUIRE(is_colored(aut));

  function<acc_t(state_t)> colors = [&](auto s){return aut.get_accs(s).front();};
  function<state_t(state_t,sym_t)> xsucc = [&](auto p, auto s){return aut.succ(p,s).front();};
  auto equiv = dfa_equivalent_states(aut.states(), colors, aut.num_syms(), xsucc);
  REQUIRE(equiv == vector<vector<state_t>>{{0},{1},{4},{2,3}});

  auto initial = aut.get_init().front();
  bool seenini = false;
  for (auto ecl : equiv) {
    auto rep = ecl.back();
    if (!seenini) {
      auto it = lower_bound(begin(ecl), end(ecl), initial);
      if (it != end(ecl)) {
        rep = initial;
        seenini = true;
        ecl.erase(it);
      }
    } else {
      ecl.pop_back();
    }
    aut.merge_states(ecl, rep);
  }
  REQUIRE(is_deterministic(aut));
  aut.normalize();
  REQUIRE(is_deterministic(aut));
  print_hoa(aut);
}

template<typename Impl>
void test_bimap_interface(bimap<string,int,Impl>& sbm) {
	SECTION("sanity-check empty bimap", "[bimap:1]") {
		REQUIRE(sbm.size()==0);
		REQUIRE(!sbm.has("hello"));
		REQUIRE(!sbm.hasi(2));
	}

	sbm.put("abc", 1);
	sbm.put("abd", 2);
	sbm.put("ax", 3);

	SECTION("put / get / has of bimap", "[bimap:2]") {
		REQUIRE(sbm.size()==3);
		REQUIRE(sbm.has("abc"));
		REQUIRE(sbm.has("abd"));
		REQUIRE(sbm.has("ax"));
		REQUIRE(sbm.hasi(1));
		REQUIRE(sbm.hasi(2));
		REQUIRE(sbm.hasi(3));
		REQUIRE(sbm.get("abc")==1);
		REQUIRE(sbm.get("abd")==2);
		REQUIRE(sbm.get("ax")==3);
		REQUIRE(sbm.geti(1)=="abc");
		REQUIRE(sbm.geti(2)=="abd");
		REQUIRE(sbm.geti(3)=="ax");
	}

	sbm.put("abc", 4);
	SECTION("Overwriting an element in bimap", "[bimap:3]") {
		REQUIRE(sbm.size()==3);
		REQUIRE(!sbm.hasi(1));
		REQUIRE(sbm.hasi(4));
		REQUIRE(sbm.geti(4)=="abc");
		REQUIRE(sbm.get("abc")==4);
	}

	SECTION("put_or_get in bimap", "[bimap:4]") {
		REQUIRE(sbm.put_or_get("abc",5) == 4);
		REQUIRE(sbm.put_or_get("abe",5) == 5);
		REQUIRE(sbm.size()==4);;
		REQUIRE(sbm.hasi(5));
		REQUIRE(sbm.geti(5)=="abe");
		REQUIRE(sbm.get("abe")==5);
	}

    // abc -> 4, abd -> 2, ax -> 3
    SECTION("erase in bimap", "[bimap:5]") {
      auto tmp = sbm.size();
      //remove non existing element -> no change
      REQUIRE(!sbm.hasi(1));
      sbm.erase(1);
      REQUIRE(sbm.size() == tmp);

      //remove existing element, others should be fine
      REQUIRE(sbm.hasi(3));
      sbm.erase(3);
      REQUIRE(sbm.size() == tmp-1);
      REQUIRE(!sbm.hasi(3));
      REQUIRE(sbm.geti(2)=="abd");
      REQUIRE(sbm.geti(4)=="abc");
    }
}

//run same tests with all implementations

TEST_CASE("Test naive bimap interface", "[bimap-interface-naive]")
{
  auto sbmp(new naive_ordered_bimap<string,int>());
  test_bimap_interface(*sbmp);
  auto sbmp2(new naive_unordered_bimap<string,int>());
  test_bimap_interface(*sbmp2);
}

TEST_CASE("Test trie bimap interface", "[bimap-interface-trie]")
{
  auto from = [](string const& s){return vector<char>(s.begin(), s.end());};
  auto to = [](vector<char> const& cs){return string(cs.begin(), cs.end());};
  auto sbmp2(new generic_trie_bimap<string, char, int>(from,to));
  test_bimap_interface(*sbmp2);
}

