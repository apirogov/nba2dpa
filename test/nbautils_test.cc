#include <vector>

#include <catch.hpp>
// #include <spdlog/spdlog.h>

#include "relorder.hh"
#include "triebimap.hh"
#include "interfaces.hh"

using namespace std;
using namespace nbautils;

TEST_CASE("Testing the relative order structure", "[relorder]")
{
    // spdlog::set_level(spdlog::level::debug);

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

	//TODO some randomized tests

}

TEST_CASE("Testing the setmap structure (adding and removing)", "[setmap]")
{
	trie_bimap<char, int> sbm;

	SECTION("Check empty bimap", "[setmap:1]") {
		REQUIRE(sbm.size()==0);
		REQUIRE(!sbm.has(vector<char>{}));
		REQUIRE(!sbm.has(vector<char>{0,1,2}));
		REQUIRE(!sbm.has(1));
		REQUIRE(!sbm.has(2));
	}

	sbm.put(vector<char>{0,1,2}, 1);
	sbm.put(vector<char>{0,1,3}, 2);
	sbm.put(vector<char>{0,2}, 3);

	SECTION("put / get / has", "[setmap:2]") {
		REQUIRE(sbm.size()==3);
		REQUIRE(sbm.has(vector<char>{0,1,2}));
		REQUIRE(sbm.has(vector<char>{0,1,3}));
		REQUIRE(sbm.has(vector<char>{0,2}));
		REQUIRE(sbm.has(1));
		REQUIRE(sbm.has(2));
		REQUIRE(sbm.has(3));
		REQUIRE(sbm.get(vector<char>{0,1,2})==1);
		REQUIRE(sbm.get(vector<char>{0,1,3})==2);
		REQUIRE(sbm.get(vector<char>{0,2})==3);
		REQUIRE(sbm.get(1)==vector<char>{0,1,2});
		REQUIRE(sbm.get(2)==vector<char>{0,1,3});
		REQUIRE(sbm.get(3)==vector<char>{0,2});
	}

	sbm.put(vector<char>{0,1,2}, 4);

	SECTION("Overwriting an element", "[setmap:3]") {
		REQUIRE(sbm.size()==3);
		REQUIRE(!sbm.has(1));
		REQUIRE(sbm.has(4));
		REQUIRE(sbm.get(4)==vector<char>{0,1,2});
		REQUIRE(sbm.get(vector<char>{0,1,2})==4);
	}

	SECTION("put_or_get function", "[setmap:4]") {
		REQUIRE(sbm.put_or_get(vector<char>{0,1,2},5) == 4);
		REQUIRE(sbm.put_or_get(vector<char>{0,1,4},5) == 5);
		REQUIRE(sbm.size()==4);;
		REQUIRE(sbm.has(5));
		REQUIRE(sbm.get(5)==vector<char>{0,1,4});
		REQUIRE(sbm.get(vector<char>{0,1,4})==5);
	}
}

TEST_CASE("Generic trie-backed bimap", "[gen_trie_bimap]")
{
  auto from = [](string const& s){return vector<char>(s.begin(), s.end());};
  auto to = [](vector<char> const& cs){return string(cs.begin(), cs.end());};

  // auto sbmp(new generic_trie_bimap<string, char, int>(from,to));
  auto sbmp(new naive_bimap<string,int>());
  auto &sbm = *sbmp;

	SECTION("Check generic empty bimap", "[gen_trie_bimap:1]") {
		REQUIRE(sbm.size()==0);
		REQUIRE(!sbm.has("hello"));
		REQUIRE(!sbm.has(2));
	}

	sbm.put("abc", 1);
	sbm.put("abd", 2);
	sbm.put("ax", 3);

	SECTION("put / get / has of generic trie bimap", "[gen_trie_bimap:2]") {
		REQUIRE(sbm.size()==3);
		REQUIRE(sbm.has("abc"));
		REQUIRE(sbm.has("abd"));
		REQUIRE(sbm.has("ax"));
		REQUIRE(sbm.has(1));
		REQUIRE(sbm.has(2));
		REQUIRE(sbm.has(3));
		REQUIRE(sbm.get("abc")==1);
		REQUIRE(sbm.get("abd")==2);
		REQUIRE(sbm.get("ax")==3);
		REQUIRE(sbm.get(1)=="abc");
		REQUIRE(sbm.get(2)=="abd");
		REQUIRE(sbm.get(3)=="ax");
	}

	sbm.put("abc", 4);
	SECTION("Overwriting an element in gen bimap", "[gen_trie_bimap:3]") {
		REQUIRE(sbm.size()==3);
		REQUIRE(!sbm.has(1));
		REQUIRE(sbm.has(4));
		REQUIRE(sbm.get(4)=="abc");
		REQUIRE(sbm.get("abc")==4);
	}

	SECTION("put_or_get function in gen bimap", "[gen_trie_bimap:4]") {
		REQUIRE(sbm.put_or_get("abc",5) == 4);
		REQUIRE(sbm.put_or_get("abe",5) == 5);
		REQUIRE(sbm.size()==4);;
		REQUIRE(sbm.has(5));
		REQUIRE(sbm.get(5)=="abe");
		REQUIRE(sbm.get("abe")==5);
	}

    //TODO: implement and test erase

}
