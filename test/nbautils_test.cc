#include <vector>

#include <catch.hpp>
// #include <spdlog/spdlog.h>

#include "relorder.hh"
#include "triebimap.hh"

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
