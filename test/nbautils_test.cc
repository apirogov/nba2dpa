#include <vector>
#include <queue>

#include <catch.hpp>
// #include <spdlog/spdlog.h>

#include "types.hh"
#include "debug.hh"
#include "io.hh"
#include "scc.hh"
#include "relorder.hh"
#include "level.hh"
#include "triebimap.hh"
#include "interfaces.hh"
#include "util.hh"

using namespace std;
using namespace nbautils;

#include <spdlog/spdlog.h>
namespace spd = spdlog;

string const filedir = "test/";

//TODO: find better way to include logger in code?
auto logger = spd::stdout_logger_mt("log");

TEST_CASE("Generic bfs", "[parse_ba]") {
  auto pba = parse_ba(filedir+"test.hoa");
  auto &ba = *pba;
  auto scci = get_scc_info(ba,true);
  auto &bai = *scci;
  trim_ba(ba, bai);

  int count = 0;
  bfs(ba.init, [&](state_t const& st, auto v){
      count ++;
      for (auto s : ba.succ(st))
        v(s);
  });
  REQUIRE(count == ba.num_states());
}

TEST_CASE("Parsing an NBA", "[parse_ba]") {
  auto pba = parse_ba("non_existing.file");
  REQUIRE(!pba);

  pba = parse_ba(filedir+"test.hoa");
  auto &ba = *pba;
  REQUIRE(ba.meta.name == "Test NBA");
  REQUIRE(ba.aps() == vector<string>{"x"});
  REQUIRE(ba.states() == vector<state_t>{0,1,2,3,4,5,6,7,8,9,10,11,12});
  REQUIRE(ba.num_states() == 13);
  REQUIRE(ba.num_syms() == 2);
  REQUIRE(!ba.has_accs(0));
  REQUIRE(ba.get_accs(10) == vector<acc_t>{0});
  REQUIRE(ba.succ_raw(8,0)==vector<state_t>{9});
  REQUIRE(ba.succ_raw(8,1)==vector<state_t>{7});
}

TEST_CASE("Level update", "[lvl_upd]") {
  auto pba = parse_ba(filedir+"test2.hoa");
  auto &ba = *pba;
  auto scci = get_scc_info(ba,true);
  auto &bai = *scci;
  // printBA(ba,bai);
  // printSCCI(bai);
  // trim_ba(ba,bai);

  LevelConfig lvc;
  lvc.aut = &ba;
  lvc.auti = &bai;
  lvc.sep_rej = false;
  lvc.sep_acc = false;

  queue<sym_t> w;
  w.push(0); w.push(1);
  Level cur;

  bool debug = false;
  auto step = [&](){
    auto s = w.front();
    w.pop();
    w.push(s);
    cur = cur.succ(lvc, s);
    cout << cur.to_string() << endl;
  };
  auto runw = [&](){
    Level ini(lvc,{(Level::state_t)ba.init});
    cout << ini.to_string() << endl;
    cur = ini;
    for (int i=0; i<8; i++)
      step();
    cout << "----" << endl;
  };
  auto run_multi = [&]() {
    runw();
    lvc.sep_rej = true;
    runw();
    lvc.sep_acc = true;
    runw();
    lvc.sep_rej = false;
    runw();
    lvc.sep_acc = false;
  };
  cout << "---- MUELLERSCHUPP ----" << endl;
  lvc.update = LevelUpdateMode::MUELLERSCHUPP;
  run_multi();
  cout << "---- SAFRA ----" << endl;
  lvc.update = LevelUpdateMode::SAFRA;
  // debug = true;
  run_multi();
}

TEST_CASE("Testing the SWA interface", "[swa]") {
  //TODO: constructor enforces creation of initial
  //TODO: aps are not allowed to be modified
  BA ba;
  ba.meta.aps = {"x"};
  ba.add_state(0);
  ba.init = 0;
  ba.add_state(1);
  ba.set_accs(1,{0});

  REQUIRE(ba.num_states() == 2);
  REQUIRE(!ba.has_accs(0));
  REQUIRE(ba.has_accs(1));
  REQUIRE(ba.get_accs(1) == vector<acc_t>{0});
  //TODO test other methods
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

template<typename Impl>
void test_bimap_interface(bimap<string,int,Impl>& sbm) {
	SECTION("sanity-check empty bimap", "[bimap:1]") {
		REQUIRE(sbm.size()==0);
		REQUIRE(!sbm.has("hello"));
		REQUIRE(!sbm.has(2));
	}

	sbm.put("abc", 1);
	sbm.put("abd", 2);
	sbm.put("ax", 3);

	SECTION("put / get / has of bimap", "[bimap:2]") {
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
	SECTION("Overwriting an element in bimap", "[bimap:3]") {
		REQUIRE(sbm.size()==3);
		REQUIRE(!sbm.has(1));
		REQUIRE(sbm.has(4));
		REQUIRE(sbm.get(4)=="abc");
		REQUIRE(sbm.get("abc")==4);
	}

	SECTION("put_or_get in bimap", "[bimap:4]") {
		REQUIRE(sbm.put_or_get("abc",5) == 4);
		REQUIRE(sbm.put_or_get("abe",5) == 5);
		REQUIRE(sbm.size()==4);;
		REQUIRE(sbm.has(5));
		REQUIRE(sbm.get(5)=="abe");
		REQUIRE(sbm.get("abe")==5);
	}

    //TODO: test erase

}

TEST_CASE("Test naive bimap interface", "[bimap-interface-naive]")
{
  auto sbmp(new naive_bimap<string,int>());
  test_bimap_interface(*sbmp);
}

TEST_CASE("Test trie bimap interface", "[bimap-interface-trie]")
{
  auto from = [](string const& s){return vector<char>(s.begin(), s.end());};
  auto to = [](vector<char> const& cs){return string(cs.begin(), cs.end());};
  auto sbmp2(new generic_trie_bimap<string, char, int>(from,to));
  test_bimap_interface(*sbmp2);
}
//TODO: add boost bimap as possibility
