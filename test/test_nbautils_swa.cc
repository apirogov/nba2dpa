#include <vector>
#include <queue>

#include <catch.hpp>

#include "types.hh"
// #include "debug.hh"
#include "io.hh"
#include "scc.hh"
#include "level.hh"

using namespace nbautils;

string const filedir = "test/";

TEST_CASE("Parsing an NBA", "[parse_ba]") {
  auto bas = parse_hoa_ba("non_existing.file");
  REQUIRE(bas.empty());

  bas = parse_hoa_ba(filedir+"test.hoa");
  REQUIRE(!bas.empty());
  auto &ba = *bas.back();

  REQUIRE(ba.get_name() == "Test NBA");
  REQUIRE(ba.get_aps() == vector<string>{"x"});
  REQUIRE(ba.states() == vector<state_t>{0,1,2,3,4,5,6,7,8,9,10,11,12});
  REQUIRE(ba.num_states() == 13);
  REQUIRE(ba.num_syms() == 2);
  REQUIRE(!ba.has_accs(0));
  REQUIRE(ba.get_accs(10) == vector<acc_t>{0});
  REQUIRE(ba.succ_raw(8,0)==vector<state_t>{9});
  REQUIRE(ba.succ_raw(8,1)==vector<state_t>{7});
}

TEST_CASE("SWA construction and simple methods", "[swa]") {
  BA aut("test",{"y","x"},{0,1});
  SECTION("initializer constructor") {
    REQUIRE(aut.get_name()=="test");

    REQUIRE(aut.get_aps()==vector<string>{"y","x"});
    REQUIRE(aut.num_syms()==4);

    REQUIRE(aut.num_states()==2);
    REQUIRE(aut.has_state(0));
    REQUIRE(aut.has_state(1));
    REQUIRE(!aut.has_state(2));

    REQUIRE(aut.get_init() == vector<state_t>{0,1});
    REQUIRE(aut.is_init(0));
    REQUIRE(aut.is_init(1));
    REQUIRE(!aut.is_init(2));
  }

  BA aut2;
  aut2.set_name(aut.get_name());
  aut2.set_init(aut.get_init());
  aut2.set_aps(aut.get_aps());
  SECTION("empty constructor") {
    REQUIRE(aut.get_name() == aut2.get_name());
    REQUIRE(aut.get_aps() == aut2.get_aps());
    REQUIRE(aut.get_init() == aut2.get_init());
  }

  aut.add_state(2);
  aut.add_state(3);
  aut.set_accs(0, {0,1});
  aut.set_accs(1, {2});
  SECTION("simple methods") {
    REQUIRE(aut.has_state(2));
    REQUIRE(aut.get_accsets() == set<acc_t>{0,1,2});
    REQUIRE(aut.has_accs(0));
    REQUIRE(aut.has_accs(1));
    REQUIRE(!aut.has_accs(2));
    REQUIRE(aut.get_accs(0)==vector<acc_t>{0,1});
    REQUIRE(aut.get_accs(1)==vector<acc_t>{2});
    REQUIRE(aut.get_accs(2)==vector<acc_t>{});
    REQUIRE(aut.states() == vector<state_t>{0,1,2});
  }
}

TEST_CASE("Level update", "[lvl_upd]") {
  auto bas = parse_hoa_ba(filedir+"test2.hoa");
  auto &ba = *bas.back();
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
    Level ini(lvc,vector<Level::state_t>(cbegin(ba.get_init()), cend(ba.get_init())));
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
  //TODO: constructor with arguments test
  BA ba;

  ba.set_aps({"x"});
  //TODO: assert throws when set aps again
  ba.add_state(0);
  ba.set_init({0});
  ba.add_state(1);
  ba.set_accs(1,{0});

  REQUIRE(ba.num_states() == 2);
  REQUIRE(!ba.has_accs(0));
  REQUIRE(ba.has_accs(1));
  REQUIRE(ba.get_accs(1) == vector<acc_t>{0});
  //TODO test other methods
}

