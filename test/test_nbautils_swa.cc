#include <vector>
#include <queue>

#include <catch.hpp>

// #include "debug.hh"
#include "swa.hh"
#include "io.hh"
#include "level.hh"

using namespace nbautils;

string const filedir = "test/";

TEST_CASE("Parsing an NBA", "[parse_ba]") {
  auto bas = parse_hoa("non_existing.file");
  REQUIRE(bas.empty());

  bas = parse_hoa(filedir+"scc_test.hoa");
  REQUIRE(!bas.empty());
  auto &ba = *bas.back();

  REQUIRE(ba.get_name() == "SCC Test NBA");
  REQUIRE(ba.get_aps() == vector<string>{"b"});
  REQUIRE(ba.states() == vector<state_t>{0,1,2,3,4,5,6,7,8,9,10,11,12});
  REQUIRE(ba.num_states() == 13);
  REQUIRE(ba.num_syms() == 2);
  REQUIRE(ba.get_init() == vector<state_t>{0});
  REQUIRE(ba.tag->geti(2)=="2");
  REQUIRE(!ba.has_accs(0));
  REQUIRE(ba.get_accs(10) == vector<acc_t>{0});
  REQUIRE(ba.succ_raw(8,0)==vector<state_t>{9});
  REQUIRE(ba.succ_raw(8,1)==vector<state_t>{7});
}

TEST_CASE("SWA construction and simple methods", "[swa]") {
  SWA<string> aut(Acceptance::BUCHI, "test",{"y","x"},{0,1});
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

  SECTION("empty constructor") {
    SWA<string> aut2;
    aut2.set_name(aut.get_name());
    for (auto s : aut.get_init())
      aut2.add_state(s);
    aut2.set_init(aut.get_init());
    aut2.set_aps(aut.get_aps());
    REQUIRE(aut.get_name() == aut2.get_name());
    REQUIRE(aut.get_aps() == aut2.get_aps());
    REQUIRE(aut.get_init() == aut2.get_init());
  }

  SECTION("acceptance sets") {
    SWA<string> aut3;
    aut3.acond = Acceptance::UNKNOWN;
    aut3.add_state(0);
    aut3.add_state(1);
    aut3.add_state(2);
    aut3.set_accs(0, {0});
    aut3.set_accs(1, {2,3});
    REQUIRE(aut3.get_accsets() == vector<acc_t>{0,2,3});
    REQUIRE(aut3.has_accs(0));
    REQUIRE(aut3.has_accs(1));
    REQUIRE(!aut3.has_accs(2));
    REQUIRE(aut3.get_accs(0)==vector<acc_t>{0});
    REQUIRE(aut3.get_accs(1)==vector<acc_t>{2,3});
    REQUIRE(aut3.get_accs(2)==vector<acc_t>{});

    SECTION("modifying acceptance sets (in)directly") {
      aut3.remove_states({0});
      REQUIRE(!aut3.has_state(0));
      aut3.set_accs(1, {});
      REQUIRE(!aut3.has_accs(1));
      aut3.set_accs(2, {4});
      REQUIRE(aut3.get_accs(2)==vector<acc_t>{4});
      REQUIRE(aut3.get_accsets() == vector<acc_t>{4});
    }
  }

  aut.add_state(2);
  aut.add_state(3);
  aut.set_succs(1, 0, {1,2});
  aut.set_succs(1, 1, {3});
  SECTION("simple methods") {
    REQUIRE(aut.has_state(2));
    REQUIRE(aut.states() == vector<state_t>{0,1,2,3});
    REQUIRE(aut.outsyms(0) == vector<sym_t>{});
    REQUIRE(aut.outsyms(1) == vector<sym_t>{0,1});
    REQUIRE(aut.state_has_outsym(1,0));
    REQUIRE(!aut.state_has_outsym(1,2));
    REQUIRE(aut.succ(0) == vector<state_t>{});
    REQUIRE(aut.succ(1) == vector<state_t>{1,2,3});
    REQUIRE(aut.succ(1,0) == vector<state_t>{1,2});
    REQUIRE(aut.succ(1,1) == vector<state_t>{3});
    REQUIRE(aut.succ(1,2) == vector<state_t>{});
  }

  aut.add_state(4);
  aut.set_succs(1,0,{1,2,4});
  aut.set_succs(1,2,{3,4});
  aut.set_succs(4,0,{0,1,2});
  aut.set_succs(4,2,{0,1,4});
  aut.set_accs(4, {0});
  aut.tag->put("1337",4);

  SECTION("removing edges and states") {
    aut.set_succs(4,1,{1,2,3});
    REQUIRE(aut.outsyms(4) == vector<sym_t>{0,1,2});
    aut.set_succs(4,1,{});
    REQUIRE(aut.outsyms(4) == vector<sym_t>{0,2});

    REQUIRE(aut.tag->geti(4)=="1337");
    REQUIRE(aut.num_states()==5);
    aut.remove_states({4});

    REQUIRE(aut.num_states()==4);
    REQUIRE(!aut.has_state(4));
    REQUIRE(!aut.tag->hasi(4));
    REQUIRE(aut.succ(1)==vector<state_t>{1,2,3});
    aut.add_state(4);
    REQUIRE(!aut.has_accs(4));
  }

  SECTION("inserting another automaton as subgraph") {
    SWA<string> aut4(Acceptance::BUCHI, "topaste", aut.get_aps(), {5});
    aut4.add_state(6);
    aut4.set_succs(5,3,{6});
    aut4.set_accs(6, {0});
    aut4.tag->put("123", 6);

    aut.insert(aut4);
    REQUIRE(aut.num_states() == 7);
    REQUIRE(aut.has_state(5));
    REQUIRE(aut.has_state(6));
    REQUIRE(!aut.tag->hasi(5));
    REQUIRE(aut.tag->hasi(6));
    REQUIRE(aut.tag->geti(6) == "123");
    REQUIRE(!aut.has_accs(5));
    REQUIRE(aut.has_accs(6));
    REQUIRE(aut.get_accs(6) == aut4.get_accs(6));
    REQUIRE(aut.succ(5,3) == vector<state_t>{6});
  }

  SWA<string> aut5(Acceptance::BUCHI, "quotienttest", aut.get_aps(), {0});
  SECTION("merging states ('quotienting')") {
    aut5.add_state(1);
    aut5.add_state(2);
    aut5.add_state(3);
    aut5.add_state(4);
    aut5.set_succs(0,0,{1});
    aut5.set_succs(0,1,{2});
    aut5.set_succs(1,0,{2,3});
    aut5.set_succs(1,1,{3});
    aut5.set_succs(2,1,{2});
    aut5.set_succs(2,2,{3});
    aut5.set_succs(2,3,{4});

    //trivial merge
    aut5.merge_states({},1);
    REQUIRE(aut5.num_states()==5);

    //merge with non-trivial edge remappings A -> (B<->C) -> D
    aut5.merge_states({2},1);
    REQUIRE(aut5.num_states()==4);
    REQUIRE(!aut5.has_state(2));
    REQUIRE(contains(aut5.succ(0,0), 1));
    REQUIRE(contains(aut5.succ(0,1), 1));
    REQUIRE(contains(aut5.succ(1,0), 1));
    REQUIRE(contains(aut5.succ(1,1), 1));
    REQUIRE(contains(aut5.succ(1,0), 3));
    REQUIRE(contains(aut5.succ(1,1), 3));
    REQUIRE(contains(aut5.succ(1,2), 3));
    REQUIRE(contains(aut5.succ(1,3), 4));


    SECTION("normalizing state numbering") {
      aut5.set_init({0,3});
      aut5.set_accs(3,{0});
      aut5.set_succs(3,1,{4});
      aut5.tag->put("123", 3);
      auto nm = aut5.normalize();
      REQUIRE(nm == map<state_t,state_t>{{0,0},{1,1},{3,2},{4,3}});
      REQUIRE(!aut5.has_state(4));
      REQUIRE(aut5.has_state(2));
      REQUIRE(aut5.get_init() == vector<state_t>{0,2});
      REQUIRE(aut5.has_accs(2));
      REQUIRE(!aut5.has_accs(3));
      REQUIRE(aut5.succ(2,1) == vector<state_t>{3});
      REQUIRE(aut5.succ(3,1) == vector<state_t>{});
      REQUIRE(!aut5.tag->hasi(3));
      REQUIRE(aut5.tag->hasi(2));
      REQUIRE(aut5.tag->geti(2) == "123");
    }
  }

  SECTION("helpers") {
    REQUIRE(powersucc(aut, vector<state_t>{}, 0) == vector<small_state_t>{});
    REQUIRE(powersucc(aut, vector<state_t>{1}, 3) == vector<small_state_t>{});
    // REQUIRE(powersucc(aut, vector<state_t>{6}, 0) == vector<small_state_t>{});
    REQUIRE(powersucc(aut, vector<state_t>{1}, 0) == to_small_state_t(aut.succ(1,0)));
    REQUIRE(powersucc(aut, vector<state_t>{1,4}, 0) == vector<small_state_t>{0,1,2,4});

    REQUIRE(!is_deterministic(aut));
    SWA<string> autdet(Acceptance::BUCHI, "detaut", aut.get_aps());
    REQUIRE(is_deterministic(autdet));
    autdet.add_state(0);
    autdet.add_state(1);
    autdet.add_state(2);
    autdet.set_succs(0,0,{1});
    autdet.set_succs(0,0,{2});
    autdet.set_succs(0,1,{2});
    REQUIRE(is_deterministic(autdet));
    REQUIRE(!is_colored(autdet));

    SWA<state_t> col(Acceptance::UNKNOWN, "paraut", aut.get_aps());
    col.add_state(0);
    col.add_state(1);
    col.add_state(2);
    col.set_accs(0,{1});
    col.set_accs(1,{2});
    REQUIRE(!is_colored(col));
    col.set_accs(2,{1});
    REQUIRE(is_colored(col));
    col.set_accs(2,{1,2});
    REQUIRE(!is_colored(col));
  }
}

TEST_CASE("Finding accepting sinks") {
  SWA<string> aut(Acceptance::BUCHI, "test", {"x"}, {0});
  aut.add_state(1);
  aut.add_state(2);
  aut.add_state(3);
  //0 has all selfloops but is not accepting
  aut.set_succs(0, 0, {0,1,2,3});
  aut.set_succs(0, 1, {0,1,2});
  //1 is accepting but has not all selfloops
  aut.set_accs(1, {0});
  aut.set_succs(1, 0, {3});
  aut.set_succs(1, 1, {0,1});

  aut.set_succs(2, 0, {0,2});
  aut.set_succs(2, 1, {1,2});
  aut.set_succs(3, 0, {2,3});
  aut.set_succs(3, 1, {1,3});

  auto aut_st = aut.states();
  function<bool(state_t)> aut_acc = [&aut](state_t v){ return aut.has_accs(v); };
  succ_sym_fun<state_t, sym_t> const aut_xsucs = [&aut](state_t v,sym_t s){ return aut.succ(v,s); };
  outsym_fun<state_t,sym_t> const aut_osyms = [&aut](state_t v){ return aut.outsyms(v); };
  REQUIRE(get_accepting_sinks(aut_st, aut.num_syms(), aut_acc, aut_osyms, aut_xsucs) == vector<unsigned>{});

  // now 2 and 3 are accepting sinks
  aut.set_accs(2, {0});
  aut.set_accs(3, {0});
  REQUIRE(get_accepting_sinks(aut_st, aut.num_syms(), aut_acc, aut_osyms, aut_xsucs) == vector<unsigned>{2,3});
}

/*
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
*/

