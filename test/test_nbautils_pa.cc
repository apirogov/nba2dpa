#include <vector>

#include <iostream>

#include <catch.hpp>

// #include "debug.hh"
#include "swa.hh"
#include "io.hh"
#include "pa.hh"
#include "common/acceptance.hh"
#include "common/util.hh"

using namespace nbautils;
using namespace std;

string const filedir = "test/";

void test_from_to_from_id(PAType from, PAType to, vector<acc_t> const& v) {
  auto f = priority_transformer(from, to, v.front(), v.back());
  // cout << seq_to_str(vec_fmap(v, f)) << endl;
  auto fmin = min(f(v.front()), f(v.back()));
  auto fmax = max(f(v.front()), f(v.back()));
  auto g = priority_transformer(to, from, fmin, fmax);
  auto h = [&](acc_t p){ return g(f(p)); };
  for (auto p : v)
    REQUIRE(p == h(p));
}

TEST_CASE("Parity condition transforms") {
  vector<vector<acc_t>> testvecs;
  testvecs.push_back({0,2,3,4,5,7});
  testvecs.push_back({1,2,4,5,6,8});
  testvecs.push_back({0,2,3,5,6,8});
  testvecs.push_back({1,2,4,5,6,7});

  SECTION("priority transformer from to . to from = id") {
    vector<PAType> pts{PAType::MIN_EVEN, PAType::MIN_ODD, PAType::MAX_EVEN, PAType::MAX_ODD};
    for (auto const from : pts) {
      for (auto const to : pts) {
        for (auto const v : testvecs) {
          /*
          cout << seq_to_str(v) << ": "<< (pa_acc_is_min(from)  ? "min" : "max") << " "
               << (pa_acc_is_even(from) ? "even" : "odd");
          cout << " -> "   << (pa_acc_is_min(to)  ? "min" : "max") << " "
                           << (pa_acc_is_even(to) ? "even" : "odd") << endl;
                           */
          test_from_to_from_id(from, to, v);
        }
      }
    }
  }

  SECTION("minimization preserves parity and relative order") {
    //TODO
  }
}

TEST_CASE("Parity Product states") {
  PAProdState ps(10, 1, 5, 11, 2, 6);
  REQUIRE(ps.a == 10);
  REQUIRE(ps.b == 11);
  REQUIRE(ps.prio == 0);
  REQUIRE(ps.priord == vector<pair<bool,int>>{{false,2},{false,4},{false,6},{true,2},{true,4},{true,6}});

  // cout << ps.to_string();
  ps = ps.succ(12, 3, 13, 2);
  // cout << ps.to_string();
  REQUIRE(ps.a == 12);
  REQUIRE(ps.b == 13);
  REQUIRE(ps.prio == 3);
  REQUIRE(ps.priord == vector<pair<bool,int>>{{false,2},{true,2},{true,4},{true,6},{false,4},{false,6}});
  ps = ps.succ(0, 4, 0, 2);
  REQUIRE(ps.prio == 4);
  REQUIRE(ps.priord == vector<pair<bool,int>>{{false,2},{true,2},{true,4},{true,6},{false,4},{false,6}});
  ps = ps.succ(0, 3, 0, 3);
  REQUIRE(ps.prio == 5);
  REQUIRE(ps.priord == vector<pair<bool,int>>{{false,2},{true,2},{false,4},{false,6},{true,4},{true,6}});
}
