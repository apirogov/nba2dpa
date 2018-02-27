#include <vector>

#include <iostream>

#include <catch.hpp>

// #include "debug.hh"
#include "swa.hh"
#include "io.hh"
#include "common/parityacc.hh"
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
