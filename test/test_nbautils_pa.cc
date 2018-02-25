#include <vector>

#include <catch.hpp>

// #include "debug.hh"
#include "swa.hh"
#include "io.hh"
#include "pa.hh"
#include "common/util.hh"

using namespace nbautils;

string const filedir = "test/";

void test_double_flip(auto const& func, vector<acc_t> const& v) {
  auto f = func(v);
  vector<acc_t> w;
  transform(begin(v), end(v), back_inserter(w), f);
  vec_to_set(w);
  auto g = func(w);
  vector<acc_t> x;
  transform(begin(w), end(w), back_inserter(x), g);
  vec_to_set(x);
  REQUIRE(v==x);
}

void test_parity_flip(vector<acc_t> const& v) {
    auto f = flip_acc_parity(v);
    vector<acc_t> w;
    transform(begin(v), end(v), back_inserter(w), f);
    for (auto i=0; i<(int)v.size(); i++) {
      REQUIRE(!same_parity(v[i], w[i]));
    }
}

void test_polarity_flip(vector<acc_t> const& v) {
    auto f = flip_acc_polarity(v);
    vector<acc_t> w;
    transform(begin(v), end(v), back_inserter(w), f);
    for (auto i=0; i<(int)v.size(); i++) {
      REQUIRE(same_parity(v[i], w[i]));
    }
}

void test_from_to_from_id(PAType from, PAType to, vector<acc_t> const& v) {
  auto f = priority_transformer(from, to, v);
  vector<acc_t> w;
  transform(begin(v), end(v), back_inserter(w), f);
  vec_to_set(w);
  auto g = priority_transformer(to, from, w);
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

  SECTION("flip . flip = id") {
    for (auto const v : testvecs) {
      test_double_flip(flip_acc_parity,   v);
      test_double_flip(flip_acc_polarity, v);
    }
  }
  SECTION("parity flip flips parity") {
    for (auto const v : testvecs) {
      test_parity_flip(v);
    }
  }
  SECTION("polarity flip does not flip parity") {
    for (auto const v : testvecs) {
      test_polarity_flip(v);
    }
  }
  SECTION("priority transformer from to . to from = id") {
    vector<PAType> pts{PAType::MIN_EVEN, PAType::MIN_ODD, PAType::MAX_EVEN, PAType::MAX_ODD};
    for (auto const from : pts) {
      for (auto const to : pts) {
        for (auto const v : testvecs) {
          test_from_to_from_id(from, to, v);
        }
      }
    }
  }

  SECTION("minimization preserves parity and relative order") {
    //TODO
  }
}
