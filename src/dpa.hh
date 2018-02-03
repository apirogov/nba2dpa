#pragma once

#include <memory>
#include <vector>

#include "level.hh"
#include "nba.hh"
#include "relorder.hh"

namespace nbautils {

class DPA {
 public:
  std::shared_ptr<NBA> const nba_;
  typedef nbautils::NBA::state_t state_t;
  typedef nbautils::NBA::sym_t sym_t;

  // initialize a new dpa backed by given NBA
  DPA(NBA const &nba);

  // get successor level by generating it from scratch
  // TODO: create lookup table and make transparent succ using this on demand
  // and another function that fully generates the graph
  // Level calc_succ(sym_t const &x);
};

}  // namespace nbautils
