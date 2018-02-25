#include <vector>
#include <cassert>
#include <functional>
#include "swa.hh"
#include "common/util.hh"

using namespace std;
using namespace nbautils;

namespace nbautils {

//the following functions return functions that can be used to transform priorieies of
//parity automata.

//normalize sets s.t. minimal prio is 1 if odd / 0 if even
/*
function<acc_t(acc_t)> normalize_priorities(vector<acc_t> const& as) {
  assert(is_set_vec(as));
  if (as.empty()) return identity<acc_t>;
  auto const minpri = as.front();

  int const shift = minpri > 1 ? (minpri/2) * 2 : 0;
  map<acc_t, acc_t> ret;

  return [shift](acc_t v){ return v-shift; };
}
*/

//bijective priority map
function<acc_t(acc_t)> flip_acc_parity(vector<acc_t> const& as) {
  assert(is_set_vec(as));
  if (as.empty()) return identity<acc_t>;

  auto const minpri = as.front();
  if (minpri==0)
    return [](acc_t v){ return v+1; };
  else
    return [](acc_t v){ return v-1; };
}

//bijective priority map
function<acc_t(acc_t)> flip_acc_polarity(vector<acc_t> const& as) {
  assert(is_set_vec(as));
  if (as.empty()) return identity<acc_t>;

  auto const maxpri = as.back();
  auto const tmp = maxpri % 2==0 ? maxpri : maxpri+1;
  return [tmp](acc_t v){ return tmp-v; };
}

// takes source and target priority acceptance conditions and a list of priorities of the
// source type, returns function mapping source pris to semantically equiv. in target cond.
function<acc_t(acc_t)> priority_transformer(PAType from, PAType to, vector<acc_t> const& as) {
  if (from==to) //nothing to do
    return identity<acc_t>;

  auto adapt_par = same_parity(from,to)   ? identity<acc_t> : flip_acc_parity(as);
  vector<acc_t> tmp;
  transform(begin(as), end(as), back_inserter(tmp), adapt_par);
  vec_to_set(tmp);
  auto adapt_pol = same_polarity(from,to) ? identity<acc_t> : flip_acc_polarity(tmp);
  return [=](acc_t v){ return adapt_pol(adapt_par(v)); };
}

}
