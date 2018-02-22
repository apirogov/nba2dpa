#include "relorder.hh"
#include <limits>
#include <list>
#include <vector>

using namespace nbautils;
using namespace std;

// use only with at most intmax elements!
RelOrder::RelOrder(int n) {
  for (int i = 0; i < n; i++) order.push_back(i);
  normalized = true;
  nextfree = order.size();
}

// just renumber from 0 upwards
void RelOrder::normalize() {
  ord_t i = 0;
  for (auto &it : order) it = i++;
  normalized = true;
  nextfree = order.size();
}

// get references to order elements with given relative order
// if element with invalid rank encountered, returns
vector<RelOrder::ordref> RelOrder::from_ranks(vector<RelOrder::ord_t> const &ranks) {
  // get all pointers
  vector<ordref> ptrs;
  auto it = order.begin();
  while (it != order.end()) ptrs.push_back(it++);

  // shuffle them according to given order
  vector<ordref> ret;
  for (auto rank : ranks) {
    if (/* rank < 0 || (rank is unsigned) */ rank >= order.size())
      return vector<ordref>();  // invalid rank for created size
    ret.push_back(ordref(ptrs[rank]));
  }
  return ret;
}

// obtain final new ordering from list of references (wraps normalizing and
// deref)
vector<RelOrder::ord_t> RelOrder::to_ranks(vector<RelOrder::ordref> const &refs) {
  if (!normalized) normalize();

  vector<ord_t> ret;
  for (auto ref : refs) ret.push_back(*ref);
  return ret;
}

RelOrder::ordref RelOrder::kill(RelOrder::ordref &ref) {
  normalized = false;
  if (nextfree >= std::numeric_limits<ord_t>::max())
    normalize();  // now we need to, otherwise the order breaks

  // just remove in list, add bigger in end
  order.erase(ref.ref);
  order.push_back(nextfree++);
  ref.ref = --order.end(); //inplace
  return ref; //and also return it
}
