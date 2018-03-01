#pragma once

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <utility>
#include <memory>

#include "common/util.hh"

namespace nbautils {
  using namespace std;

//TODO: this is similar to what happens in levels!

//partition refinement structure for a fixed number of elements
template<typename T>
class PartitionRefiner {
public:
  using vec_it = typename vector<T>::iterator;
  using bounds = pair<vec_it, vec_it>;
  using sym_set = typename list<bounds>::iterator;

private:
  vector<T> elements; //states are grouped by sets
  list<bounds> sets;  //set boundaries over elements as iterator pairs
  map<T, sym_set> set_of; //map from states back to iterator pair

 public:
  size_t num_sets() const { return sets.size(); }
  size_t get_set_size(sym_set const& ref) const { return ref->second - ref->first; }

  //start with existing sets (or singleton set)
  PartitionRefiner(vector<vector<T>> const& startsets) {
    // cout << "copy" << endl;
    for (auto s : startsets)
      copy(cbegin(s),cend(s), back_inserter(elements));

    auto l = elements.begin();
    for (auto s : startsets) {
      auto r = l+s.size();

      // cout << "sort" << endl;
      sort(l, r);
      // cout << "add pair" << endl;
      sets.emplace_back(make_pair(l,r));
      // cout << "add backmap" << endl;
      auto sym = sets.end();
      --sym;
      for (auto el : s)
        set_of[el] = sym;

      l = r;
    }
    // cout << "constr done" << endl;
  };

  //a set is identified by the iterator to the bound pair
  vector<sym_set> get_set_ids() {
    vector<sym_set> ret;
    for (auto it=sets.begin(); it!=sets.end(); ++it)
      ret.push_back(it);
    return ret;
  }
  //convenience function - returns all sets
  vector<vector<T>> get_refined_sets() {
    vector<vector<T>> ret;
    for (auto it=begin(sets); it!=end(sets); ++it) {
      ret.push_back(get_elements_of(it));
    }
    return ret;
  }

  vector<T> get_elements_of(sym_set const& ref) const {
    vector<T> ret;
    for (auto it=ref->first; it!=ref->second; ++it)
      ret.push_back(*it);
    return ret;
  }

  sym_set get_set_of(T const& el) const {
    assert(map_has_key(set_of, el));
    return set_of.at(el);
  }

  //separate given set into satisfying and not satisfying predicate
  //if both are nonempty, returns token of second set, otherwise returns nullptr
  unique_ptr<sym_set> separate(sym_set& set, function<bool(T)> const& pred) {
    auto mid = stable_partition(set->first, set->second, pred);
    //nonempty intersection and difference -> split
    if (mid != set->first && mid != set->second) {
      auto newset = make_unique<sym_set>(sets.insert(set, make_pair(set->first, mid)));
      for (auto it=(*newset)->first; it!=(*newset)->second; ++it)
        set_of[*it] = *newset; //update element -> set mapping
      set->first = mid;
      return newset;
    }
    return nullptr;
  }
};

}  // namespace nbautils
