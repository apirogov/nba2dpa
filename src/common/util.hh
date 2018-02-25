#pragma once

#include <algorithm>
#include <functional>
#include <queue>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <sstream>
#include <cassert>

template<typename T>
T identity(T t){ return t; };

// is sorted + unique vector?
template<typename T>
bool is_set_vec(std::vector<T> const& v) {
  std::set<T> s(std::cbegin(v),std::cend(v));
  return std::equal(std::cbegin(v), std::cend(v), std::cbegin(s));
}

template<typename T>
std::string seq_to_str(T const& s, std::string const& sep=",") {
  std::stringstream ss;
  auto last = --std::end(s);
  for (auto it = std::cbegin(s); it!=std::cend(s); ++it) {
    ss << std::to_string(*it);
    if (it != last) {
      ss << sep;
    }
  }
  return ss.str();
}

// sort + make unique inplace
template <typename T>
void vec_to_set(std::vector<T>& v) {
  std::sort(std::begin(v),std::end(v));
  v.erase(std::unique(std::begin(v),std::end(v)), std::end(v));
}

// monomorphic inplace fmap
// template<template<typename> typename C, typename T, typename F>
// inline void fmap_inplace(C<T>& v, F const& w) {
//   transform(begin(v), end(v), begin(v), w);
// }
// template<template<typename> typename C, typename T, typename F>
// inline C<T> fmap(C<T>& v, F const& w) {
//   auto ret = v;
//   transform(begin(ret), end(ret), begin(ret), w);
//   return ret;
// }

//TODO: make one abstract for these

//TODO: make set_intersect_empty without full calc.
//and replace all set_intersect(..).empty()

//TODO: group_by for quotienting

template<typename T>
inline std::vector<T> set_intersect(std::vector<T> const& v,
                                    std::vector<T> const& w) {
  assert(is_set_vec(v)); assert(is_set_vec(w));
  std::vector<T> ret;
  std::set_intersection(std::cbegin(v), std::cend(v),
                        std::cbegin(w), std::cend(w), std::back_inserter(ret));
  return ret;
}
template<typename T>
inline std::vector<T> set_merge(std::vector<T> const& v,
                                std::vector<T> const& w) {
  assert(is_set_vec(v)); assert(is_set_vec(w));
  std::vector<T> ret;
  std::set_union(std::cbegin(v), std::cend(v),
                 std::cbegin(w), std::cend(w), std::back_inserter(ret));
  return ret;
}
template<typename T>
inline std::vector<T> set_diff(std::vector<T> const& v,
                               std::vector<T> const& w) {
  assert(is_set_vec(v)); assert(is_set_vec(w));
  std::vector<T> ret;
  std::set_difference(std::cbegin(v), std::cend(v),
                      std::cbegin(w), std::cend(w), std::back_inserter(ret));
  return ret;
}


// collect all keys in a map, returns sorted vector
template <typename K, typename V>
std::vector<K> map_get_keys(std::map<K,V> const& m) {
  std::vector<K> ret;
  for (auto &it : m)
    ret.push_back(it.first);
  return ret;
}

template <typename K, typename V>
std::vector<K> map_get_vals(std::map<K,V> const& m) {
  std::vector<K> ret;
  for (auto &it : m)
    ret.push_back(it.second);
  return ret;
}

// check whether a map has a key using .find()
template <typename K, typename V>
inline bool map_has_key(std::map<K,V> const& m, K const& k) {
  return m.find(k) != end(m);
}

// check whether a container has a value using .find()
template <typename C, typename V>
inline bool contains(C const& c, V const& v) {
  return std::find(std::cbegin(c), std::cend(c), v) != std::cend(c);
}

//generic bfs. input: start node, function that takes current node and a pusher function
//the visit function just does whatever needed with current node and calls
//pusher function on all successors that also need to be visited.
//bfs keeps track that each node is visited once in bfs order automatically.
template <typename T, typename F>
void bfs(T const& start, F visit) {
  std::queue<T> bfsq;
  std::set<T> visited;

  auto pusher = [&](T const& el){ bfsq.push(el); };
  auto checker = [&](T const& el){ return contains(visited, el); };
  pusher(start);

  while (!bfsq.empty()) {
    auto const st = bfsq.front();
    bfsq.pop();
    if (visited.find(st) != visited.end()) continue;  // have visited this one
    visited.emplace(st);

    visit(st, pusher, checker);
  }
}
