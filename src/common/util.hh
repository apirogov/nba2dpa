#pragma once

#include <algorithm>
#include <numeric>
#include <functional>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <sstream>
#include <cassert>

auto identity = [](auto const& t){ return t; };

auto const_true = [](auto const&){ return true; };

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

// vector fmap
template<typename A, typename F>
auto vec_fmap(std::vector<A> const& v, F const& f) {
  std::vector<decltype(f(v.front()))> ret;
  transform(std::cbegin(v), std::end(v), std::back_inserter(ret), f);
  return ret;
}

template<typename A, typename B>
B vec_fold_mapped(std::vector<A> const& v, B const& def,
  std::function<B(A)> lift, std::function<B(B,B)> const& op) {
  return std::accumulate(cbegin(v), cend(v), def,
      [&](B const& accum, A const& val){ return op(accum, lift(val)); });
}

template<typename T, typename F>
std::vector<T> vec_filter(std::vector<T> const& v, F const& f) {
  std::vector<T> ret;
  copy_if(std::cbegin(v), std::cend(v), std::back_inserter(ret), f);
  return ret;
}

//TODO: make one abstract for these

//TODO: make set_intersect_empty without full calc.
//and replace all set_intersect(..).empty()

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
template <typename K>
inline std::set<K> mapbool_to_set(std::map<K,bool> const& m) {
  std::set<K> ret;
  for (auto &it : m)
    if (it.second)
      ret.emplace(it.first);
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

//TODO: group_by for quotienting
template <typename T, typename F>
std::vector<std::vector<T>> group_by(std::vector<T> v, F const& f) {
  if (v.empty())
      return {};
}
