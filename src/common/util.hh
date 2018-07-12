#pragma once

#include <algorithm>
#include <iostream>
#include <numeric>
#include <functional>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <string>
#include <sstream>
#include <cassert>
#include <bitset>

// gives identity function for any type
auto identity = [](auto const& t){ return t; };

// gives function that returns fixed value for arbitrary parameter
// template <typename T, T V>
// auto return_const = [](auto const&){ return V; };
auto const_true = [](auto const&){ return true; };

// useful to compare on some mapping function, e.g.
// sort(some_seq, on(less<int>(), [](auto x){return -x;}));
auto on = [](auto f, auto g){ return [f, g](auto a, auto b){ return f(g(a),g(b)); }; };

// is sorted + unique vector?
template<typename T>
bool is_set_vec(std::vector<T> const& v) {
  std::vector<T> s(std::cbegin(v),std::cend(v));
  sort(std::begin(s), std::end(s));
  return std::equal(std::cbegin(v), std::cend(v), std::cbegin(s));
}

// sort + make unique inplace
template <typename T>
void vec_to_set(std::vector<T>& v) {
  std::sort(std::begin(v),std::end(v));
  v.erase(std::unique(std::begin(v),std::end(v)), std::end(v));
}

// arbitrary sequence (with iterators) to string, intercalated with separator
template<typename T>
std::string seq_to_str(T const& s, std::string const& sep=",") {
  std::stringstream ss;
  auto last = --std::end(s);
  for (auto it = std::cbegin(s); it!=std::cend(s); ++it) {
    ss << *it;
    if (it != last) {
      ss << sep;
    }
  }
  return ss.str();
}

//take sequence containing "small numbers", return small set
//numbers must be compatible with maximum size!
template<typename B, typename I>
B to_bitset(I s) {
  B ret = 0;
  for (auto const it : s)
    ret[it] = 1;
  return ret;
}

//takes a bitset and a container inserter. stores bits set to 1
template<typename B, typename O>
void from_bitset(B const& s, O out) {
  for (int i=0; i<(int)s.size(); i++) {
    if (s[i])
      *out = i;
    ++out;
  }
}

template<typename B, typename T>
std::vector<T> vec_from_bitset(B const& s) {
  std::vector<T> tmp;
  from_bitset<B>(s, back_inserter(tmp));
  return tmp;
}

template<typename B>
std::string pretty_bitset(B const& s) {
  return "{"+seq_to_str(vec_from_bitset<B,unsigned>(s))+"}";
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

template<typename T>
inline bool set_intersect_empty(std::vector<T> const& v,
                                          std::vector<T> const& w) {
  // assert(is_set_vec(v)); assert(is_set_vec(w));
  auto i = v.cbegin();
  auto j = w.cbegin();
  while (i != v.cend() && j != w.cend()) {
    if (*i == *j)
      return false;
    else if (*i < *j)
      ++i;
    else
      ++j;
  }
  return true;
}

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

/*
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
*/

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
template <template<typename,typename,typename...> class M, typename K, typename V, typename ... A>
inline bool map_has_key(M<K,V,A...> const& m, K const& k) {
  return m.find(k) != end(m);
}

// check whether a container has a value using .find()
template <typename C, typename V>
inline bool contains(C const& c, V const& v) {
  return std::find(std::cbegin(c), std::cend(c), v) != std::cend(c);
}

// use binary search to check containment
// (for random access iterators with ineff. find)
template <typename C, typename V>
inline bool sorted_contains(C const& c, V const& v) {
  auto it = std::lower_bound(std::cbegin(c), std::cend(c), v);
  return it != std::cend(c) && *it == v;
}

//group_by for quotienting
template <typename T, typename F>
std::vector<std::vector<T>> group_by(std::vector<T> v, F const& f) {
  if (v.empty())
    return {};
  if (v.size()==1)
    return {v};

  std::vector<std::vector<T>> ret;
  auto l=cbegin(v);
  ret.push_back({*l});
  auto r=l; ++r;
  do {
    if (!f(*l,*r))
      ret.push_back({});
    ret.back().push_back(*r);
    ++l; ++r;
  } while (r != cend(v));
  return ret;
}

//generic bfs. input: start node,
//function that takes current node,
//a function to schedule a visit
//a has_visited and has_discovered predicate
//the visit function just does whatever needed with current node and calls
//pusher function on all successors that also need to be visited.
//bfs keeps track that each node is visited once in bfs order automatically.
//TODO: maybe make something with "process_edge, give_edges" ?
template <typename Node, typename F>
void bfs(Node const& start, F visit) {
  std::queue<Node> bfsq;
  std::set<Node> visited;
  std::set<Node> discovered;

  auto pusher = [&](Node const& st){
    if (!contains(discovered, st)) {
      discovered.emplace(st);
      bfsq.push(st);
    }
  };
  auto visited_f = [&](Node const& el){ return contains(visited, el); };
  // auto discovered_f = [&](T const& el){ return contains(visited, el); };

  pusher(start);
  while (!bfsq.empty()) {
    auto const st = bfsq.front();
    bfsq.pop();
    if (visited.find(st) != visited.end()) continue;  // have visited this one
    visited.emplace(st);

    visit(st, pusher, visited_f /*, discovered_f */);
  }
}

template<typename T>
auto default_printer() {
  return [](std::ostream& out, T const& val){ out << val; };
}
