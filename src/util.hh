#pragma once

#include <functional>
#include <queue>
#include <vector>
#include <map>
#include <set>

struct Unit { bool operator==(Unit const&) const { return true; } };

template <typename K, typename V>
std::vector<K> map_get_keys(std::map<K,V> const& m) {
  std::vector<K> ret;
  for (auto &it : m)
    ret.push_back(it.first);
  return ret;
}

template <typename K, typename V>
inline bool map_has_key(std::map<K,V> const& m, K const& k) {
  return m.find(k) != end(m);
}

template <typename C, typename V>
inline bool contains(C const& c, V const& v) {
  return c.find(v) != end(c);
}

template <typename T, typename F>
void bfs(T const& start, F visit) {
  std::queue<T> bfsq;
  std::set<T> visited;

  auto pusher = [&](T const& el){ bfsq.push(el); };
  pusher(start);

  while (!bfsq.empty()) {
    auto const st = bfsq.front();
    bfsq.pop();
    if (visited.find(st) != visited.end()) continue;  // have visited this one
    visited.emplace(st);

    visit(st, pusher);
  }
}
