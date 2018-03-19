#pragma once
#include <functional>
#include <map>
#include <unordered_map>
#include <memory>
#include <iostream>
#include "triebimap.hh"

namespace nbautils {
using namespace std;

//bimap interface, as required for graph node tagging
//K=tag, V=node id
template <typename K, typename V, class Impl>
struct bimap {
  typedef unique_ptr<bimap<K, V, Impl>> ptr;

  Impl* impl() { return static_cast<Impl*>(this); }

  size_t size() { return impl()->size(); }
  bool has(K const& k) { return impl()->has(k); }
  bool hasi(V const& v) { return impl()->hasi(v); }
  V get(K const& k) { return impl()->get(k); }
  K geti(V const& v) { return impl()->geti(v); }
  V put_or_get(K const& k, V const& v) { return impl()->put_or_get(k, v); }
  void put(K const& k, V const& v) { impl()->put(k, v); }
  void erase(V const& v) { impl()->erase(v); }
};

// the concrete instances:

template <typename T, typename K, typename V>
class generic_trie_bimap : public bimap<T, V, generic_trie_bimap<T, K, V>> {
  typedef function<vector<K>(T const&)> fun_from;
  typedef function<T(vector<K> const&)> fun_to;

  fun_from serialize;
  fun_to deserialize;
  trie_bimap<K, V> m;

 public:
  generic_trie_bimap(fun_from from, fun_to to) : serialize(from), deserialize(to) {}

  size_t size() const { return m.size(); }
  bool has(T const& k) { return m.has(serialize(k)); }
  bool hasi(V const& v) { return m.has(v); }
  V get(T const& k) { return m.get(serialize(k)); }
  T geti(V const& v) { return deserialize(m.get(v)); }
  V put_or_get(T const& k, V const& v) { return m.put_or_get(serialize(k), v); }
  void put(T const& k, V const& v) { m.put(serialize(k), v); }
  void erase(V const& v) { m.erase(v); }
};

template <typename K, typename V, template <typename... Args> class M>
class naive_bimap : public bimap<K, V, naive_bimap<K, V, M>> {
  map<K, V> ktov;
  map<V, K> vtok;

 public:
  size_t size() const { return ktov.size(); }
  bool has(K const& k) const { return ktov.find(k) != end(ktov); }
  bool hasi(V const& v) const { return vtok.find(v) != end(vtok); }
  V get(K const& k) const { return ktov.at(k); }
  K geti(V const& v) const { return vtok.at(v); }

  V put_or_get(K const& k, V const& v) {
    if (has(k)) return ktov.at(k);
    ktov[k] = v;
    vtok[v] = k;
    return v;
  }
  void put(K const& k, V const& v) {
    erase(v);
    if (has(k)) erase(get(k));
    ktov[k] = v;
    vtok[v] = k;
  }
  void erase(V const& v) {
    if (!hasi(v)) return;
    auto const k = geti(v);
    ktov.erase(ktov.find(k));
    vtok.erase(vtok.find(v));
  }
};

template <typename K, typename V>
using naive_ordered_bimap = naive_bimap<K,V,map>;
template <typename K, typename V>
using naive_unordered_bimap = naive_bimap<K,V,unordered_map>;

//TODO: add boost bimap as possibility?

};  // namespace nbautils
