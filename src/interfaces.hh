#pragma once
#include <functional>
#include <map>
#include <memory>

#include "triebimap.hh"

namespace nbautils {
using namespace std;

template <typename K, typename V, class Impl>
struct bimap {
  typedef unique_ptr<bimap<K, V, Impl>> ptr;

  size_t size() const { return static_cast<Impl const*>(this)->size(); }

  bool has(K const& k) { return static_cast<Impl*>(this)->has(k); }
  bool has(V const& v) { return static_cast<Impl*>(this)->has(v); }
  V get(K const& k) { return static_cast<Impl*>(this)->get(k); }
  K get(V const& v) { return static_cast<Impl*>(this)->get(v); }
  V put_or_get(K const& k, V const& v) { static_cast<Impl*>(this)->put_or_get(k, v); }
  void put(K const& k, V const& v) { static_cast<Impl*>(this)->put(k, v); }
  void erase(V const& v) const { static_cast<Impl*>(this)->erase(v); }
};

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
  bool has(V const& v) { return m.has(v); }
  V get(T const& k) { return m.get(serialize(k)); }
  T get(V const& v) { return deserialize(m.get(v)); }
  V put_or_get(T const& k, V const& v) { return m.put_or_get(serialize(k), v); }
  void put(T const& k, V const& v) { m.put(serialize(k), v); }
  void erase(V const& v) { /* not implemented */
  }
};

template <typename K, typename V>
class naive_bimap : public bimap<K, V, naive_bimap<K, V>> {
  map<K, V> ktov;
  map<V, K> vtok;

 public:
  size_t size() const { return ktov.size(); }
  bool has(K const& k) const { return ktov.find(k) != end(ktov); }
  bool has(V const& v) const { return vtok.find(v) != end(vtok); }
  V get(K const& k) const { return ktov.at(k); }
  K get(V const& v) const { return vtok.at(v); }

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
    if (!has(v)) return;
    auto const k = get(v);
    ktov.erase(ktov.find(k));
    vtok.erase(vtok.find(v));
  }
};

};  // namespace nbautils
