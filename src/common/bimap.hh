#pragma once
#include <functional>
#include <map>
#include <unordered_map>
#include <memory>
#include <iostream>

namespace nbautils {
using namespace std;

//bimap interface, as required for graph node tagging
//K=tag, V=node id
template <typename K, typename V, class Impl>
struct bimap {
  Impl* impl() { return static_cast<Impl*>(this); }

  size_t size() { return impl()->size(); }

  bool has(K const& k) { return impl()->has(k); }
  bool hasi(V const& v) { return impl()->hasi(v); }

  V get(K const& k) { return impl()->get(k); }
  K geti(V const& v) { return impl()->geti(v); }

  V put_or_get(K const& k, V const& v) { return impl()->put_or_get(k, v); }
  void put(K const& k, V const& v) { impl()->put(k, v); }

  void erase(K const& k) { impl()->erase(k); }
  void erasei(V const& v) { impl()->erasei(v); }
};

// the concrete instances:

template <typename K, typename V, template <typename... Args> class M>
class naive_bimap : public bimap<K, V, naive_bimap<K, V, M>> {
  M<K, V> ktov;
  M<V, K> vtok;

 public:

  size_t size() const { return ktov.size(); }

  bool has(K const& k) const { return ktov.find(k) != end(ktov); }
  bool hasi(V const& v) const { return vtok.find(v) != end(vtok); }
  V get(K const& k) const { return ktov.at(k); }
  K geti(V const& v) const { return vtok.at(v); }

  //if key has other value, return existing value
  V put_or_get(K const& k, V const& v) {
    if (has(k)) return ktov.at(k);
    ktov[k] = v;
    vtok[v] = k;
    return v;
  }

  //if key or value associated otherwise, remove old link
  void put(K const& k, V const& v) {
    erasei(v);
    if (has(k)) erasei(get(k));
    ktov[k] = v;
    vtok[v] = k;
  }

  //remove value
  void erasei(V const& v) {
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

};  // namespace nbautils
