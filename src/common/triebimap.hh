#pragma once

#include "util.hh" //map_get_keys
#include <map>
#include <memory>
#include <vector>

namespace nbautils {

template <typename K, typename V>
struct trie_bimap_node {
  using node_ptr = std::unique_ptr<trie_bimap_node<K, V>>;
  trie_bimap_node* parent = nullptr;
  K key = 0;
  std::unique_ptr<V> value = nullptr;

  std::map<K, node_ptr> suc;
};

template <typename K, typename V>
class trie_bimap {
  using node_t = trie_bimap_node<K, V>;
  node_t root;
  std::map<V,node_t*> revmap;

  // traverse trie to given node and return it. when create=false and it does
  // not exist, return null pointer
  node_t* traverse(std::vector<K> const &ks, bool create = false) {
    auto *curr = &root;
    for (int i = ks.size() - 1; i >= 0; i--) {
      if (curr->suc.find(ks[i]) == curr->suc.end()) {
        if (create) {
          curr->suc[ks[i]] = std::make_unique<node_t>();
          curr->suc[ks[i]]->parent = curr;
          curr->suc[ks[i]]->key = ks[i];
        } else
          return nullptr;
      }
      curr = curr->suc[ks[i]].get();
    }
    return curr;
  }

 public:
  size_t size() const { return revmap.size(); }

  std::vector<V> values() const { return map_get_keys(revmap); }

  // puts a set,value pair
  void put(std::vector<K> const &ks, V val) {
    auto curr = traverse(ks, true);
    if (curr->value) revmap.erase(revmap.find(*(curr->value)));
    curr->value = std::make_unique<V>(val);
    revmap[val] = curr;
  }

  V put_or_get(std::vector<K> const &ks, V val) {
    auto *curr = traverse(ks, true);
    if (!curr->value) {
      curr->value = std::make_unique<V>(val);
      revmap[val] = curr;
    }
    return *(curr->value);
  }

  bool has(std::vector<K> const &ks) {
    auto *curr = traverse(ks, false);
    if (!curr || !curr->value) return false;
    return true;
  }

  bool has(V val) { return revmap.find(val) != revmap.end(); }

  V get(std::vector<K> const &ks) {
    auto const *curr = traverse(ks, false);
    return *(curr->value);
  }

  std::vector<K> get(V val) {
    auto const *curr = revmap[val];
    std::vector<K> ret;
    while (curr->parent) {
      ret.push_back(curr->key);
      curr = curr->parent;
    }
    return ret;
  }

  //TODO: clean up branches when no leaves are below? (if efficiently possible)
  void erase(V val) {
    auto it = revmap.find(val);
    if (it == revmap.end())
      return;
    node_t* curr = it->second;
    revmap.erase(it);
    curr->value = nullptr;
  }

};
}  // namespace nbautils
