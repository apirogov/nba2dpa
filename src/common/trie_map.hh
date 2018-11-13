#pragma once

#include <unordered_map>
#include <memory>
#include <vector>

namespace nbautils {

using namespace std;

template <typename K, typename V>
struct trie_node {
  using node_ptr = unique_ptr<trie_node<K, V>>;

  // trie_node* parent = nullptr;
  K key = 0;
  unique_ptr<V> value = nullptr;

  unordered_map<K, node_ptr> suc;
};

template <typename K, typename V>
class trie_map {
  using node_t = trie_node<K, V>;
  node_t root;
  size_t sz = 0;

 public:

  // traverse trie to given node and return it.
  // when create=false and it does not exist, return null pointer
  node_t* traverse(vector<K> const &ks, bool create = false) {
    // if (k<0) k=ks.size();
    auto *curr = &root;
    for (int i = 0; i < (int)ks.size(); i++) {
      if (curr->suc.find(ks[i]) == curr->suc.end()) {
        if (create) {
          curr->suc[ks[i]] = make_unique<node_t>();
          // curr->suc[ks[i]]->parent = curr;
          curr->suc[ks[i]]->key = ks[i];
        } else
          return nullptr;
      }
      curr = curr->suc[ks[i]].get();
    }
    return curr;
  }

  size_t size() const { return sz; }

  // puts a set,value pair
  void put(vector<K> const &ks, V val) {
    auto curr = traverse(ks, true);
    if (curr->value)
      --sz;
    curr->value = make_unique<V>(val);
    ++sz;
  }

  bool has(vector<K> const &ks, bool subtree=false) {
    auto const *curr = traverse(ks, false);
    if (!curr)
      return false;
    if (!subtree)
      return curr->value!=nullptr;
    else
      return curr->value!=nullptr || curr->suc.size()>0;
  }

  V get(vector<K> const &ks, bool any=false) {
    auto *curr = traverse(ks, false);
    if (any) {
      while (!curr->value) {
        curr = curr->suc.begin()->second.get();
      }
    }
    return *(curr->value);
  }

};

}  // namespace nbautils
