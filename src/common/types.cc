#include <iostream>
#include <algorithm>
#include <stack>
#include "types.hh"

namespace nbautils {
using namespace std;

std::ostream& operator<<(std::ostream& os, ranked_slice const& rs) {
  for (int i=0; i<(int)rs.size(); ++i)
    os << pretty_bitset(rs[i].first) << ":" << rs[i].second << (i!=(int)rs.size()-1 ? ", " : "");
  return os;
}

std::ostream& operator<<(std::ostream& os, vector<ranked_slice> const& rss) {
  for (int i=0; i<(int)rss.size(); ++i)
    os << rss[i] << (i!=(int)rss.size()-1 ? " | " : "");
  return os;
}

std::ostream& operator<<(std::ostream& os, tree_history const& th) {
  for (int i=0; i<(int)th.size(); ++i)
    os << pretty_bitset(th[i]) << (i!=(int)th.size()-1 ? ", " : "");
  return os;
}

// ----

// takes a ranked slice, returns a "unpruned" tuple,
// i.e., labels contain states of whole subtree
// unpruned nodes in rank order uniquely determine a rank slice and are useful for
// storing ranked slices in k-equiv-aware lookup table (trie)
ranked_slice unprune(ranked_slice const& rslice) {
  int const n = rslice.size();
  ranked_slice ret(n);
  stack<pair<nba_bitset,pri_t>> s;
  for (int i=0; i<n; i++) {
    nba_bitset tmp = rslice[i].first;
    while (!s.empty() && rslice[i].second < s.top().second) {
      tmp |= s.top().first;
      s.pop();
    }
    auto const el = make_pair(tmp, rslice[i].second);
    ret[i] = el;
    s.push(el);
  }
  return ret;
}

// take unpruned tuple, reverse operation
ranked_slice prune(ranked_slice const& rslice) {
  int const n = rslice.size();
  ranked_slice ret(n);
  stack<pair<nba_bitset,pri_t>> s;
  for (int i=0; i<n; i++) {
    nba_bitset tmp = rslice[i].first;
    while (!s.empty() && rslice[i].second < s.top().second) {
      tmp &= ~s.top().first;
      s.pop();
    }
    ret[i] = make_pair(tmp, rslice[i].second);
    s.push(rslice[i]);
  }
  return ret;
}

//sort some ranked slice by rank, drop ranks.
//invertible if ranked slice was unpruned
tree_history in_rank_order(ranked_slice const& uprslice) {
  auto rng = uprslice;
  sort(rng.begin(), rng.end(), [](auto const &a, auto const &b){ return a.second < b.second; });
  tree_history ret;
  transform(rng.begin(), rng.end(), back_inserter(ret), [](auto const e){ return e.first; });
  return ret;
}

// given ranked slice, convert to tree history
tree_history slice_to_history(ranked_slice const& rs) {
  return in_rank_order(unprune(rs));
}

// given tree-history, convert to ranked slice
ranked_slice history_to_slice(tree_history const& th) {
  ranked_slice res(th.size());
  for (int i=0; i<(int)res.size(); ++i)
    res[i] = make_pair(th[i], i+1);

  // reconstruct correct sorting order from ranks and label/Safra properties
  // ( a subset b, b subset a or a cap b = emptyset )
  sort(res.begin(), res.end(), [](auto const &a, auto const &b){
    if ((a.first&~b.first)==0)
      return true;
    else if ((b.first&~a.first)==0)
      return false;
    else
      return a.second < b.second;
  });

  return prune(res);
}

// ----

//rs1 <= rs2 if rs2 consists of merged neighboring sets of rs1
//assumes valid ranked slices with same total state set
bool finer_or_equal(ranked_slice const& rs1, ranked_slice const& rs2) {
  if (rs1.size() == 0 && rs2.size() == 0)
    return true;
  if ((rs1.size()==0) != (rs2.size()==0))
    return false;

  nba_bitset pref1 = 0;
  nba_bitset pref2 = 0;
  int i=0;
  int j=0;
  while (j<(int)rs2.size()) {
    pref2 |= rs2[j].first;
    // cout << "j=" << j << " " << pretty_bitset(pref2) << endl;

    while (i<(int)rs1.size() && ((pref1|rs1[i].first) & ~pref2) == 0) {
      pref1 |= rs1[i].first;
      // cout << "i=" << i << " " << pretty_bitset(pref1) << endl;

      i++;
    }

    if (pref1 != pref2)
        return false;

    j++;
  }
  return ((i==(int)rs1.size()) == (j==(int)rs2.size()));
}

pair<nba_bitset,vector<nba_bitset>> kcut_mask(tree_history const& th, int k) {
  nba_bitset forbidden=0;
  vector<nba_bitset> masks;
  nba_bitset tmp=0;
  for (int i=k; i<(int)th.size(); ++i) {
    tmp |= th[i];
    masks.push_back(tmp);
  }
  forbidden = th[0] & ~tmp;
  return make_pair(forbidden, masks);
}

// k equiv if have same k prefix in trie branch
bool kequiv(tree_history const& th1, tree_history const& th2, int k) {
  if (k>=(int)th1.size() || k>=(int)th2.size())
    return false;
  for (int i=0; i<=k; ++i)
    if (th1[i] != th2[i])
      return false;
  return true;
}

// k cut not worse in t1 compared to t2
// if in trie nodes >= k never appear states of t2 that never appear >= k in t2
// and all states in t1 appear in nodes >= k not later than in t2
bool not_worse(tree_history const& th1, tree_history const& th2, int k) {
  if (!kequiv(th1, th2, k))
    return false;

  auto const msk = kcut_mask(th2, k);
  nba_bitset tmp = 0;

  for (int i=k; i<(int)th1.size(); ++i) {
    if ((th1[i] & msk.first) != 0)
      return false;

    tmp |= th1[i];
    if ((msk.second[i-k] & ~tmp) != 0)
      return false;
  }
  return true;
}

// ----

//returns Parent and Left-border relationship for a list of ranks
//(left border = left sibling for the last one popped in the inner while loop)
pair<vector<int>,vector<int>> unflatten(ranked_slice const& rank) {
  int const n = rank.size();
  pair<vector<int>,vector<int>> ret = make_pair(vector<int>(n),vector<int>(n));
  stack<int> s;
  for (int i=n-1; i>=0; i--) {
    while (!s.empty() && rank[i].second < rank[s.top()].second) {
      ret.second[s.top()] = i;
      s.pop();
    }
    ret.first[i] = s.empty() ? -1 : s.top();
    s.push(i);
  }
  while (!s.empty()) {
    ret.second[s.top()] = -1;
    s.pop();
  }
  return ret;
}


}
