#pragma once

#include <vector>
#include <map>
#include <set>
#include <algorithm>

//input: list of sets. output: a hitting set
//(i.e., set containing at least one element from each set in the list)
template <typename A>
std::set<A> greedy_hitting_set(std::map<A,std::set<A>> const& s2u) {
  // auto s2u = sets; //copy, as we will modify it
  // now we need the reverse mapping from "universe" to sets
  // (we conceptually have a bipartite graph)

  //map elements of universe to sets they appear in
  std::map<A,std::set<A>> u2s;
  for (auto const& s : s2u) {
    for (auto const e : s.second) {
      u2s[e].emplace(s.first);
    }
  }

  std::set<A> h; //hitting set

  //each loop iteration starts with remaining uncovered sets on the left as reachable sets from right
  //and useful hitting set candidates on the right which reach at least one uncovered set on left
  //after an iteration, covered sets are unreachable and isolated elements on the right removed.
  //invariant: el has edge to set -> set has edge to el, and every remaining el has an edge to some set
  while (!u2s.empty()) {
    //sort candidates by how many yet uncovered sets they would cover
    std::vector<A> cands;
    for (auto const& it : u2s)
      cands.push_back(it.first);

    std::sort(std::begin(cands), std::end(cands), [&u2s](A a, A b){
        auto const sza = u2s.at(a).size();
        auto const szb = u2s.at(b).size();
        //sort desc. by # of touched sets, then asc. by val
        return  sza > szb || (sza == szb && a < b);
    });

    //take greedily a candidate that covers as many uncovered sets as possible
    auto const newel = cands.front();
    h.emplace(newel); //add it to hitting set
    auto const curcovered = u2s.at(newel); //get newly covered sets

    // cerr << "picked " << newel << " freshly covering " << seq_to_str(curcovered) << endl;

    //make each covered set remove itself from all its elements, i.e. kill incoming edges
    //this makes those sets unreachable from right side
    for (auto const sidx : curcovered) {
      for (auto const el : s2u.at(sidx)) {
        // the set sidx was uncovered before, i.e. all its elements are still candidates
        // hence there are no el in s2u[sidx] that are not in u2s anymore
        // those candidates must also still have an edge back to set sidx
        // hence no need to check for null/end
        u2s.at(el).erase(u2s.at(el).find(sidx));
      }
    }

    //conceptually, we could remove covered sets
    //practically, it suffices that they are unreachable from the elements now
    /*
    for (auto const sidx : curcovered) {
      auto it = s2u.find(sidx);
      if (it != end(s2u))
        s2u.erase(it);
    }
    */

    //remove the now useless elements (not in any remaining uncovered set,
    //i.e. isolated nodes without outgoing edges),
    //including the new hitting set element, from graph
    auto it = std::begin(u2s);
    while (it != std::end(u2s)) {
      auto const tmp = it++;
      if (tmp->second.empty()) {
        // cerr << "remove now useless " << it->first << endl;
        u2s.erase(tmp);
      }
    }
  }

  return h;
}

