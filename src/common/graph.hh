#pragma once

#include "common/util.hh"

#include <functional>
#include <queue>
#include <set>

namespace nbautils {

using namespace std;
using namespace nbautils;

template <typename Node>
using succ_fun = function<vector<Node>(Node)>;

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

template <typename Node>
vector<Node> reachable_states(Node from, succ_fun<Node> get_succ) {
  std::set<Node> reached;
  bfs(from, [&](Node const& st, auto const& pusher, auto const&) {
      reached.emplace(st);
      for (auto sucst : get_succ(st))
        pusher(sucst);
  });
  return vector<Node>(cbegin(reached), cend(reached));
}

template <typename Node>
vector<Node> unreachable_states(std::vector<Node> const& allstates, Node from,
    succ_fun<Node> get_succ) {
  assert(is_set_vec(allstates));
  return set_diff(allstates, reachable_states(from, get_succ));
}

}
