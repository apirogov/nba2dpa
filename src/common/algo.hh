#pragma once

#include <iostream>
#include "common/util.hh"

#include <functional>
#include <queue>
#include <stack>
#include <set>
#include <map>

namespace nbautils {
  using namespace std;
  using namespace nbautils;

  template <typename T>
  using succ_fun = function<vector<T>(T)>;

  using succ_scc_fun = succ_fun<unsigned>;

  template <typename T,typename S>
  using outsym_fun = function<vector<S>(T)>;

  template <typename T,typename S>
  using succ_sym_fun = function<vector<T>(T,S)>;


// datatype generic graph algorithms
// (the price is that all required operations must be passed as lambdas explicitly)

//generic bfs. input: start node, function that takes current node,
//a function to schedule a visit and a visited and discovery check
//the visit function just does whatever needed with current node and calls
//pusher function on all successors that also need to be visited.
//can use visited function to check for already visited states
//can use discovered function to check for states already in the visit pipeline
//bfs keeps track that each node is visited once in bfs order automatically.
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

// ----------------------------------------------------------------------------

template <typename Node>
struct SCCDat {
  std::vector<std::vector<Node>> sccs;
  std::map<Node, unsigned> scc_of;
};

// https://en.wikipedia.org/wiki/Path-based_strong_component_algorithm with extensions
// takes list of all states of graph we want to have an scc for
// a function that supplies successors of a state
// a function that is called for each new found SCC and returns whether search has to go on
// performs an SCC DFS traversal.
// returns list of SCCs such that later SCCs can not reach earlier SCCs
template <typename T, typename F, typename G>
SCCDat<T> get_sccs(std::vector<T> const& states, F get_succs, G scc_handler) {
  SCCDat<T> ret;

  std::stack<T> call;  // dfs call stack
  std::stack<T> reps;  // scc representative stack
  std::stack<T> open;  // not yet fully completed vertex stack

  int count = 0;
  std::map<T, unsigned> order;  // first visit order

  // schedule all states to be called
  for (auto const& v : states)
    call.push(v);

  while (!call.empty()) {
    auto const v = call.top();

    if (!map_has_key(order,v)) {  // dfs just "called" with current node
      // std::cout << "discover " << v << std::endl;

      order[v] = count++;         // assign visit order
      reps.push(v);  // SCC representative candidate (for now)
      open.push(v);  // this node is "pending" (not completely discovered from here)

      auto const sucs = get_succs(v);

      for (auto const w : sucs) {  // process edges with any label
        if (!map_has_key(order, w)) {
          call.push(w);  // recursively explore nodes that have not been visited yet
        } else if (!map_has_key(ret.scc_of, w)) {
          // if already visited, but not with assigned scc, we have found a loop
          // -> drop candidates, keep oldest on this loop as SCC representative
          while (order.at(reps.top()) > order.at(w)) reps.pop();
        }
      }

    } else {
      if (map_has_key(ret.scc_of, v)) {
        //this node is already completed and uselessly visited
        call.pop();
        continue;
      }

      // returned from recursive calls
      // std::cout << "return to " << v << std::endl;

      // is still rep. -> we found an SCC
      if (reps.top() == v) {
        // drop states up to the current state, they are done and part of the SCC
        std::vector<T> scc_states;
        T tmp;
        do {
          tmp = open.top();
          open.pop();
          ret.scc_of[tmp] = ret.sccs.size();
          scc_states.push_back(tmp);
        } while (tmp != v);

        ret.sccs.push_back({});
        swap(ret.sccs.back(), scc_states);

        bool go_on = scc_handler(ret);
        if (!go_on) //can be used to abort on-the-fly scc search
          return ret;

        // current SCC is done
        reps.pop();
      }

      // std::cout << "done " << v << std::endl;
      call.pop();  // current node done -> dfs "returns" to previous caller
    }
  }

  return ret;
}

// takes SCC information and an SCC, returns successor SCCs
template <typename T>
std::vector<unsigned> succ_sccs(SCCDat<T> const& scci, unsigned const& num,
    succ_fun<T> const& get_suc_st) {
  std::set<unsigned> sucsccs;
  bfs(scci.sccs.at(num).front(), [&](auto const& st, auto const& visit, auto const&) {
    for (auto const& sucst : get_suc_st(st)) {
      auto sucscc = scci.scc_of.at(sucst);
      if (sucscc == num)
        visit(sucst);
      else
        sucsccs.emplace(sucscc);
    }
  });
  return vector<unsigned>(cbegin(sucsccs), cend(sucsccs));
}

template <typename T>
std::set<unsigned> trivial_sccs(SCCDat<T> const& scci, succ_fun<T> const& get_succs) {
  set<unsigned> ret;
  unsigned i=0;
  for (auto const& scc : scci.sccs) {
    // single state with no self-loop?
    bool const trivacc = scc.size() == 1;
    bool const noselfloop = !contains(get_succs(scc.front()), scc.front());
    if (trivacc && noselfloop)
      ret.emplace(i);

    ++i;
  }
  return ret;
}

// given sccs and map of states to semigroup together with semigroup operation
// and default value, return map from scc to semigroup element
// (used to get accepting, rejecting sccs of NBA and dominant SCC priority of DPAs)
template <typename T, typename V>
map<unsigned, V> fold_sccs(SCCDat<T> const& scci,
    V def, function<V(T)> lift, function<V(V,V)> op) {
  map<unsigned, V> ret;
  unsigned i=0;
  for (auto const& scc : scci.sccs) {
    V accum = def;
    for (auto v : scc) {
      accum = op(accum, lift(v));
    }
    ret[i] = accum;
    ++i;
  }
  return ret;
}

// ----------------------------------------------------------------------------

struct BaSccClassification {
  set<unsigned> accepting;  // tags an scc as fully accepting
  set<unsigned> rejecting;  // tags an scc as fully rejecting
};

// classify sccs as accepting or rejecting
template <typename T>
BaSccClassification ba_classify_sccs(SCCDat<T> const& scci,
    function<bool(T)> const& is_acc) {
  BaSccClassification ret;

  auto const conj = [](bool a, bool b){return a&&b;};
  ret.accepting = mapbool_to_set(fold_sccs<T,bool>(scci, true, is_acc, conj));
  ret.rejecting = mapbool_to_set(fold_sccs<T,bool>(scci, true,
        [&is_acc](T v){ return !is_acc(v); }, conj));

  return ret;
}

//returns sorted list of accepting sinks (acc. states with self-loop for each sym)
template <typename Node, typename Sym>
vector<Node> ba_get_accepting_sinks(vector<Node> const& states, int num_syms,
    function<bool(Node)> is_acc,
    outsym_fun<Node,Sym> get_outsyms,
    succ_sym_fun<Node,Sym> get_xsuccs) {
  vector<Node> ret;
  for (auto const& v : states) {
    auto outsyms = get_outsyms(v);
    //must be accepting and have successors for each symbol
    bool accsink = is_acc(v) && (int)outsyms.size() == num_syms;
    for (auto i=0; i<num_syms; i++) {
        // must have self-loop for each symbol
        auto xsucs = get_xsuccs(v, i);
        if (!contains(xsucs, v))
          accsink = false;
    }
    if (accsink)
      ret.push_back(v);
  }
  return ret;
}

inline void mark_dead_sccs(set<unsigned> const& rejecting, set<unsigned> const& trivial,
    succ_scc_fun const& get_suc_sccs, map<unsigned,bool>& dead, unsigned num) {
  if (map_has_key(dead, num))  // done already
    return;

  // std::cout << "scc " << num << std::endl;
  auto sucsccs = get_suc_sccs(num);
  // mark children first
  for (auto sucnum : sucsccs)
    mark_dead_sccs(rejecting, trivial, get_suc_sccs, dead, sucnum);

  // if we are rejecting and trivial, assume we're dead
  bool isdead = contains(rejecting, num) || contains(trivial, num);
  // check children and try to falsify
  for (auto sucscc : sucsccs) isdead = isdead && dead[sucscc];
  // if still dead, we're really dead.
  dead[num] = isdead;
}

// run dfs that marks dead sccs (assuming sccs 0..num_sccs-1 exist)
inline set<unsigned> ba_get_dead_sccs(int num_sccs,
    set<unsigned> const& rejecting, set<unsigned> const& trivial,
    succ_scc_fun const& get_succs) {
  map<unsigned, bool> dead;
  for (int i=0; i<num_sccs; i++)
    mark_dead_sccs(rejecting, trivial, get_succs, dead, i);

  set<unsigned> ret;
  for (auto const& it : dead)
    if (it.second)
      ret.emplace(it.first);
  return ret;
}


// ----------------------------------------------------------------------------

template <typename T>
unsigned max_chain(function<unsigned(T)> const& oldpri, map<unsigned,unsigned>& newpri,
    vector<T> const& p, succ_fun<T> const& get_succs) {
  unsigned maxlen = 0;
  // cout << "max_chain " << seq_to_str(p) << endl;

  // maximal essential subsets = non-trivial SCCs in restricted graph
  succ_fun<T> succs_in_p = [&](auto v) { return set_intersect(get_succs(v), p); };
  auto scci = get_sccs(p, succs_in_p, const_true);
  auto triv = trivial_sccs(scci, get_succs);

  // lift priority map to state sets
  auto strongest_oldpri = fold_sccs<T,unsigned>(scci, 0, oldpri, [](auto a, auto b){return max(a,b);});

  unsigned i=0;
  for (auto const& scc : scci.sccs) {
    // cout << "scc " << seq_to_str(scc) << endl;

    if (contains(triv,i)) {
      ++i;
      continue; //non-essential
    }

    auto const scc_pri = strongest_oldpri.at(i);

    // cout << "scc pri " << scc_pri << endl;

    auto deriv_scc = vec_filter(scc, [&](auto v){return oldpri(v) < scc_pri;});
    vec_to_set(deriv_scc);
    auto const not_deriv_scc = vec_filter(scc, [&](auto v){return oldpri(v) >= scc_pri;});

    unsigned m = 0;
    if (!deriv_scc.empty() && scc_pri > 0) {
      m = max_chain(oldpri, newpri, deriv_scc, succs_in_p);

      if ((scc_pri - m) % 2 == 1) //alternation -> requires new priority
        m++;
    }

    for (auto s : not_deriv_scc) {
      newpri[oldpri(s)] = m;
    }

    maxlen = max(maxlen, m);

    ++i;
  }
  // cout << "max_chain " << seq_to_str(p)  << " done" << endl;
  return maxlen;
}

// takes max odd priorities (max odd parity), minimizes them
// "Computing the Rabin Index of a parity automaton"
// takes all states, normal successor function and old priority map function
// priorities must be from a max odd acceptance
// returns new priority map (old priority to new)
template <typename T>
map<T, unsigned> minimize_priorities(vector<T> const& states,
    succ_fun<T> const& get_succs, function<unsigned(T)> const& oldpri) {
  map<unsigned, unsigned> primap;
  for (auto const s : states)
    primap[oldpri(s)] = 0;
  max_chain(oldpri, primap, states, get_succs);
  return primap;
}

}
