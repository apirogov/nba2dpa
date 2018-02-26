#include <iostream>
#include "io.hh"
#include "common/util.hh"

using namespace std;
using namespace nbautils;

// datatype generic graph algorithms

template <typename T>
struct SCCDat {
  std::vector<std::vector<T>> sccs;
  std::map<T, unsigned> scc_of;
};

// https://en.wikipedia.org/wiki/Path-based_strong_component_algorithm with extensions
// takes list of all states of graph we want to have an scc for
// a function that supplies successors of a state
// a function that takes a resulting SCC
// performs an SCC DFS traversal.
// returns list of SCCs such that later SCCs can not reach earlier SCCs
template <typename T, typename F>
SCCDat<T> get_sccs(std::vector<T> const& states, F get_succs) {
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

        // current SCC is done
        reps.pop();
      }

      // std::cout << "done " << v << std::endl;
      call.pop();  // current node done -> dfs "returns" to previous caller
    }
  }

  return ret;
}

struct BaSccClassification {
  set<unsigned> accepting;  // tags an scc as fully accepting
  set<unsigned> rejecting;  // tags an scc as fully rejecting
  set<unsigned> trivial;    // tags an scc as trivial (single state, no loop)
};

// classify sccs as accepting, rejecting or trivial
template <typename T>
BaSccClassification classify_sccs(SCCDat<T> const& scci,
    function<vector<T>(T)> const& get_succs, function<bool(T)> const& is_acc) {
  BaSccClassification ret;

  unsigned i=0;
  for (auto const& scc : scci.sccs) {
    // assume SCC is accepting and rejecting, then falsify
    bool accscc = true;
    bool rejscc = true;
    // identify acc./rej. sccs
    for (auto v : scc) {
      auto actmp = is_acc(v);
      accscc = accscc && actmp;
      rejscc = rejscc && !actmp;
    }

    // mark completed SCC as accepting/rejecting
    if (accscc)
      ret.accepting.emplace(i);
    if (rejscc)
      ret.rejecting.emplace(i);

    // trivial scc with no self-loop?
    bool const trivacc = scc.size() == 1;
    bool const noselfloop = !contains(get_succs(scc.front()), scc.front());
    if (trivacc && noselfloop)
      ret.trivial.emplace(i);

    ++i;
  }

  return ret;
}

template <typename T>
std::vector<unsigned> succ_sccs(SCCDat<T> const& scci, unsigned const& num,
    function<vector<T>(T)> const& get_suc_st) {
  std::set<nbautils::scc_t> sucsccs;
  bfs(scci.sccs.at(num).front(), [&](auto const& st, auto const& visit, auto const&) {
    for (auto const& sucst : get_suc_st(st)) {
      auto sucscc = scci.scc_of.at(sucst);
      if (sucscc == num)
        visit(sucst);
      else
        sucsccs.emplace(sucscc);
    }
  });
  return vector<nbautils::scc_t>(cbegin(sucsccs), cend(sucsccs));
}


template <typename T>
void mark_dead_sccs(SCCDat<T> const& scci, BaSccClassification const& scccl,
    function<vector<unsigned>(unsigned)> const& get_suc_sccs,
    map<unsigned,bool>& dead, unsigned num) {
  if (map_has_key(dead, num))  // done already
    return;

  // std::cout << "scc " << num << std::endl;
  auto sucsccs = get_suc_sccs(num);
  // mark children first
  for (auto sucnum : sucsccs) mark_dead_sccs(scci, scccl, get_suc_sccs, dead, sucnum);

  // if we are rejecting and trivial, assume we're dead
  bool isdead = contains(scccl.rejecting, num) || contains(scccl.trivial, num);
  // check children and try to falsify
  for (auto sucscc : sucsccs) isdead = isdead && dead[sucscc];
  // if still dead, we're really dead.
  dead[num] = isdead;
}

// run dfs that marks dead sccs
template <typename T>
set<unsigned> get_dead_sccs(SCCDat<T> const& scci, BaSccClassification const& scccl,
    function<vector<unsigned>(unsigned)> const& get_succs) {
  map<unsigned, bool> dead;
  for (int i=0; i<(int)scci.sccs.size() ; i++)
    mark_dead_sccs(scci, scccl, get_succs, dead, i);

  set<unsigned> ret;
  for (auto const& it : dead)
    if (it.second)
      ret.emplace(it.first);
  return ret;
}

//this is for trying stuff out. do whatever you want here.
int main(int argc, char *argv[]) {
  string file="";
  if (argc>1)
    file = argv[1];
  auto bas = parse_hoa_ba(file);
  if (bas.empty())
    cout << "something went wrong!" << endl;

  auto &aut = bas.front();

  cout << "starting get_sccs" << endl;

  function<vector<state_t>(state_t)> const sucs = [&aut](state_t v){ return aut->succ(v); };
  function<bool(state_t)> const ac = [&aut](state_t v){ return aut->has_accs(v); };

  auto const scci = get_sccs(aut->states(), sucs);
  auto const bascl = classify_sccs(scci, sucs, ac);
  auto const sucsccs = [&](unsigned num){ return succ_sccs(scci, num, sucs); };
  auto const badead = get_dead_sccs(scci, bascl, sucsccs);

  for (int i=0; i<(int)scci.sccs.size(); i++) {
    cout << seq_to_str(scci.sccs.at(i));
    if (contains(bascl.accepting, i))
      cout << " (acc)";
    if (contains(bascl.rejecting, i))
      cout << " (rej)";
    if (contains(bascl.trivial, i))
      cout << " (triv)";
    if (contains(badead, i))
      cout << " (dead)";
    cout << endl;
  }

  /*
  if (argc>2){
    auto ps = powerset_product(*ba);
    print_hoa(*ps);
  } else {
    auto ps = powerset_construction(*ba);
    print_hoa(*ps);
  }
  */

}
