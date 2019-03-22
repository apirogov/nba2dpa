#pragma once

#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <sstream>
// #include <string>

#include "common/types.hh"

#include <process.h>

std::vector<int> solve_maxsat(std::string formula) {
  using namespace std;

  //TODO: remove this wrapper, implement functionality here
  procxx::process maxsat{"pipemaxsat.rb"};
  maxsat.exec();
  maxsat << formula << endl;
  maxsat.close(procxx::pipe_t::write_end());

  string linestr;
  while (std::getline(maxsat.output(), linestr)) {
    istringstream line(linestr);
    string tok;

    line >> tok;
    if (tok == "c")
      continue;
    if (tok == "s")
      continue;
    if (tok == "o")
      continue;
    if (tok == "v") {
      //collect variable assignment
      vector<int> ret;
      int tmp;
      while (line >> tmp) {
        ret.push_back(tmp);
      }
      return ret;
    }
  }
  return {};
}

//input: list of constraints (each state+sym must have an edge from the set in the map)
//each value in sets on the right appears also on the left
//assuming the graph is a single SCC
//output: minimal nontrivial strongly connected subgraph satisfying constraints
//
//this is an NP complete problem, simple reduction from Set Cover
//(build bipartite graph with elements on left and sets on right, edges correspond to
//inclusion, the edges between sets from a single node on left = constraints,
//add dummy, connect right nodes to dummy + edges from dummy to left nodes (fixed edges))
std::vector<nbautils::state_t> altmap_to_maxsat(
    std::map<nbautils::state_t,std::map<nbautils::sym_t, std::vector<nbautils::state_t>>> const& altmap) {
  using namespace std;
  using namespace nbautils;

  std::stringstream ss;

  //build reverse mapping: predecessors of each state
  std::map<state_t, std::set<state_t>> preds;
  for (auto const& it : altmap) {
    state_t const cur = it.first;
    for (auto const& s2e : it.second)
      for (auto const suc : s2e.second)
        preds[suc].emplace(cur);
  }

  //we build a CNF for a partial MaxHornSAT instance

  //each state is mapped to a variable
  std::map<state_t, int> s2v;
  for (auto const& it : altmap)
    s2v[it.first] = s2v.size();
  //inverse mapping
  std::map<int, state_t> v2s;
  for (auto const& it : s2v)
    v2s[it.second] = it.first;

  //hard CNF clauses that may not be violated, contain variables x_i with i>0, -i = neg. literal
  std::vector<std::vector<int>> hclauses;

  //hard clause: at least one state is preserved ("not all states killed")
  std::vector<int> atleastone;
  for (auto const& it : s2v)
    atleastone.push_back(-it.second);
  hclauses.push_back(atleastone);

  //hard clause: for all q,a, if all successors on q,a dead -> q dead
  for (auto const& it : altmap) {
    state_t const cur = it.first;
    for (auto const& s2e : it.second) {
      if (s2e.second.empty()) //no successor for that letter -> add no constraint
        continue;

      std::vector<int> tmp;

      for (auto const suc : s2e.second)
        tmp.push_back(-s2v[suc]);
      tmp.push_back(s2v[cur]);

      hclauses.push_back(tmp);
    }
  }

  //hard clause: if all predecessors dead -> cur dead (useless state, not part of bottom SCC)
  for (auto const& it : preds) {
    state_t const cur = it.first;
    std::vector<int> tmp;

    for (auto const& p : it.second)
      tmp.push_back(-s2v[p]);
    tmp.push_back(s2v[cur]);

    hclauses.push_back(tmp);
  }

  //we actually optimize number of set variables encoded as partial max(Horn)SAT problem
  //hard constraints have topval weight, topval > sum of weights of soft constr.
  //hard constraints: consistency stuff for graph
  //soft constraints: trivial singleton variables (we want to set (=kill) as many as possible)
  //we set soft weights to 1, i.e. sum of states < topval
  int numvars = s2v.size();
  int topval = numvars + 1;
  ss << "p wcnf " << numvars << " " << (numvars+hclauses.size()) << " " << topval << endl;
  for (auto const& clause : hclauses) {
    ss << topval << " " << seq_to_str(clause, " ") << " 0" << endl;
  }
  for (auto const& it : s2v)
    ss << "1 " << it.second << " 0" << endl;

  // cerr << "CNF:" << endl << ss.str() << endl;
  auto const sol = solve_maxsat(ss.str());

  //start with all states
  set<state_t> keep;
  for (auto const& it : s2v)
    keep.emplace(it.first);

  for (int const var : sol) {
    if (var <= 0) //negative = kept
      continue;
    //positively set variable = killed state -> remove from kept
    state_t const st = v2s.at(var);
    keep.erase(keep.find(st));
  }

  vector<state_t> ret(begin(keep),end(keep));

  //sort by "usefulness" (how many others can point to that state)
  //this slightly helps the Hopcroft minimization
  sort(begin(ret), end(ret), [&](auto a, auto b){
      vector<nbautils::state_t> tmpa;
      set_intersection(begin(preds.at(a)),end(preds.at(a)),
                       begin(keep),end(keep),back_inserter(tmpa));
      vector<nbautils::state_t> tmpb;
      set_intersection(begin(preds.at(b)),end(preds.at(b)),
                       begin(keep),end(keep),back_inserter(tmpb));
      auto const sza = tmpa.size();
      auto const szb = tmpb.size();
      //sort desc. by # of touched sets, then asc. by val
      return  sza > szb || (sza == szb && a < b);
  });

  //return states to keep, in order of "usefulness"
  return ret;
}


