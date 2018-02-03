#include "nba.hh"
#include <algorithm>
#include <memory>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;

#include <spot/parseaut/public.hh>
#include <spot/twaalgos/sccinfo.hh>

string const err_notsba = "It does not look like this is really a state-based NBA!";

namespace nbautils {

// load state-based NBA from file in a format supported by spot
spot::twa_graph_ptr spot_nba_from_file(string const &file) {
  spot::parsed_aut_ptr pa = parse_aut(file, spot::make_bdd_dict());

  if (pa->format_errors(std::cerr))  // prints parser diagnostic errors
    return nullptr;

  if (pa->aborted) {
    cerr << "--ABORT-- read" << endl;
    return nullptr;
  }
  if (!pa->errors.empty() || pa->aut == nullptr) {  // had serious errors
    cerr << "Serious errors occured reading " << pa->filename
         << ". Please check your input!" << endl;
    return nullptr;
  }

  if (!pa->aut->is_sba().is_true()) {
    cerr << err_notsba << endl;
    return nullptr;
  }

  return pa->aut;
}

// calculate bdds for each combination of atomic proposition
// useful for quicker iterated membership query of letters in edges
vector<bdd> calc_bdd_syms(spot::const_twa_ptr aut) {
  auto const ap = aut->ap();  // vector of atomic propositions of automaton
  auto const dict =
      aut->get_dict();  // map from atomic propositions to corresponding bdd variables

  auto const numsyms = 1 << ap.size();
  auto ret = vector<bdd>();
  // calculate a bdd for each conjunction of all variable literals (positive or negated)
  // each combination represents an individual actual symbol read by the automaton
  for (auto sym = 0; sym < numsyms; sym++) {
    auto tmp = sym;
    auto bsym = bddtrue;
    for (auto p : ap) {
      auto bvar = dict->var_map[p];
      bsym &= (tmp & 1 ? bdd_ithvar(bvar) : !bdd_ithvar(bvar));
      tmp >>= 1;
    }
    ret.push_back(bsym);
  }
  return ret;
}

// get a graph from a spot graph, eliminating bdd stuff
NBA::NBA(spot::twa_graph_ptr aut) {
  if (!aut->is_sba().is_true()) {
    throw runtime_error(err_notsba);
  }

  auto allaps = aut->ap();  // backup original aps
  aut->remove_unused_ap();  // clean up aps not used in labels of automaton
  num_syms =
      1 << (aut->ap().size());  // so many different automaton symbols we really need
  auto uaps = aut->ap();        // grab reduced aps

  // resulting order of aps is used in order like in the twa_graph, then the unused
  sort(allaps.begin(), allaps.end());
  sort(uaps.begin(), uaps.end());
  aps = aut->ap();  // grab the effective ones
  set_difference(allaps.begin(), allaps.end(), uaps.begin(), uaps.end(),
                 std::back_inserter(aps));  // add useless to back

  // removes unreachable or who cannot be part of inf. run
  aut->purge_dead_states();

  // get scc info on cleaned automaton
  scci = make_shared<spot::scc_info>(spot::scc_info(aut));

  // calculate bdds for edge queries
  auto const syms = calc_bdd_syms(aut);
  // get initial state number
  auto const initst = aut->get_init_state();
  init = aut->state_number(initst);
  // and then use calculated bdds to discover all "real" edges in given graph
  queue<const spot::state *> bfsq;
  bfsq.push(initst);
  while (!bfsq.empty()) {
    auto const curr = bfsq.front();
    auto const currnum = aut->state_number(curr);
    bfsq.pop();

    if (states.find(currnum) != states.end()) continue;  // have visited this one

    // mark as visited by putting into state set
    states.emplace(currnum);
    // mark as accepting, if appropriate
    if (aut->state_is_accepting(curr)) acc.emplace(currnum);

    // first add all successors to queue
    for (auto const s : aut->succ(curr)) bfsq.push(s->dst());

    // next fill out all successors of state separated by letter
    for (auto i = 0; i < num_syms; i++) {
      set<state_t> succs;  // set because could be added multiple times
      for (auto const s : aut->succ(curr)) {
        if (syms[i] == (s->cond() & syms[i])) {
          succs.emplace(aut->state_number(s->dst()));
        }
      }
      // add all found successors for letter
      std::copy(succs.begin(), succs.end(), std::back_inserter(adj[currnum][i]));
    }
  }
}

std::vector<NBA::state_t> NBA::powersucc(std::vector<state_t> const &ps, sym_t const &x) {
  set<state_t> suc;
  for (auto const p : ps) {
    auto const qs = adj[p][x];
    std::copy(qs.begin(), qs.end(), std::inserter(suc, suc.end()));
  }
  return vector<state_t>(suc.begin(), suc.end());
}

}  // namespace nbautils
