#include <iostream>
#include "io.hh"
#include "common/algo.hh"
#include "common/part_refinement.hh"
#include "pa.hh"

#include <range/v3/all.hpp>

using namespace std;
using namespace nbautils;

//this is for trying stuff out. do whatever you want here.
int main(int argc, char *argv[]) {
  string file="";
  if (argc>1)
    file = argv[1];
  auto bas = parse_hoa_ba(file);
  if (bas.empty())
    cout << "something went wrong!" << endl;

  auto &aut = bas.front();

  // function<vector<state_t>(state_t)> sucs(&aut->succ); //[&aut](state_t v){ return aut->succ(v); };
  /*
  function<vector<state_t>(state_t,sym_t)> const xsucs = [&aut](state_t v,sym_t s){ return aut->succ(v,s); };
  function<vector<sym_t>(state_t)> const outsyms = [&aut](state_t v){ return aut->outsyms(v); };
  function<bool(state_t)> const ac = [&aut](state_t v){ return aut->has_accs(v); };
  auto const states = aut->states();

  auto const unreach = unreachable_states(states, aut->get_init().front(), sucs);
  cout << "unreachable from initial: " << seq_to_str(unreach) << endl;

  auto const accsinks = ba_get_accepting_sinks(states, aut->num_syms(),
      ac, outsyms, xsucs);
  cout << "accepting sinks: " << seq_to_str(accsinks) << endl;

  auto const scci = get_sccs(aut->states(), sucs, const_true);
  auto const trivial = trivial_sccs(scci, sucs);
  auto const bascl = ba_classify_sccs(scci, ac);
  auto const sucsccs = [&](unsigned num){ return succ_sccs(scci, num, sucs); };
  auto const badead = ba_get_dead_sccs(scci.sccs.size(), bascl.rejecting, trivial, sucsccs);

  for (int i=0; i<(int)scci.sccs.size(); i++) {
    cout << seq_to_str(scci.sccs.at(i));
    if (contains(bascl.accepting, i))
      cout << " (acc)";
    if (contains(bascl.rejecting, i))
      cout << " (rej)";
    if (contains(trivial, i))
      cout << " (triv)";
    if (contains(badead, i))
      cout << " (dead)";
    cout << endl;
  }
  */

  /*
  SWA<Acceptance::PARITY, string> aut;
  aut.set_patype(PAType::MIN_EVEN);
  aut.add_state(0);
  aut.set_accs(0,{4});
  aut.add_state(1);
  aut.set_accs(1,{2});
  aut.add_state(2);
  aut.set_accs(2,{3});
  aut.add_state(3);
  aut.set_accs(3,{3});
  aut.add_state(4);
  aut.set_accs(4,{1});
  aut.add_state(5);
  aut.set_accs(5,{0});
  aut.set_aps({"x"});
  aut.set_succs(0,0,{1});
  aut.set_succs(1,0,{2});
  aut.set_succs(2,0,{0});
  aut.set_succs(2,1,{3});
  aut.set_succs(3,0,{4});
  aut.set_succs(4,0,{3});
  aut.set_succs(4,1,{5});
  aut.set_succs(5,0,{3});

  // function<vector<state_t>(state_t)> const sucs = [&](state_t v){ return aut.succ(v); };
  // function<unsigned(state_t)> const get_pri = [&](state_t v){return aut.get_accs(v).front();};
  auto const newpris = minimize_priorities(aut); //aut.states(), sucs, get_pri);
  for (auto p : aut.get_accsets()) {
    cout << p << " -> " << newpris(p) << endl;
  }
  */

  auto const accsinks = get_accepting_sinks(*aut);
  auto ps = powerset_construction(*aut, accsinks);
  print_hoa(*ps);

}
