/**
*	\file compl_constr2.hh
*	\author Lasse Nitz
*
*	This file contains methods that are used for the Construction 2, that creates a complementary
*	Büchi automaton, constisting of a powerset-automaton (PS) and rankings.
*	The resulting automaton is the complementary automaton @f$ \overline{\mathcal{A}} @f$, as defined
*	in the thesis "Optimization for the Complementation of Büchi Automata".
*	The top-level method that constructs the automaton @f$ \overline{\mathcal{A}} @f$ is compl_construction().
*/

#pragma once

// C++ Standard Library
#include <bitset>
#include <vector>

// nbautils
#include "aut.hh"
#include "common/scc.hh"

// Compl
#include "compl/compl_constr1.hh"
#include "compl/compl_print.hh"
#include "compl/compl_tag.hh"



namespace cmpl{

using namespace std;
using namespace nbautils;
using namespace cmpl;

using ComplAut = Aut<ComplTag>;




/////////// COMPL-POWERSET AUTOMATON ///////////////////////////////////////////////////////////
/**
*	\brief	Creates a powerset-automaton that uses ComplTags.
*
*	This is a mofified version for the powerset construction found in ps.hh.
*
*	\param nba	The automaton, for which the powerset-automaton should be constructed.
*	\param mat	The adjacency-matrix of nba.
*	\param sinks	Defines sinks of the input-automaton and does not define any sinks if the parameter is undefined. This parameter is only used as a parameter for the function powersucc() in aut.hh.
*
*	\return	The powerset-ComplTag-automaton for nba.
*/
ComplAut compl_ps_construction(auto const& nba, adj_mat const& mat, nba_bitset const& sinks=0) {
  assert(nba.is_buchi());

  // Create aut, add initial state and associate with initial states in original automaton
  state_t const myinit = 0;
  auto ps = ComplAut(true, nba.get_name(), nba.get_aps(), myinit);
  ps.tag_to_str = [](ostream& out, ComplTag const& t){
    print_compl_tag(out, t);
  };

  ComplTag ct1;					//'p', nba_bitset(1<<nba.get_init()), nba_bitset(1<<nba.get_init()));
  ct1.stateType = 'p';			// This will be a state in stage 1
  ct1.stateSet = nba_bitset(1<<nba.get_init());
  ps.tag.put(ct1, myinit);


  bfs(myinit, [&](auto const& st, auto const& visit, auto const&) {
    // Get inner states of current ps state
    auto const curset = ps.tag.geti(st).stateSet;
    // Calculate successors and add to graph
    for (auto const i : ps.syms()) {
      auto const sucset = powersucc(mat, curset, i, sinks);

		bool tagKnown = false;			// Is used to determine whether a tag with a certain stateSet exists
		unsigned int tagFoundAt = 0;	// If a tag with a certain stateSet exists, this variable stores the according state

		// Since tag stores pointers as values, this ugly approach is used instead of put_or_get()
		for(auto j : ps.states()){
			if(ps.tag.geti(j).stateSet == sucset){
				tagKnown = true;
				tagFoundAt = j;
				break;
			}
		}


		// If the stateSet sucset is unknown, add a state with sucset as its stateSet
		if (!tagKnown){
			// Add a new state
			ps.add_state(ps.num_states());

			// Create and set a tag for the new state
			ComplTag ct2;
			ct2.stateType = 'p';		// This will be a state in stage 1
			ct2.stateSet = sucset;
			ps.tag.put(ct2, ps.num_states()-1);
		}


		// Add edge & schedule bfs visit of successor
		if(tagKnown){			// Destination-state of transition already existed
          if (!ps.has_edge(st, i, tagFoundAt)) {
			ps.add_edge(st, i, tagFoundAt);
			visit(tagFoundAt);
          }
		}
		else{					// Destination-state of transition has just been added
			ps.add_edge(st,i,ps.num_states()-1);
			visit(ps.num_states()-1);
		}

    }
  });

  return ps;
}









/**
*	\brief	Copies a powerset-ComplTag-automaton and adds type-2 transitions from powerset-states to ranking-states.
*
*	\param ps	The powerset-ComplTag-automaton of nba.
*	\param nba	The nba, for which the powerset-automaton was created. In the context of Buechi-complementation, this is the input-automaton.
*
*	\return	A ComplTag-automaton, that consists of the powerset-automaton and type-2 transitions, including the according ranking states.
*/
ComplAut add_type2_trans_ps(ComplAut ps, auto const& nba){

	// Create STS to find out which slices are reachable
	auto sts = sts_construction(nba, get_adjmat(nba));

	state_t stateRange = ps.num_states();

	// Add type 2 transitions and according states to ps
	for(state_t p = 0; p < stateRange; p++){		// For each ps-state

		for(auto const x : ps.syms()){				// For each transition-type

			for(auto& succ : ps.succ(p, x)){
				//Determine successor in sts and add ranking for it to ps
				for(auto& s : sts.states()){//(auto& suc : sts.succ(s, x)){
					if(ps.tag.geti(succ).stateSet == sts.tag.geti(s).stateSet){
						vector<signed long> ranking = slice_to_rank(nba, sts.tag.geti(s).slice);
						bool tagKnown = false;
						unsigned long tagFoundAt;

						// Check if the ranking is already known (i.e., if the according state already exists)
						for(auto& j : ps.states()){
							if(ps.tag.geti(j).stateType == 'r' && ps.tag.geti(j).ranking == ranking){
								tagKnown = true;
								tagFoundAt = j;
								break;
							}
						}

						if(tagKnown){		// Ranking is already known, only add edge
                          if (!ps.has_edge(p,x,tagFoundAt))
							ps.add_edge(p, x, tagFoundAt);
						}
						else{				// Ranking is new, add state, tag and edge to new state
							ComplTag tag;
							tag.stateType = 'r';

							tag.stateSet = ps.tag.geti(succ).stateSet;
							tag.ranking = ranking;
							tag.obligationSet = 0;

							ps.add_state(ps.num_states());
							ps.add_edge(p,x,ps.num_states()-1);
							ps.tag.put(tag, ps.num_states()-1);
						}
					}
				}
			}
		}
	}

	return ps;
}





/////////// COMPLEMENTARY AUTOMATON A_comp ///////////////////////////////////////////////////////////
/**
*	\brief	Constructs and returns the complementary automaton @f$ \overline{\mathcal{A}} @f$ via Construction 2.
*
*	\param nba	Input-NBW that should be complemented.
*	\param mat	Adjacency-Matrix of the input-NBW
*
*	\return	The complementary NBW @f$ \overline{\mathcal{A}} @f$ for the input-NBW nba, constructed via Construction 2.
*/
ComplAut compl_construction(auto const& nba, adj_mat const& mat){
	assert(nba.is_buchi());

	// Create the Powerset-automaton for the given input-automaton nba
	auto PS = compl_ps_construction(nba, mat);

	// Create automaton from PS that contains transitions of type 2 (stateset to ranking) and the according ranking-states
	auto temp = add_type2_trans_ps(PS, nba);

	// Add transitions between rankings and the according states
	auto Acomp = add_type3_trans(temp, nba);

	// Add acceptance-property to states (final or non-final)
	add_acceptance(&Acomp);

	// Set name
	string name = Acomp.get_name();
	name += "_compl_2";
	Acomp.set_name(name);

	return Acomp;
}


}	// End of namespace cmpl
