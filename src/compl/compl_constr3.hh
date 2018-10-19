/**
*	\file compl_constr3.hh
*	\author Lasse Nitz
*
*	This file contains methods that are used for the Construction 3, that creates an optimized complementary
*	Büchi automaton, constisting of a powerset-automaton (PS) and rankings.
*	The resulting automaton is the complementary automaton @f$ \overline{\mathcal{A}}_{SCC} @f$, as defined
*	in the thesis "Optimization for the Complementation of Büchi Automata".
*	The top-level method that constructs the automaton @f$ \overline{\mathcal{A}}_{SCC} @f$ is compl_construction_opt().
*/

#pragma once

// C++ Standard Library
#include <bitset>
#include <map>
#include <vector>

// nbautils
#include "aut.hh"
#include "common/scc.hh"

// Compl
#include "compl/compl_constr1.hh"
#include "compl/compl_constr2.hh"
#include "compl/compl_print.hh"
#include "compl/compl_tag.hh"


namespace cmpl{

using namespace std;
using namespace nbautils;
using namespace cmpl;

using ComplAut = Aut<ComplTag>;


/**
*	\brief Costum comparator to compare bitsets.
*
*	For larger values of max_nba_states (see aut.hh), the transformation from bitsets of size max_nba_states
*	to unsigned long long variables is not possible, due to a lack of bits in that type.
*	This comparator offers the possibility to use bitsets as keys in maps to avoid a limitation regarding the
*	size of input-automata by the unsigned long long type.
*
*/
struct BitsetComp{	// Custom comparator to compare bitsets
	bool operator() (const bitset<max_nba_states> &b1, const bitset<max_nba_states> &b2) const{
		// max_nba_states might be larger than the amount of bits for unsigned long long,
		// thus this approach has to be chosen
		for(signed long i = max_nba_states-1; i >= 0; i--){
			if(b1[i] && !b2[i]){ return false;}
			if(!b1[i] && b2[i]){ return true;}
		}
		return false;
	}
};



/**
*	\brief Determines unreachable states in a ComplTag-aut.
*
*	This method assumes that the initial state is 0.
*	The returned vector is always sorted ascendingly.
*
*	\param aut	The automaton for which the unreachable states should be found.
*
*	\return	An ascendingly sorted vector of all unreachable states in the given automaton aut.
*/
vector<state_t> get_unreachable_states(auto const& aut){	// Currently not used

	// Find reachable states
	map<state_t, bool> isReachable;
	isReachable[0] = true;
	for(auto p : isReachable){

		for(auto x : aut.syms()){
			for(auto s : aut.succ(p.first, x)){
				if(isReachable[s] != true){
					isReachable[s] = true;
				}
			}
		}
	}

	// Create vector with unreachable states
	vector<state_t> res;
	for(auto i : aut.states()){
		if(isReachable[i] != true){
			res.push_back(i);
		}
	}

	return res;
}




/**
*	\brief Determines recursively whether a certain SCC is a last SCC with a certain property.
*
*	\param sts	The automaton, of which the SCC-structure is taken.
*	\param vec	A vector that contains all SCCs from sts that fulfill a certain property. In the context of Buechi-complementation, this vector contains all SCCs that contain at least one state with a certain stateSet.
*	\param scc	The SCC, for which it should be determined whether another SCC in vec is reachable from it. This scc should appear in in the vector vec.
*	\param sccDat	SCC-data of sts.
*
*	\return	True, if no other SCC in vec is reachable from the given scc False otherwise.
*/
bool is_last_scc(ComplAut const& sts, vector<unsigned> vec, unsigned scc, SCCDat sccDat){
	auto succ = succ_sccs(aut_succ(sts), sccDat,scc);

	// Check direct successor-SCCs
	for(unsigned i = 0; i < vec.size(); i++){
		if(find(succ.begin(), succ.end(), vec[i]) != succ.end()){
			return false;
		}
	}

	// Recursively check further away SCCs
	bool res = true;
	for(auto const& i : succ){
		res = res && is_last_scc(sts, vec, i, sccDat);
	}

	return res;
}






/**
*	\brief	Returns a vector of "last SCCs" for a given stateSet.
*
*	\param sts		The ComplTag-automaton from which the SCC-structure is taken.
*	\param stateSet	The stateSet, for which the last SCCs (in which it appears) should be determined.
*
*	\return	A vector of all SCCs from the automaton sts such that from these SCCs no other SCC containing a state with the given stateSet is reachable.
*/
vector<unsigned> last_sccs_for_stateset(ComplAut const& sts, bitset<max_nba_states> stateSet){

	vector<unsigned> vec;		// Will constain all SCCs in which a slice with the given stateSet appears
	vector<unsigned> res;		// Will contain all "last SCCs" in which a slice with the given stateSet appears

	SCCDat sccDat = get_sccs(sts.states(), aut_succ(sts), true);

	// Fill vec
	for(auto st : sts.states()){
		unsigned temp = sccDat.scc_of[st];
		if(sts.tag.geti(st).stateSet == stateSet && find(vec.begin(), vec.end(), temp) == vec.end()){
			vec.push_back(temp);
		}
	}

	// Fill res
	for(unsigned scc = 0; scc < vec.size(); scc++){
		if(is_last_scc(sts, vec, vec[scc], sccDat)){
			res.push_back(vec[scc]);
		}
	}

	return res;
}




/**
*	\brief	Copies ps and adds adds type-2 transitions for construction 3.
*
*	This method adds transitions from powerset-states to ranking-states, with the
*	optimization presented in "Optimization for the Complementation of Buechi Automata".
*
*	\param ps	The powerset-ComplTag-automaton of nba.
*	\param nba	An automaton. In the context of Buechi-complementation, this is the input-automaton.
*	\param heuristic	Determines which heuristic should be used to choose a last SCC in the STS of nba for each state-set appearing in ps.
*					0: The last SCC with the smallest number assigned to it in the SCC-Data of the STS. This heuristic is the fastest.
*					1: The last SCC consisting of the least amount of states. If this last SCC is not unique, the one with the smallest number assigned to it by the SCC-Data is chosen.
*					2: The last SCC containing the least amount of states with the focused stateSet is chosen.
*					else:	Randomly chosen
*
*	\return	A ComplTag-automaton consisting of the powerset-automaton and type-2 transitions (via construction 3).
*/
ComplAut add_type2_trans_opt(ComplAut const& ps, auto const& nba, unsigned short heuristic){

	ComplAut res = ps;

	ComplAut sts = sts_construction(nba, get_adjmat(nba));
	SCCDat sccDat = get_sccs(sts.states(), aut_succ(sts), true);		// SCC-Data for sts

	map<bitset<max_nba_states>, vector<unsigned>, BitsetComp> lastSCCs;	// Key: stateSet, Value: All SCCs in the STS that are a last SCC for the stateSet
	map<bitset<max_nba_states>, unsigned, BitsetComp> chosenLastSCC;	// Key: stateSet, Value: One "last SCC" that is chosen from the vector of "last SCCs"

	// Fill lastSCCs
	for(auto st : sts.states()){

		auto b = sts.tag.geti(st).stateSet;

		if(lastSCCs.find(b) == lastSCCs.end()){	// Check if this stateSet has not been added as a key yet
			lastSCCs[b] = last_sccs_for_stateset(sts, b);
		}
	}

	// Fill chosenLastSCC by choosing one "last SCC" from the vector of each stateSet
	for(auto const& k : lastSCCs){

		// HEURISTICS //
		switch(heuristic){

		case 0:{	// Heuristic 0: Take the first SCC from the vector
			chosenLastSCC[k.first] = k.second[0];		// Take the first SCC from the vector
			break;}

		case 1:{	// Heuristic 1: Choose smallest SCC regarding the amount of states in the SCC
			auto sccVec = k.second;
			int smallestSCC = -1;
			for(auto i : sccVec){
				if(smallestSCC == -1 || sccDat.sccs[i].size() < sccDat.sccs[smallestSCC].size()){
					smallestSCC = i;
				}
			}
			chosenLastSCC[k.first] = smallestSCC;
			break;}

		case 2:{	// Heuristic 2: Choose SCC with the least appearances if the according stateSet
			auto sccVec = k.second;
			int smallestSCC = -1;
			int curStateSetCount = 0;
			int smallestStateSetCount = -1;
			for(auto i : sccVec){
				for(auto j : sccDat.sccs[i]){
					if(k.first == sts.tag.geti(j).stateSet){
						curStateSetCount++;
					}
				}
				if(smallestStateSetCount == -1 || curStateSetCount < smallestStateSetCount){
					smallestStateSetCount = curStateSetCount;
					smallestSCC = i;
				}
				curStateSetCount = 0;
			}
			chosenLastSCC[k.first] = smallestSCC;
			break;}

		default:{	// Randomly choose a last SCC
			auto randEntry = rand() % k.second.size();
			chosenLastSCC[k.first] = k.second[randEntry];
			break;}
		}
	}



	// Add type-2-transitions and according states
	for(auto st : ps.states()){		// Outgoing edges of all 'p'-type states (thus iterating over states of ps is ok)

		for(auto const& x : ps.syms()){		// Symbols are the same for ps, res and sts

			for(auto suc : ps.succ(st, x)){	// Contains only one element, since ps should be deterministic

				unsigned tempSCC = chosenLastSCC[ps.tag.geti(suc).stateSet];
				vector<unsigned> sucVec = sccDat.sccs[tempSCC];		// Contains all slice-states in the STS that are in the chosen "last SCC"

				for(auto const& s : sucVec){

					if(sts.tag.geti(s).stateSet != ps.tag.geti(suc).stateSet){continue;}	// Skip slices with a different stateSet than the x-successor in ps

					vector<signed long> ranking = slice_to_rank(nba, sts.tag.geti(s).slice);

					bool tagKnown = false;
					unsigned long tagFoundAt;

					// Check if the ranking is already known
					for(auto j : res.states()){
						if(res.tag.geti(j).stateType == 'r' && res.tag.geti(j).ranking == ranking){
							tagKnown = true;
							tagFoundAt = j;
							break;
						}
					}

					if(tagKnown){		// Tag is known, only add edge
                      if (!res.has_edge(st,x,tagFoundAt))
						res.add_edge(st, x, tagFoundAt);
					}
					else{				// Tag is unknown, add new state, edge and tag
						ComplTag tag;
						tag.stateType = 'r';
						tag.stateSet = ps.tag.geti(suc).stateSet;
						tag.ranking = ranking;
						tag.obligationSet = 0;

						res.add_state(res.num_states());
						res.add_edge(st, x, res.num_states()-1);
						res.tag.put(tag, res.num_states()-1);
					}

				}
			}
		}
	}


	return res;
}





/**
*	\brief	Copies comp and adds type-3 transitions according to construction 3 to it.
*
*	All successor-rankings of ranking-states in comp and newly added ranking-states are determined, with the
*	optimization defined in "Optimization for the Complementation of Buechi Automata".
*
*	\param comp	A ComplTag-automaton that consists of type-1 and type-2 transitions of construction 3.
*	\param nba	An automaton. In the context of Buechi-complementation, this is the input-automaton.
*
*	\return	A ComplTag automaton consisting of the automaton given via comp, with type-3 transitions and according state.
*/
ComplAut add_type3_trans_opt(ComplAut const& comp, auto const& nba){


	ComplAut res = comp;

	//// OPT ///////////////////////////////////////////////////////////////////////
	// Construct the powerset-automaton for nba, and get the according SCC-data
	auto ps = compl_ps_construction(nba, get_adjmat(nba));
	SCCDat psSCC = get_sccs(ps.states(), aut_succ(ps), true);

	// Key: a ComplTag-stateSet; Value: the SCC of this stateSet in the powerset-automaton
	map<bitset<max_nba_states>, state_t, BitsetComp> scc_of_stateset;

	// Fill scc_of_stateset
	for(auto p: ps.states()){
		scc_of_stateset[ps.tag.geti(p).stateSet] = psSCC.scc_of[p];
	}
	////////////////////////////////////////////////////////////////////////////////

	bfs(res.get_init(), [&](auto const& st, auto const& visit, auto const&) {
		if(st == comp.get_init()){for(auto& i : comp.states()){visit(i);}}		// Add all existing states to visit

		if(res.tag.geti(st).stateType != 'r'){return;}		// Only ranking-states are of interest

		// Generate the x-successor of state st for each x in the alphabet
		for(auto const x : res.syms()){

			auto succ_tag = succ_ranking(nba, res.tag.geti(st), x);		// Get tag for x-successor of st

			///// OPTIMIZATION /////////////////////////////////////////////////////////////////////
			if(scc_of_stateset[res.tag.geti(st).stateSet] != scc_of_stateset[succ_tag.stateSet]){continue;}	// OPT
			////////////////////////////////////////////////////////////////////////////////////////

			// Check if succ_tag does already exist
			bool tagKnown = false;
			unsigned long tagFoundAt;
			for(auto& i : res.states()){
				if(res.tag.geti(i).ranking == succ_tag.ranking && res.tag.geti(i).obligationSet == succ_tag.obligationSet){
					tagKnown = true;
					tagFoundAt = i;
					break;
				}
			}

			// Modify output-automaton
			if(tagKnown){			// Ranking is already known, only add edge
              if (!res.has_edge(st,x,tagFoundAt)) {
				res.add_edge(st, x, tagFoundAt);
				visit(tagFoundAt);
              }
			}
			else{					// Ranking is new, add state, tag and edge to new state
				res.add_state(res.num_states());
				res.add_edge(st, x, res.num_states()-1);
				res.tag.put(succ_tag, res.num_states()-1);
				visit(res.num_states()-1);
			}
		}

	});

	return res;
}







/////////// COMPLEMENTARY AUTOMATON A_comp_scc ///////////////////////////////////////////////////////
/**
*	\brief	Constructs and returns the complementary automaton @f$ \overline{\mathcal{A}}_{SCC} @f$ via Construction 3.
*
*	\param nba	Input-NBW that should be complemented.
*	\param mat	Adjacency-Matrix of the input-NBW
*
*	\return	The complementary NBW @f$ \overline{\mathcal{A}}_{SCC} @f$ for the input-NBW nba, constructed via Construction 3.
*/
ComplAut compl_construction_opt(auto const& nba, adj_mat const& mat){
	assert(nba.is_buchi());

	// Create the Powerset-automaton for the given input-automaton nba
	auto PS = compl_ps_construction(nba, mat);

	// Create automaton from PS that contains transitions of type 2 (stateset to ranking) and the according ranking-states
	auto temp = add_type2_trans_opt(PS, nba, 1);

	// Add transitions between rankings and the according states
	auto Acomp = add_type3_trans_opt(temp, nba);

	// Add acceptance-property to states (final or non-final)
	add_acceptance(&Acomp);

	// Set name
	string name = Acomp.get_name();
	name += "_compl_3";
	Acomp.set_name(name);

	return Acomp;
}


}	// End of namespace cmpl
