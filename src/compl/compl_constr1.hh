/**
*	\file compl_constr1.hh
*	\author Lasse Nitz
*
*	This file contains methods that are used for the Construction 1, that creates a complementary
*	Buechi automaton, constisting of a Slice-Transition-System (STS) and rankings.
*	The resulting automaton is the complementary automaton @f$ A_L @f$, as defined on page 17 of
*	"Unifying Buechi Complementation Constructions" (https://arxiv.org/pdf/1302.2675.pdf).
*	This implementation is limited to a single initial state.
*	The top-level method that constructs the automaton @f$ A_L @f$ is al_construction().
*/

#pragma once

// C++ Standard Library
#include <bitset>
#include <vector>

// nbautils
#include "aut.hh"

// Compl
#include "compl/compl_print.hh"
#include "compl/compl_tag.hh"






namespace cmpl{

using namespace std;
using namespace nbautils;
using namespace cmpl;

using ComplAut = Aut<ComplTag>;



/**
*	\brief	Returns the largest slice-value for x-predecessors of a given state.
*
*	This method returns the largest value in a slice-vector that is assigned to the x-predecessor in nba of the given state.
*
*	\param nba	The automaton, whose states define the elements of &preSlice and in which the x-predecessors are defined. In the context of Buechi-complementation, this is the input-automaton.
*	\param state	The state, for which the largest x-predecessor should be determined.
*	\param x		The letter, via which the predecessors are connected to focused state in the automaton nba.
*	\param preSlice	Pointer to the slice-vector that contains the values, from which the largest slice-value for the x-predecessors of the given state is chosen.
*
*	\return		The largest value in &preSlice for an x-predecessor of the given state.
*/
//TODO: this should be replaced with precomputed predecessors
signed long get_largest_predecessor(auto const& nba, state_t state, sym_t x, vector<signed long>* preSlice){
	signed long largest = -1;

	for(auto i : nba.states()){

		// Go through all successors of i for symbol x
		for(auto& scope : nba.succ(i, x)){

			if(scope == state && preSlice->operator[](i) != -1 && largest < preSlice->operator[](i)){
				largest = preSlice->operator[](i);
			}
		}
	}
	return largest;
}





// Currently not needed, only for mirrored version of the slice-reduction!
signed long get_lowest_predecessor(auto const& nba, state_t state, sym_t x, vector<signed long>* preSlice){
	signed long lowest = -1;

	for(auto i : nba.states()){

		// Go through all successors of i for symbol x
		for(auto& scope : nba.succ(i, x)){

			if(scope == state && preSlice->operator[](i) != -1 && (lowest == -1 || lowest > preSlice->operator[](i))){
				lowest = preSlice->operator[](i);
			}
		}
	}
	return lowest;
}


/**
*	\brief	Returns a normalized a slice-vector.
*
*	Normalized a slice-vector such that every value between 0 and the highest value
*	does appear in the vector.
*	Since slices define pre-orders, different vectors representing the same pre-order
*	result in the same normalized slice-vector.
*	This normalization allows to easily check whether two given slice-vectors represent the same pre-order:
*	If they do, then their normalized slice-vectors are identical.
*
*	\param slice	A pointer to the slice-vector that should be normalized.
*
*	\return	The normalized slice-vector.
*/
vector<signed long> normalize_slice(vector<signed long>* slice){

	vector<signed long> resSlice;
	long offset = 0;

	// Store the offset for each value at its position in a vector (will be substracted later)
	vector<signed long> offsetVector;
	for(auto i = 0; i <= *max_element(slice->begin(), slice->end()); i++){
		if(find(slice->begin(), slice->end(), i) == slice->end()){
			offset++;
		}
		offsetVector.push_back(offset);
	}

	// "Tighten" the slice to normalize it unique by substracting offset
	for(state_t i = 0; i < slice->size(); i++){
		if(slice->operator[](i) == -1){
			resSlice.push_back(-1);
		}
		else{
			resSlice.push_back(slice->operator[](i) - offsetVector[slice->operator[](i)]);
		}
	}

	return resSlice;
}




/////////// STS CONSTRUCTION /////////////////////////////////////////////////////////////////////////
/**
*	\brief	Constructs a Slice-Transition-System (STS).
*
*	\param nba		The automaton, for which the STS should be constructed. In the context of Buechi-complementation, this is the input-automaton.
*	\param mat		Adjacency-matrix of nba.
*	\param sinks	Defines sinks of the input-automaton and does not define any sinks if the parameter is undefined. This parameter is only used as a parameter for the function powersucc() in aut.hh.
*
*	\return The STS as a ComplTag-automaton.
*/
ComplAut sts_construction(auto const& nba, adj_mat const& mat, nba_bitset const& sinks=0){
	assert(nba.is_buchi());

	// Create automaton, add initial state, and associate with initial states in original aut
	state_t const myinit = 0;
	auto sts = ComplAut(true, nba.get_name(), nba.get_aps(), myinit);
	sts.tag_to_str = [](ostream& out, ComplTag const& t){
		print_compl_tag(out, t);
	};

	// Create tag for initial state
	ComplTag initTag;
	initTag.stateType = 's';
	initTag.stateSet = nba_bitset(1<<nba.get_init());
	// Fill slice of initial state
	for(auto i : nba.states()){
		if(i == nba.get_init())	{initTag.slice.push_back(0);}		// Only initial state is in the initial slice
		else					{initTag.slice.push_back(-1);}		// Give dummy value to all non-slice states
	}
	// Add initial tag
	sts.tag.put(initTag, 0);


   bfs(myinit, [&](auto const& st, auto const& visit, auto const&) {
    // get inner states of current sts state
    auto const curset = sts.tag.geti(st).stateSet;
	vector<signed long> curSlice = sts.tag.geti(st).slice;
    // calculate successors and add to graph
    for (auto const i : sts.syms()) {		// Go through the letters in the alphabet
		auto const sucset = powersucc(mat, curset, i, sinks);

		// Calculate successor-slice
		vector<signed long> tempSlice;

		// Known: curset, curSlice, sucset
		// Now compute sucSlice!

		for(auto x : nba.states()){
			if(sucset[x] == 0){			// State not part of the slice
				tempSlice.push_back(-1);
			}
			else{						// State is part of the slice, determine its position in the preorder
				auto largest = get_largest_predecessor(nba, x, i, &curSlice);

				assert(largest != -1);

				if(nba.get_pri(x) == 0){	// Final state
					tempSlice.push_back(2*largest + 1);
				}
				else{						// Non-final state
					tempSlice.push_back(2*largest);
				}
			}
		}

		// "Normalize" the slice to make it unique; necessary for testing if a certain slice already exists
		auto sucSlice = normalize_slice(&tempSlice);



        //TODO: this should be done using tag.put_or_get(...)
		// Check if sucSlice does already exist
		bool tagKnown = false;			// Is used to determine whether a tag with a certain stateSet exists
		unsigned int tagFoundAt = 0;	// If a tag with a certain stateSet exists, this variable stores the according state

		for(auto j : sts.states()){
			if(sts.tag.geti(j).slice == sucSlice){
				tagKnown = true;
				tagFoundAt = j;
				break;
			}
		}


	  // If the stateSet sucset is unknown, add a state with sucset as its stateSet
      if (!tagKnown){
		// Add a new state
        sts.add_state(sts.num_states());

		// Create and set a tag for the new state
		ComplTag ct2;
		ct2.stateType = 's';		// This will be a state in stage 1
		ct2.stateSet = sucset;
		ct2.slice = sucSlice;
		sts.tag.put(ct2, sts.num_states()-1);
	  }


      // Add edge & schedule bfs visit of successor
	  if(tagKnown){			// Destination-state of transition already existed
        if (!sts.has_edge(st, i, tagFoundAt)) {
          sts.add_edge(st, i, tagFoundAt);
          visit(tagFoundAt);
        }
	  }
	  else{					// Destination-state of transition has just been added
		sts.add_edge(st,i,sts.num_states()-1);
		visit(sts.num_states()-1);
	  }

    }
  });

  return sts;
}





/**
*	\brief	Converts a slice to a ranking.
*
*	Convert a given slice to a ranking via "torank"-function from
*	the paper "Unifying Buechi Complementation Constructions" (see page 16).
*	The input-automaton nba is necessary to be able to distinguish between final and non-final states.
*
*	\param nba	The automaton, whose states are used as elements of the given slice. In the context of Buechi-complementation, this is the input-automaton.
*	\param slice	The slice-vector that should be transformed to a ranking.
*
*	\return	A ranking-vector for the given slice-vector.
*/
vector<signed long> slice_to_rank(auto const& nba, vector<signed long> slice){

	vector<signed long> rank;		// Empty rank-vector

	// Fill rank-vector with the according ranking-values for each state.
	for(unsigned long i = 0; i < slice.size(); i++){

		signed long beta = 0;		// The beta-value for state i, initially 0

		// Test if the current state is actually in slice
		if(slice[i] == -1){
			rank.push_back(-1);
			continue;
		}

		// Find the actual beta-value for the current state i:
		// For each position greater than the one of i, look if there is a non-final state at this position.
		// If so, increase the beta-value of i by 1.
		for(signed long x = *max_element(slice.begin(), slice.end()); x > slice[i]; x--){

			for(unsigned long j = 0; j < slice.size(); j++){

				if(slice[j] == x && nba.get_pri(j) == 1){
					beta++;
					break;
				}

			}
		}

		// Current state is part of slice, now give it its rank
		if(nba.get_pri(i) == 0){	// Final state
			rank.push_back(2*beta);
		}
		else{						// Non-final state
			rank.push_back(2*beta+1);
		}

	}

	return rank;
}















/**
*	\brief	For a given ComplTag-automaton, a ComplTag-automaton with type-2 transitions and according states is created and returned.
*
*	Adds transitions from slice-states to ranking-states, as defined on pages 16 and 17 of "Unifying Buechi Complementation Constructions".
*	If the receiving ranking-state does not exist already, it is created.
*
*	\param sts	A ComplTag-automaton with Slice-states. In the context of Buechi-complementation, this should be the Slice-Transition-System for nba (stage 1).
*	\param nba	In the context of Buechi-complementation, this should be the input-automaton.
*
*	\return
*/
ComplAut add_type2_trans_sts(ComplAut const& sts, auto const& nba){

	ComplAut res = sts;

	for(auto& i : sts.states()){

		for(auto const x : sts.syms()){		// Go through all elements in alphabet

			for(auto& suc : sts.succ(i, x)){		// Go through all x-successors (only one, since STS is deterministic)

				// Generate ranking for that slice
				vector<signed long> ranking = slice_to_rank(nba, sts.tag.geti(suc).slice);
				bool tagKnown = false;
				unsigned long tagFoundAt;

				// Test if this ranking does already exist
				for(auto& j : res.states()){
					if(res.tag.geti(j).stateType == 'r' && res.tag.geti(j).ranking == ranking){
						tagKnown = true;
						tagFoundAt = j;
						break;
					}
				}

				// Modify output-automaton
				if(tagKnown){			// Ranking is already known, only add edge
                  if (!res.has_edge(i,x,tagFoundAt))
					res.add_edge(i, x, tagFoundAt);
				}
				else{					// Ranking is new, add state, tag and edge to new state
					ComplTag tag;
					tag.stateType = 'r';

					tag.stateSet = res.tag.geti(suc).stateSet;
					tag.ranking = ranking;
					tag.obligationSet = 0;

					res.add_state(res.num_states());
					res.add_edge(i,x,res.num_states()-1);
					res.tag.put(tag, res.num_states()-1);
				}
			}
		}
	}
	return res;
}




/**
*	\brief	Tightens a ranking-vector.
*
*	Tightens a ranking according to the textual description of the tighten-function on page 17 of
*	the paper "Unifying Buechi Complementation Constructions".
*	For the formal definition, refer to "Optimization for the Complementation of Buechi Automata".
*	The nba as input-automaton is primarily used to determine which states are final.
*
*	\param nba	The automaton, for whose states the given ranking is defined.
*	\param ranking	The pointer to the ranking-vector, that should be tightened.
*
*	\return	A tight ranking-vector.
*/
vector<signed long> tighten_ranking(auto const& nba, vector<signed long>* ranking){

	vector<signed long> tightRanking;

	for(auto& i : nba.states()){

		// Case 1: If the ranking of state i is -1 in the input ranking, it also is -1 in the tight ranking
		if(ranking->operator[](i) == -1){
			tightRanking.push_back(-1);
			continue;
		}

		// Determine gamma-value for each state i of the input-automaton
		unsigned long gamma = 0;
		for(signed long j = 1; j < ranking->operator[](i); j += 2){
			if(find(ranking->begin(), ranking->end(), j) != ranking->end()){
				gamma++;
			}
		}

		// Add ranking value
		if(nba.get_pri(i) == 0 || ranking->operator[](i) % 2 == 0){	// Case 2: i is a final state
			// The 2nd case is not defined clearly in the source-paper, but solves the correctness problem
			tightRanking.push_back(2*gamma);
		}
		else{						// Case 3: i is a non-final state
			tightRanking.push_back(2*gamma + 1);
		}
	}

	return tightRanking;
}





/**
*	\brief	Generates the successor-ComplTag of a ranking-ComplTag, automaton and letter.
*
*	Since ranking-states constist of both ranking and obligation-set, this method returns a ComplTag, and not
*	just the representation of the ranking-function itselt.
*	Additionally to the ranking and obligationSet, the variable stateSet of the returned ComplTag defines,
*	for which states the ranking is defined, i.e. is not equal to -1.
*	For details of a successor-ranking and successor-obligation-sets, refer to page 17 of the paper "Unifying Buechi Complementation Constructions".
*
*	\param nba	The automaton, from which the states for the ranking originate. In the context of Buechi-complementation, this is the input-NBW that should be complemented.
*	\param curTag	The ComplTag that includes the ranking, for which the successor should be computed.
*	\param x	The letter that defines which successor-ranking should be determined for curTag.
*
*	\return	The x-successor ComplTag of the ranking-ComplTag curTag.
*/
ComplTag succ_ranking(auto const& nba, ComplTag curTag, sym_t x){
	assert(curTag.stateType == 'r');

	ComplTag tag;
	vector<signed long> succ_ranking;
	bitset<max_nba_states> succ_obSet;

	// Generate ranking-value for each state of the given nba
	for(auto& i : nba.states()){

		vector<state_t> pred;	// Vector of x-predecessors of state i with a ranking different from -1

		// Fill pred
		for(auto& j : nba.states()){
			// Is there a way to directly check if an element is in a given range?
			for(auto const p : nba.succ(j, x)){
				if(i == p && curTag.ranking[j] != -1){
					pred.push_back(j);
				}
			}
		}

		// Calculate ranking value for state i
		if(pred.size() == 0){			// Case 1 (pred is empty)
			succ_ranking.push_back(-1);
			continue;
		}
		else{							// Case 2 and 3 (pred is not empty)
			signed long smallestRanking = -1;
			// Find smallest ranking among predecessors
			for(unsigned long p = 0; p < pred.size(); p++){
				if(smallestRanking == -1 || smallestRanking > curTag.ranking[pred[p]]){
					smallestRanking = curTag.ranking[pred[p]];
				}
			}
			// Case 2: Substract 1, if i is final state and smallestRanking is odd
			if(nba.get_pri(i) == 0 && smallestRanking % 2 == 1){
				smallestRanking--;
			}
			succ_ranking.push_back(smallestRanking);
		}
	}

	// Tighten the successor-ranking
	auto tight_succ = tighten_ranking(nba, &succ_ranking);

	// Calculate obligation-set for the successor-ranking
	if(curTag.obligationSet != 0){		// Non-empty Obligation-Set (Case 1)
		// New obligation set is Delta(curTag.obligationSet, x) without odd(succ_ranking)
		succ_obSet = powersucc(get_adjmat(nba), curTag.obligationSet, x);
		for(auto& i : nba.states()){
			if(tight_succ[i] % 2 == 1 && tight_succ[i] != -1){
				succ_obSet.reset(i);
			}
		}
	}
	else{								// Empty Obligation-Set (Case 2)
		// Add all states with even rank to the new obligation set
		for(auto& i : nba.states()){
			if(tight_succ[i] % 2 == 0){
				succ_obSet.set(i);
			}
		}
	}

	// Fill output-tag with information
	tag.stateType = 'r';
	tag.stateSet = powersucc(get_adjmat(nba), curTag.stateSet, x);
	tag.ranking = tight_succ;
	tag.obligationSet = succ_obSet;

	return tag;
}




/**
*	\brief For a given ComplTag-automaton, a ComplTag-automaton with type-3 transitions and according states is created and returned.
*
*	Type-3 transitions (transitions between rankings) according to the a-successor definition of rankings are added, including
*	the possibly new states associated with these transitions.
*
*	\param comp	The ComplTag-automaton, on which the returned automaton is based. It should contain type-1 and type-2 transitions, if this method is called in the context of Buechi-complementation.
*	\param nba 	The input-automaton that should be complemented.
*
*	\return	A ComplTag automaton with type-3 transitions, i.e. all ranking-successors have been added for already existing ranking-states and ranking-states that were created in the process.
*/
ComplAut add_type3_trans(ComplAut const& comp, auto const& nba){

	ComplAut res = comp;

	bfs(res.get_init(), [&](auto const& st, auto const& visit, auto const&) {
		if(st == comp.get_init()){for(auto& i : comp.states()){visit(i);}}		// Add all existing states to visit

		if(res.tag.geti(st).stateType != 'r'){return;}		// Only ranking-states are of interest

		// Generate the x-successor of state st for each x in the alphabet
		for(auto const x : res.syms()){

			auto succ_tag = succ_ranking(nba, res.tag.geti(st), x);		// Get tag for x-successor of st

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
				res.add_edge(st, x, tagFoundAt);
				visit(tagFoundAt);
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





/**
*	\brief Adds acceptance (final/non-final) to all states in the ComplTag-automaton that is given via its pointer.
*
*	A state of the automaton is final if, and only if, the according ComplTag is of stateType 'r' and has an empty
*	obligationSet (i.e., no bit of the obligationSet is set to 1).
*	This implements the acceptance for Constructions 1, 2 and 3.
*
*	\param aut	Pointer to the ComplTag-automaton, for which acceptance should be added.
*/
void add_acceptance(Aut<ComplTag>* aut){
	for(auto& i : aut->states()){
		if(aut->tag.geti(i).stateType == 'r' && aut->tag.geti(i).obligationSet == 0){
			aut->set_pri(i, 0);		// Ranking-states with empty Obligation-Set are final
		}
		else{
			aut->set_pri(i, 1);		// All other states are non-final
		}
	}
}





/////////// COMPLEMENTARY AUTOMATON A_L //////////////////////////////////////////////////////////////
/**
*	\brief	Constructs and returns the complementary automaton @f$ A_L @f$ via Construction 1.
*
*	\param nba	Input-NBW that should be complemented.
*	\param mat	Adjacency-Matrix of the input-NBW
*
*	\return	The complementary NBW @f$ A_L @f$ for the input-NBW nba, constructed via Construction 1.
*/
ComplAut al_construction(auto const& nba, adj_mat const& mat){
	assert(nba.is_buchi());

	// Create the STS for the given input-automaton nba
	auto STS = sts_construction(nba, mat);

	// Create automaton from STS that contains transitions of type 2 (slice to ranking) and the according ranking-states
	auto temp = add_type2_trans_sts(STS, nba);

	// Add transitions between rankings and the according states
	auto AL = add_type3_trans(temp, nba);

	// Add acceptance-property to states (final or non-final)
	add_acceptance(&AL);

	// Set name
	string name = AL.get_name();
	name += "_compl_1";
	AL.set_name(name);

	return AL;
}


}	// End of namespace cmpl
