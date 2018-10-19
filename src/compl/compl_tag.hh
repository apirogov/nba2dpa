/**
*	\file compl_tag.hh
*	\author Lasse Nitz
*
*	This file contains the custom tag used for the complementation of Buechi-automata
*	and defines according methods like a comparator and a hash-function.
*/

#pragma once

// C++ Standard Library
#include <bitset>
#include <vector>

// nbautils
#include "aut.hh"

namespace cmpl{

using namespace std;
using namespace nbautils;
using namespace cmpl;

/**
*	\brief	This class is used to represent tags for the states in the complementary Buechi-automaton.
*/
class ComplTag{
public:
	// Constructor
	ComplTag(){};

	// Currently (17.09.2018): nbautils::max_nba_states = 256 in aut.hh
	
	// Data
	char stateType;						/// Defines the type of the state. p: state on powerset-stage, r: ranking, s: slice	
	bitset<max_nba_states> stateSet;	/// Represents a stateset of the input-automaton for complementation.
	vector<signed long> slice;			/// Vector that represents a slice. For more information, see the documentation of print_compl_slice in compl_print.hh
	vector<signed long> ranking;		/// Vector that represents a ranking. For more information, see the documentation of print_compl_ranking() in compl_print.hh
	bitset<max_nba_states> obligationSet; /// Represents a stateset of the input-automaton, used as obligation-set.
	
	// Comparator
	bool operator==(const ComplTag &toComp) const{
		return stateType == toComp.stateType
			&& stateSet == toComp.stateSet
			&& slice == toComp.slice
			&& ranking == toComp.ranking
			&& obligationSet == toComp.obligationSet;
	}
};

}	// End of namespace cmpl


namespace std{

	using namespace cmpl;	// To find ComplTag
	
/**
*	\brief Definition of a hash-function for objects of type ComplTag.
*
* 	Hash-function to hash ComplTag-objects such that they can be stored in 
* 	an unordered map (used for managing tags of automata).
* 	It is constructed similar to the hash functions found in pa.hh.
*/
	template <>
	struct hash<ComplTag> {
	  size_t operator()(ComplTag const& k) const {
		// Compute individual hash values for each sub-type
		// http://stackoverflow.com/a/1646913/126995
		size_t res = 17;
		res = res * 31 + hash<char>()( k.stateType );
		res = res * 31 + hash<bitset<max_nba_states>>()( k.stateSet );
		
		for(auto i : k.slice){
			res = res * 31 + (hash<signed long>() ( k.slice[i] ));
		}
		
		for(auto i : k.ranking){
			res = res * 31 + (hash<signed long>() ( k.ranking[i] ));
		}
		
		res = res * 31 + hash<bitset<max_nba_states>>()( k.obligationSet );

		return res;
	  }
	};
}	// End of namespace std