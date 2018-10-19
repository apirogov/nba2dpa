/**
*	\file compl_print.hh
*	\author Lasse Nitz
*
*	This file contains methods that print the custom ComplTag tag-type, as defined in compl_tag.hh.
*	Additionally, it defines some print-functions for data-types like SCCDat, which were used for
*	debugging the compl-module.
*/


#pragma once

// C++ Standard Library
#include <iomanip>
#include <iostream>
#include <bitset>

// Compl
#include "common/scc.hh"
#include "compl/compl_tag.hh"

namespace cmpl{

using namespace std;
using namespace nbautils;
using namespace cmpl;


/**
*	\brief	Method that prints bitsets of size max_nbs_states (see aut.hh) as statesets to a defined ostream.
*
*	In the context of the compl-module, this method is used to print the stateSet and obligationSet of a ComplTag.
*
*	\param bs	Pointer to the bitset of size max_nba_states that should be printed as a stateset.
*	\param out 	The ostream, to which the result should be output. If not defined, cout is used.
*/
void print_compl_bitset(const bitset<max_nba_states>* bs, ostream &out = cout){
	
	bool first = true;		// is true until the first positions with 1-bit has been found
	out << "{";
	for(size_t i = 0; i < bs->size(); i++){
		if(bs->operator[](i) == 1){
		
			if(!first){
				out << ", ";
			}
			
			out << "q" << i;
			first = false;
		}
	}
	out << "}";
}


/**
*	\brief Prints a vector<signed long> as a slice to a defined ostream.
*
*	In the context of the compl-module, this method is used to print the slice-variable of a ComplTag.
*	The value n saved at position x of the vector defines the equivalence class, in which state x is.
*	The larger n is, the larger the equivalence class is in the pre-order.
*	If n is -1, state x is not part of the pre-order. A value should under no circumstances be smaller than -1.
*
*	\param vec	Pointer to a vector of signed longs that should be printed as a slice.
*	\param out 	The ostream, to which the result should be output. If not defined, cout is used.
*/
void print_compl_slice(const vector<signed long>* vec, ostream &out = cout){
	
	bool first = true;			// Remembers if there is an element in the currently checked equivalence class
	bool firstEqClass = true;	// Remembers if there was an element that was not -1
	size_t count = 0;
	signed long curValue;
	
	if(*max_element(vec->begin(), vec->end()) == -1){	// Empty slice
		out << "[]";
		return;
	}
	
	for(auto j = 0; j <= *max_element(vec->begin(), vec->end()); j++){		// Go through eq.-classes, start with smallest one		
		
		for(size_t i = 0; i < vec->size(); i++){	// Find all states that belong to the currently focused eq.-class
			
			if(vec->operator[](i) == -1){continue;}
			
			curValue = vec->operator[](i);
			
			if(curValue == j){
				
				if(!firstEqClass && first){
					out << "] < ";
				}				
			
				if(!first){
					out << ", ";
				}
				else{
					out << "[";
				}
				out << "q" << i;
				count++;
				first = false;
				firstEqClass = false;
			}
		}
		
		first = true;
	}
	
	if(!firstEqClass){
		out << "]";
	}
}


/**
*	\brief Prints a vector<signed long> as a ranking to a defined ostream.
*
*	In the context of the compl-module, this method is used to print the ranking-variable of a ComplTag.
*	The value n stored at position x of the vector denotes the rank of state x.
*	If n is -1, state x is not part of the ranking, i.e. has the ranking "bottom" for well-defined rankings.
*	No value of the input-vector should be below -1.
*
*	\param vec	Pointer to a vector of signed longs that should be printed as a ranking.
*	\param out 	The ostream, to which the result should be output. If not defined, cout is used.
*/
void print_compl_ranking(const vector<signed long>* vec, ostream &out = cout){
	
	bool first = true;		// is true until the first positions with 1-bit has been found
	out << "";
	for(size_t i = 0; i < vec->size(); i++){
		
		if(!first){
			out << ", ";
		}
		
		out << "q" << i << ":" << setw(3) << vec->operator[](i);
		first = false;
	}
	out << "";
}


/**
* 	\brief	Method that prints the name for a ComplTag, depending on state-type.
*	
*	Depending on the type of the given ComplTag, different output-formats are used.
*	'p'-states are interpreted as powerset-states, thus bitset stateSet is printed as a set of states.
*	's'-states are interpreted as slice-states, thus the vector slice-vector is printed in slice-notation.
*	'r'-states are interpreted as ranking-states, thus the vector ranking is printed as a ranking, 
*	and additionally the bitset obligationSet is printed as a set of states.
*
*	\param cout	The ostream, to which the result should be output.
*	\param t	The ComplTag, that should be printed.
*/
void print_compl_tag(ostream& cout, ComplTag const& t){

	char stateType = t.stateType;		// t is pointer, -> accesses member of pointer
	
	switch(stateType){
		case 'p': 	// Return the set as state-name for a powerset-state
					cout << "Stateset: ";
					print_compl_bitset(&t.stateSet, cout);
					break;
					
		case 's': 	// Return the slice
					cout << "Slice: ";
					print_compl_slice(&t.slice, cout);
					break;
					
		case 'r': 	// Return the ranking
					cout << "Ranking: ";
					print_compl_ranking(&t.ranking, cout);
					cout << "	   Obligation-Set: ";
					print_compl_bitset(&t.obligationSet, cout);
					break;
					
		default:	// Undefined case, return error
					cout << "State-Type error. Type: " << stateType;
					break;
	}
}


/**
*	\brief	Prints the SCCs and the states they contain for a given SCCDat-pointer via cout.
*
*	\param sccDat Pointer to the SCCDat that should be printed.
*/
void print_SCC(SCCDat* sccDat){
	cout << endl;
	for(state_t i = 0; i < sccDat->sccs.size(); i++){	// For each SCC i
		cout << "SCC " << i << ":   ";
		for(state_t st = 0; st < sccDat->sccs[i].size() ; st++){	// For each state st in SCC i
			cout << "q" << sccDat->sccs[i][st] << " ";
		}
		cout << endl;
	}
}


/**
*	\brief	Prints the elemets of a vector of signed longs with a space between them to a given ostream.
*
*	\param vec	Pointer to the vector<signed long> that contains the elements that should be printed.
*	\param out 	The ostream, to which the result should be output. If not defined, cout is used.
*/
void print_vector(const vector<signed long>* vec, ostream &out = cout){

	for(size_t i = 0; i < vec->size(); i++){
		out << vec->operator[](i) << " ";
	}
}


}	// End of namespace cmpl