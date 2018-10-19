/**
*	\file compl_tests.hh
*	\author Lasse Nitz
*
*	This file contains methods that were used to test different parts of the nbautils-project 
*	and the compl-module.
*	No guarantee is given regarding the functionality and correctness of these test-methods, some might be out-dated.
*/

#pragma once

#include <bitset>
#include <iostream>
#include <string>

#include "aut.hh"
#include "io.hh"
#include "ps.hh"		// For powerset-construction
#include "common/scc.hh"
#include "compl/compl_constr1.hh"
#include "compl/compl_constr2.hh"
#include "compl/compl_constr3.hh"
#include "compl/compl_print.hh"
#include "compl/compl_tag.hh"

using namespace std;
using namespace nbautils;
using namespace cmpl;




void automaton_test1(){

  ////////////////////////////////////////
  // Create an automaton incrementally	//
  ////////////////////////////////////////
  
  cout << "\033[1;4;34m";		// Set custom print color
  cout << "Creating a test-automaton" << endl << endl;
  cout << "\033[0m";			// Reset to default print
  

  // Create an automaton named TestAutomaton with 2 different
  // transition-types (i.e. x, ~x), and initial state 0
  nbautils::Aut<string> aut(true, "TestAutomaton", {"x"}, 0);
  aut.set_pri(0, 1);	// explicitly set initial state 0 to be non-accepting by setting its priority to 1
  
  // Define the generic output function for printing tags
  // To output the tag for a specific state x, call aut.print_state_tag(cout, x)
  aut.tag_to_str = [](ostream& cout, string const& t){
	cout << t;			// t is in this context aut.tag.geti(x) for state x
  };
  
  
  aut.tag.put("q0", 0);	// Give state 0 the name "q0"
  
  cout << "New automaton '" << aut.get_name() << "' with initial state '" << aut.tag.geti(0) << "' created." << endl;  
  cout << "#States: " << aut.num_states() << endl << endl;

  
  
  // Add another state (state 1) and make it a non-final state
  aut.add_state(1);
  aut.set_pri(1, 1);	// Make state 1 a non-final state (by setting priority to 1)
  aut.tag.put("q1", 1);	// Give state 1 the name "q1"
  
  cout << "State '" << aut.tag.geti(1) << "' added." << endl;
  cout << "#States: " << aut.num_states() << endl << endl;

  
  // Add another state (state 2) and make it a final state
  aut.add_state(2);
  aut.set_pri(2, 0);	// Make state 2 a final state (by setting priority to 0)
  aut.tag.put("q2", 2);
  
  cout << "State '" << aut.tag.geti(2) << "' added." << endl;
  cout << "#States: " << aut.num_states() << endl << endl;
  
  

  // Add some transitions, Parameters: (fromState, letter, toState)
  aut.add_edge(0, 0, 0);	// Transition from state 0 to state 0, taking letter 0
  aut.add_edge(0, 1, 0);	// Transition from state 0 to state 0, taking letter 1
  aut.add_edge(0, 1, 1);	// Transition from state 0 to state 1, taking letter 1
  aut.add_edge(1, 1, 2);	// Transition from state 1 to state 2, taking letter 1
  
  cout << "Transitions added." << endl << endl;

  
  
  
  
  ////////////////////////////////////////
  // Property-checks 					//
  ////////////////////////////////////////
  
  // Check if the created automaton is a Buechi-automaton
  cout << "Is the created automaton a Buechi-automaton? - ";
  if(aut.is_buchi()){		cout << "YES" << endl;	}
  else{						cout << "NO" << endl;	}
  
  
  // Check if the created automaton is a deterministic
  cout << "Is the created automaton deterministic? - ";
  if(aut.is_deterministic()){		cout << "YES" << endl;	}
  else{								cout << "NO" << endl;	}
  
  
  // Check if the created automaton is a complete
  cout << "Is the created automaton complete? - ";
  if(aut.is_complete()){	cout << "YES" << endl;	}
  else{						cout << "NO" << endl;	}
  
 
 
  
  
  ////////////////////////////////////////
  // HOA-Output		 					//
  ////////////////////////////////////////
  
  // Output automaton in HOA-format
  cout << "\033[1;4;34m";		// Set custom print color
  cout << endl << endl << "Automaton '" << aut.get_name() <<"' in HOA" << endl <<endl;
  cout << "\033[0m";			// Reset to default print
  print_aut(aut);
  
  
  
  
  
  ////////////////////////////////////////
  // Powerset-Automaton					//
  ////////////////////////////////////////
  
  cout << "\033[1;4;34m";		// Set custom print color
  cout << endl << endl << endl << "Construction of powerset-automaton" << endl << endl;
  cout << "\033[0m";			// Reset to default print
  
  // Construct powerset-automaton of 'TestAutomaton'
  auto psAut = powerset_construction(aut, get_adjmat(aut));
  psAut.set_name("TestAutomaton PS");
  
  cout << "Powerset-automaton created." << endl;
  cout << "#States: " << psAut.num_states() << endl;
  
  cout << "State-list: ";
  for( uint8_t i = 0; i < psAut.num_states(); i++){
	psAut.print_state_tag(cout, i);
	cout << "   ";
  }
  cout << endl;
  
  cout << endl << endl;
  print_aut(psAut);

}




void automaton_test2(){


  // Store pointers to objects of class ComplTag as tags
  nbautils::Aut<ComplTag> aut(true, "TestAutomaton", {"x"}, 0);
  aut.set_pri(0, 1);	// explicitly set initial state 0 to be non-accepting by setting its priority to 1
  
  
  
  ComplTag tag0;
  tag0.stateType = 'p';
  tag0.stateSet.set(0, 1);			// Set bit at position 0 to value 1
  
  // Print out ComlTag, for which the address is put in tag
  cout << "\033[1;4;34m" 	<< "Initial ComplTag, for which the address is saved" << "\033[0m" << endl << endl;
  cout << "StateType: " 	<<  	   tag0.stateType << endl << endl;
  cout << "StateSet: " 		<< endl << tag0.stateSet << endl << endl;
  cout << "ObligationSet: "	<< endl << tag0.obligationSet << endl << endl << endl;
  
  
  aut.tag.put(tag0, 0);		// Save address of tag0 as the tag of state 0
  
  
  ComplTag returnedTag = aut.tag.geti(0);	// Get the ComplTag 
  
  
  
  // Print out ComlTag, for which the address was received from tag
  cout << "\033[1;4;34m" 	<< "ComplTag, for which the address was saved" << "\033[0m" << endl << endl;
  cout << "StateType: " 	<< 		   returnedTag.stateType << endl << endl;
  cout << "StateSet: " 		<< endl << returnedTag.stateSet << endl << endl;
  cout << "ObligationSet: "	<< endl << returnedTag.obligationSet << endl << endl << endl;
  
  
  
  
  // Define output for state-name
  aut.tag_to_str = [](ostream& cout, ComplTag const& t){
	print_compl_tag(cout, t);
  };
  
  
  print_aut(aut);
  
}





void input_test(){

  // Import HOA-file and construct the automaton as datastructure
  nbautils::AutStream<nbautils::Aut<string>> autstream("TestAut.hoa");
  
  //nbautils::MyConsumer<nbautils::Aut<string>> cons();
  
  
  if(autstream.has_next()){				// Check, if there is a next automaton in the HOA-file
	cout << "Next element exists." << endl << endl;
	print_aut(autstream.parse_next());	// Print the next automaton
  }
  else{
	cout << "Nope, garbage!" << endl;
  }
}




void mod_ps_test(){

	nbautils::AutStream<nbautils::Aut<string>> autstream("TestAut.hoa");
	
	if(!autstream.has_next()){
		cout << "\033[31m" << "No next element in inputfile. Aborting test procedure." << "\033[0m" << endl << endl;
		return;
	}
	
	auto aut = autstream.parse_next();

	auto modPsAut = compl_ps_construction(aut, get_adjmat(aut));
	
	/*
	// Set priorities of powerset-state to 1, thus making them non-accepting for a Buechi-automaton
	for(unsigned int i = 0; i < modPsAut.num_states(); i++){
		modPsAut.set_pri(i, 1);
	}
	*/
	
	
	/*
	//////// TEMP TEST///////////////////////////////////////////////////////////////////////////
	list<ComplTag>::iterator it;
	for(it = tagList.begin(); it != tagList.end(); ++it){
		cout << it->stateType << "      " << it->stateSet << endl;
	}
	cout << endl << endl;
	for(unsigned int i = 0; i < modPsAut.num_states(); i++){
		cout << "yello" << endl;
		cout << modPsAut.tag.geti(i)->stateType << endl;//<< "      " << aut.tag.geti(i).stateSet << endl;
		print_compl_bitset(&modPsAut.tag.geti(i)->stateSet);
		cout << endl << endl;
	}
	////////////////////////////////////////////////////////////////////////////////////////////
	*/
	
	cout << "\033[1;34m" 	<< "Powerset-automaton via modified construction: " << "\033[0m" << endl;
	print_aut(modPsAut);
	
	auto psAut = powerset_construction(aut, get_adjmat(aut));
	cout << endl << endl << "\033[1;34m" 	<< "Reference automaton via old construction: " << "\033[0m" << endl;
	print_aut(psAut);
}





void slice_test(){

	nbautils::AutStream<nbautils::Aut<string>> autstream("STS_Aut.hoa");
	
	if(!autstream.has_next()){
		cout << "\033[31m" << "No next element in inputfile. Aborting test procedure." << "\033[0m" << endl << endl;
		return;
	}
	
	auto aut = autstream.parse_next();
	cout << "\033[1;34m" 	<< "Input automaton: " << "\033[0m" << endl;
	print_aut(aut);
	cout << endl;

	auto sliceAut = sts_construction(aut, get_adjmat(aut));

	cout << "\033[1;34m" 	<< "STS: " << "\033[0m" << endl;
	print_aut(sliceAut);
	cout << endl;
	
}




void normalize_slice_test(){
	
	vector<signed long> vec;
	vec.push_back(-1);
	vec.push_back(5);
	vec.push_back(9);
	vec.push_back(2);
	vec.push_back(-1);
	vec.push_back(0);
	auto normVec = normalize_slice(&vec);
	
	cout << endl;
	cout << "\033[1;34m" 	<< "Tightening slice test: " << "\033[0m" << endl;
	print_vector(&vec); 		cout << "		"; print_compl_slice(&vec);
	cout << endl << endl;
	print_vector(&normVec); 	cout << "		"; print_compl_slice(&normVec);
}




void slice_to_rank_test(){
	
	nbautils::Aut<string> aut(true, "TestAutomaton", {"x"}, 0);
	aut.set_pri(0, 1);							// State 0 is 	non-final
	aut.add_state(1);	aut.set_pri(1, 0);		// State 1 is		final
	aut.add_state(2);	aut.set_pri(2, 0);		// State 2 is		final
	aut.add_state(3);	aut.set_pri(3, 1);		// State 3 is	non-final
	aut.add_state(4);	aut.set_pri(4, 0);		// State 4 is		final
	aut.add_state(5);	aut.set_pri(5, 1);		// State 5 is	non-final
	
	
	vector<signed long> slice1;
	slice1.push_back(-1);			// Value for state 0
	slice1.push_back(0);			// Value for state 1
	slice1.push_back(1);			// Value for state 2
	slice1.push_back(2);			// Value for state 3
	slice1.push_back(3);			// Value for state 4
	slice1.push_back(4);			// Value for state 4
	
	auto rank1 = slice_to_rank(aut, slice1);
	
	
	
	
	vector<signed long> slice2;
	slice2.push_back(2);			// Value for state 0
	slice2.push_back(1);			// Value for state 1
	slice2.push_back(1);			// Value for state 2
	slice2.push_back(3);			// Value for state 3
	slice2.push_back(4);			// Value for state 4
	slice2.push_back(0);			// Value for state 4
	
	auto rank2 = slice_to_rank(aut, slice2);
	
	
	
	
	// HERE: TEST SLICE TO RANK
	cout << "\033[1;34m" 	<< "Slice to rank test: " << "\033[0m" << endl << endl;
	
	cout << "Non-final states: ";
	for(auto& i : aut.states()){
		if(aut.get_pri(i) == 1)
			cout << "q" << i << "   ";
	}
	cout << endl << endl;
	
	print_compl_slice(&slice1);
	cout << " ------> ";
	print_vector(&rank1);
	cout << endl << endl;
	
	print_compl_slice(&slice2);
	cout << " ------> ";
	print_vector(&rank2);
	cout << endl << endl;
	
	
	// Example from thesis
	aut.set_pri(0, 0);			// State 0 is 		final
	aut.set_pri(1, 1);			// State 1 is 	non-final
	aut.set_pri(2, 1);			// State 2 is 	non-final
	aut.set_pri(3, 0);			// State 3 is 		final
	aut.set_pri(4, 0);			// State 4 is 		final
	aut.set_pri(5, 1);			// State 5 is 	non-final
	
	cout << endl << endl << "\033[1;34m" << "Example from thesis" << "\033[0m" << endl << endl << "Non-final states: ";
	for(auto& i : aut.states()){
		if(aut.get_pri(i) == 1)
			cout << "q" << i << "   ";
	}
	cout << endl << endl;
	
	vector<signed long> slice3;
	slice3.push_back(0);			// Value for state 0
	slice3.push_back(1);			// Value for state 1
	slice3.push_back(1);			// Value for state 2
	slice3.push_back(2);			// Value for state 3
	slice3.push_back(3);			// Value for state 4
	slice3.push_back(4);			// Value for state 4
	
	auto rank3 = slice_to_rank(aut, slice3);
	
	print_compl_slice(&slice3);
	cout << " ------> ";
	print_vector(&rank3);
	cout << endl << endl;
}





void tighten_ranking_test(){
	nbautils::Aut<string> aut(true, "TestAutomaton", {"x"}, 0);
	aut.set_pri(0, 1);							// State 0 is 	non-final
	aut.add_state(1);	aut.set_pri(1, 0);		// State 1 is		final
	aut.add_state(2);	aut.set_pri(2, 0);		// State 2 is		final
	aut.add_state(3);	aut.set_pri(3, 1);		// State 3 is	non-final
	aut.add_state(4);	aut.set_pri(4, 0);		// State 4 is		final
	aut.add_state(5);	aut.set_pri(5, 1);		// State 5 is	non-final
	
	
	vector<signed long> rank1;
	rank1.push_back(-1);		// Value for state 0
	rank1.push_back(0);			// Value for state 1
	rank1.push_back(5);			// Value for state 2
	rank1.push_back(4);			// Value for state 3
	rank1.push_back(3);			// Value for state 4
	rank1.push_back(2);			// Value for state 4
	
	auto tRank1 = tighten_ranking(aut, &rank1);
	
	
	
	cout << "\033[1;34m" 	<< "Tighten ranking test: " << "\033[0m" << endl << endl;
	
	cout << "Non-final states: ";
	for(auto& i : aut.states()){
		if(aut.get_pri(i) == 1)
			cout << "q" << i << "   ";
	}
	cout << endl << endl;
	
	print_vector(&rank1);
	cout << " ------> ";
	print_vector(&tRank1);
	cout << endl << endl;
}




void successor_ranking_test(){

	nbautils::AutStream<nbautils::Aut<string>> autstream("STS_Aut.hoa");
	
	if(!autstream.has_next()){
		cout << "\033[31m" << "No next element in inputfile. Aborting test procedure." << "\033[0m" << endl << endl;
		return;
	}
	
	auto aut = autstream.parse_next();
	cout << "\033[1;34m" 	<< "Input automaton: " << "\033[0m" << endl;
	print_aut(aut);
	cout << endl;
	
	vector<signed long> ranking1;
	ranking1.push_back(1);
	ranking1.push_back(0);
	ranking1.push_back(0);
	
	bitset<max_nba_states> obSet;
	obSet.set(1);
	obSet.set(2);
	
	ComplTag tag;
	tag.ranking = ranking1;
	tag.obligationSet = obSet;
	
	auto sucTag = succ_ranking(aut, tag, 1);
	
	cout << "\033[1;34m" 	<< "Successor ranking test: " << "\033[0m" << endl << endl;
	
	cout << endl;
	print_compl_ranking(&sucTag.ranking);
	cout << endl;
	print_compl_bitset(&sucTag.obligationSet);
}





void scc_test(){
	nbautils::AutStream<nbautils::Aut<string>> autstream("SCC_Aut.hoa");
	
	if(!autstream.has_next()){
		cout << "\033[31m" << "No next element in inputfile. Aborting test procedure." << "\033[0m" << endl << endl;
		return;
	}
	
	auto aut = autstream.parse_next();
	cout << "\033[1;34m" 	<< "Input automaton: " << "\033[0m" << endl;
	print_aut(aut);
	cout << endl;
	
	auto aut_comp = compl_construction(aut, get_adjmat(aut));
	cout << "\033[1;34m" 	<< "Complementary automaton A_comp: " << "\033[0m" << endl;
	print_aut(aut_comp);
	cout << endl;
	
	// Find SCCs in the input automaton
	SCCDat sccs = get_sccs(aut_comp.states(), aut_succ(aut_comp) , true);
	print_SCC(&sccs);
	cout << endl;
	
}





void last_scc_test(){

	cout << "\033[1;34m" 	<< "Last SCC test: " << "\033[0m" << endl;

	nbautils::Aut<ComplTag> aut(true, "TestAutomaton", {"x"}, 0);
	aut.set_pri(0, 1);	// explicitly set initial state 0 to be non-accepting by setting its priority to 1

  aut.tag_to_str = [](ostream& out, ComplTag const& t){
    print_compl_tag(out, t);
  };

  // Add states
  aut.add_state(1);
  aut.set_pri(1, 0);
  
  aut.add_state(2);
  aut.set_pri(2, 0);	
  
  aut.add_state(3);
  aut.set_pri(3, 1);	
  
  aut.add_state(4);
  aut.set_pri(4, 1);	 
  
  aut.add_state(5);
  aut.set_pri(5, 1);	
  
  aut.add_state(6);
  aut.set_pri(6, 1);	
  
  aut.add_state(7);
  aut.set_pri(7, 1);	
  
  aut.add_state(8);
  aut.set_pri(8, 0);	
  
  
  ComplTag tag0;
  ComplTag tag1;
  ComplTag tag2;
  ComplTag tag3;
  ComplTag tag4;
  ComplTag tag5;
  ComplTag tag6;
  ComplTag tag7;
  ComplTag tag8;
  
  tag0.stateSet = 0;	//0
  tag1.stateSet = 1;		//1
  tag2.stateSet = 3;				//3
  tag3.stateSet = 2;			//2
  tag4.stateSet = 3;				//3
  tag5.stateSet = 3;				//3
  tag6.stateSet = 0;	//0
  tag7.stateSet = 2;			//2
  tag8.stateSet = 3;				//3
  
  tag0.stateType = 'p';
  tag1.stateType = 'p';
  tag2.stateType = 'p';
  tag3.stateType = 'p';
  tag4.stateType = 'p';
  tag5.stateType = 'p';
  tag6.stateType = 'p';
  tag7.stateType = 'p';
  tag8.stateType = 'p';
  
	
  aut.tag.put(tag0, 0);
  aut.tag.put(tag1, 1);		
  aut.tag.put(tag2, 2);		
  aut.tag.put(tag3, 3);	
  aut.tag.put(tag4, 4);	
  aut.tag.put(tag5, 5);	
  aut.tag.put(tag6, 6);		
  aut.tag.put(tag7, 7);		
  aut.tag.put(tag8, 8);
  

  
  
  // Add transitions; a=0, b=1
  aut.add_edge(0, 0, 5);
  aut.add_edge(0, 1, 1);
  
  aut.add_edge(1, 0, 2);
  aut.add_edge(1, 1, 1);
  
  aut.add_edge(2, 0, 3);
  aut.add_edge(2, 1, 1);
  
  aut.add_edge(3, 0, 4);
  aut.add_edge(3, 1, 4);
  
  aut.add_edge(4, 0, 4);
  aut.add_edge(4, 1, 4);
  
  aut.add_edge(5, 0, 6);
  aut.add_edge(5, 1, 7);
  
  aut.add_edge(6, 0, 6);
  aut.add_edge(6, 1, 8);
  
  aut.add_edge(7, 0, 8);
  aut.add_edge(7, 1, 7);
  
  aut.add_edge(8, 0, 8);
  aut.add_edge(8, 1, 8);
  
  
  print_aut(aut);
  
  SCCDat sccDat = get_sccs(aut.states(), aut_succ(aut) , true);
  print_SCC(&sccDat);
  
  
  for(unsigned i = 0; i < 4; i++){
	cout << endl << "Last SCCs for " << i <<":  ";
	auto vec = last_sccs_for_stateset(aut, i);
	for(unsigned j = 0; j < vec.size(); j++){
		cout << vec[j] << "   ";
	}
	cout << endl;
  }
}






void write_to_file_test(){
	// Read automaton from file
	nbautils::AutStream<nbautils::Aut<string>> autstream("SCC_Aut.hoa");
	
	if(!autstream.has_next()){
		cout << "\033[31m" << "No next element in inputfile. Aborting test procedure." << "\033[0m" << endl << endl;
		return;
	}
	
	auto aut = autstream.parse_next();
	auto aut_comp = compl_construction(aut, get_adjmat(aut));
	
	// Output A_compl to "Output.hoa"
	ofstream output;
	output.open("Output.hoa");
	print_aut(aut_comp, output);
	output.close();
}




void kill_test(){

	nbautils::Aut<string> aut(true, "TestAutomaton", {"x"}, 0);
	aut.set_pri(0, 1);	// explicitly set initial state 0 to be non-accepting by setting its priority to 1

  aut.tag_to_str = [](ostream& cout, string const& t){
	cout << t;			// t is in this context aut.tag.geti(x) for state x
  };
  
  
  aut.tag.put("q0", 0);	// Give state 0 the name "q0"

  
  // Add states
  aut.add_state(1);
  aut.set_pri(1, 0);
  aut.tag.put("q1", 1);	
  
  aut.add_state(2);
  aut.set_pri(2, 0);	
  aut.tag.put("q2", 2);
  
  aut.add_state(3);
  aut.set_pri(3, 1);	
  aut.tag.put("q3", 3);
  
  cout << endl << "BEFORE DELETION:" << endl;
  print_aut(aut);
  
  vector<state_t> killVec;
  killVec.push_back(2);
  
  cout << endl;
  for(unsigned i : aut.states()){
	cout << "q" << i << "  ";
  }
  
  aut.remove_states(killVec);
  
  cout << endl << "Correct: ";
  for(unsigned i : aut.states()){
	cout << "q" << i << "  ";
  }
  
  cout << endl << "Incorrect: ";
  for(unsigned i = 0; i < aut.num_states(); i++){
	cout << "q" << i << "  ";
  }
  
  cout << endl << "Number of states: " << aut.num_states() << endl;
  
  cout << "AFTER DELETION: " << endl;
  print_aut(aut);

}






void compl_AL_test(){
	nbautils::AutStream<nbautils::Aut<string>> autstream("STS_Aut.hoa");
	
	if(!autstream.has_next()){
		cout << "\033[31m" << "No next element in inputfile. Aborting test procedure." << "\033[0m" << endl << endl;
		return;
	}
	
	auto aut = autstream.parse_next();
	cout << "\033[1;34m" 	<< "Input automaton: " << "\033[0m" << endl;
	print_aut(aut);
	cout << endl;

	auto autAL = al_construction(aut, get_adjmat(aut));

	cout << "\033[1;34m" 	<< "Complementary automaton AL: " << "\033[0m" << endl;
	print_aut(autAL);
	cout << endl;
}





void compl_Acomp_test(){
	nbautils::AutStream<nbautils::Aut<string>> autstream("STS_Aut.hoa");
	
	if(!autstream.has_next()){
		cout << "\033[31m" << "No next element in inputfile. Aborting test procedure." << "\033[0m" << endl << endl;
		return;
	}
	
	auto aut = autstream.parse_next();
	//cout << "\033[1;34m" 	<< "Input automaton: " << "\033[0m" << endl;
	//print_aut(aut);
	cout << endl;

	auto aut_comp = compl_construction(aut, get_adjmat(aut));

	cout << "\033[1;34m" 	<< "Complementary automaton A_comp: " << "\033[0m" << endl;
	print_aut(aut_comp);
	cout << endl;
}






void create_some_aut(){
	nbautils::Aut<string> aut(true, "TestAutomaton", {"x"}, 0);
	aut.set_pri(0, 1);	// explicitly set initial state 0 to be non-accepting by setting its priority to 1

  aut.tag_to_str = [](ostream& cout, string const& t){
	cout << t;			// t is in this context aut.tag.geti(x) for state x
  };
  
  
  aut.tag.put("q0", 0);	// Give state 0 the name "q0"

  
  // Add states
  aut.add_state(1);
  aut.set_pri(1, 0);
  aut.tag.put("q1", 1);	
  
  aut.add_state(2);
  aut.set_pri(2, 0);	
  aut.tag.put("q2", 2);
  
  aut.add_state(3);
  aut.set_pri(3, 1);	
  aut.tag.put("q3", 3);
  
  aut.add_state(4);
  aut.set_pri(4, 1);	
  aut.tag.put("q4", 4);
  
  
  aut.add_state(5);
  aut.set_pri(5, 1);	
  aut.tag.put("q5", 5);
  
  aut.add_state(6);
  aut.set_pri(6, 1);	
  aut.tag.put("q6", 6);
  
  aut.add_state(7);
  aut.set_pri(7, 1);	
  aut.tag.put("q7", 7);
  
  aut.add_state(8);
  aut.set_pri(8, 0);	
  aut.tag.put("q8", 8);
  
  
  

  // Add transitions; a=0, b=1
  aut.add_edge(0, 0, 5);
  aut.add_edge(0, 1, 1);
  
  aut.add_edge(1, 0, 2);
  aut.add_edge(1, 1, 1);
  
  aut.add_edge(2, 0, 3);
  aut.add_edge(2, 1, 1);
  
  aut.add_edge(3, 0, 4);
  aut.add_edge(3, 1, 4);
  
  aut.add_edge(4, 0, 4);
  aut.add_edge(4, 1, 4);
  
  aut.add_edge(5, 0, 6);
  aut.add_edge(5, 1, 7);
  
  aut.add_edge(6, 0, 6);
  aut.add_edge(6, 1, 8);
  
  aut.add_edge(7, 0, 8);
  aut.add_edge(7, 1, 7);
  
  aut.add_edge(8, 0, 8);
  aut.add_edge(8, 1, 8);
  
  
  print_aut(aut);
}














void buchi_compl(){	

	// Set up input-stream for a HOA-file that contains the input automata
	nbautils::AutStream<nbautils::Aut<string>> autstream("TestAut.hoa");
	
	// Check if there does exist a first element in the given file. Abort, if no first element exists.
	if(!autstream.has_next()){
		cout << "\033[31m" << "No elements in inputfile. Aborting test procedure." << "\033[0m" << endl << endl;
		return;
	}
	
	ofstream output_AL;
	output_AL.open("Output_AL.hoa");
	
	ofstream output_A_comp;
	output_A_comp.open("Output_A_comp.hoa");
	
	ofstream output_A_opt;
	output_A_opt.open("Output_A_opt.hoa");
	
	unsigned long c = 0;
	
	while(autstream.has_next()){
		
		// Input automaton
		auto aut = autstream.parse_next();
		auto adj = get_adjmat(aut);
		
		// Create complementary automata
		auto AL = 		al_construction			(aut, adj);
		auto A_comp = 	compl_construction		(aut, adj);
		auto A_opt = 	compl_construction_opt	(aut, adj);
		
		// Give names to complementary automata for increased readability regarding console output
		AL.set_name("AL #" + to_string(c));
		A_comp.set_name("A_comp #" + to_string(c));
		A_opt.set_name("A_opt #" + to_string(c));
		
		// Print to output files
		print_aut(AL, output_AL);
		print_aut(A_comp, output_A_comp);
		print_aut(A_opt, output_A_opt);
		
		// Print to console
		cout << endl << endl << "\033[1;34m" << "Next Automaton read" << "\033[0m" << endl << endl;
		print_aut(AL);
		cout << endl << endl;
		print_aut(A_comp);
		cout << endl << endl;
		print_aut(A_opt);
		cout << endl;
		
		c++;
	}

}





void compl_test(){
	// Test Methods 
	
	//automaton_test1();		// Test method that builds an automaton and performs some operations on it
	//automaton_test2();		// Test method that builds an automaton using custom ComplTag-type
	//input_test();			// Test method reading input-file containing HOA-parity-automatons
	//mod_ps_test();			// Test method that creates a powerset-automaton with the modified method
	//slice_test();				// Test method that creates a Slice-Transition-System for some automaton
	//normalize_slice_test();	// Creates a slice and performs the tighten function on it to make it unique for that order
	//slice_to_rank_test();		// Test that generates rankings from slices
	//tighten_ranking_test();	// Test method to check whether a tight ranking is generated correctly for an arbitrary given ranking
	//successor_ranking_test();
	compl_AL_test();
	compl_Acomp_test();
	scc_test();
	write_to_file_test();
	//kill_test();				// Shows problems of deleting states from an automaton without normalizing the result
	last_scc_test();			// Determines "last SCCs" for statesets appearing in an automaton using ComplTags
}





// Written to find a bug
void rankTEST(){

	// Set up input-stream for a HOA-file that contains the input automata
	nbautils::AutStream<nbautils::Aut<string>> autstream("error.hoa");
	
	// Check if there does exist a first element in the given file. Abort, if no first element exists.
	if(!autstream.has_next()){
		cout << "\033[31m" << "No elements in inputfile. Aborting test procedure." << "\033[0m" << endl << endl;
		return;
	}
	
	auto aut = autstream.parse_next();
	
	ComplTag tag;
	tag.stateType = 'r';
	
	tag.stateSet = 3;
	
	vector<signed long> ranking;
	ranking.push_back(0);
	ranking.push_back(1);
	
	tag.ranking = ranking;
	tag.obligationSet = 0;
	
	
	
	auto aTrans = succ_ranking(aut, tag, 0);
	auto bTrans = succ_ranking(aut, tag, 1);
	
	auto tight_aTrans = tighten_ranking(aut, &aTrans.ranking);
	auto tight_bTrans = tighten_ranking(aut, &bTrans.ranking);
	
	print_compl_ranking(&aTrans.ranking);
	cout << endl << endl;
	print_compl_ranking(&bTrans.ranking);
	cout << endl << endl;
	print_compl_ranking(&tight_aTrans);
	cout << endl << endl;
	print_compl_ranking(&tight_bTrans);
	cout << endl << endl;
	
	
}






