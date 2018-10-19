/**
*	\file compl.cc
*	\author Lasse Nitz
*
*	This file contains the main-method of the compl-module, which deals with
*	the complementation of Buechi-automata.
*
*	It also deals with the processing of commandline-arguments.
*/

// C++ Standard Library
#include <cstdint>
#include <iostream>
#include <string>

// Third-party
#include <args.hxx>

// nbautils
#include "aut.hh"
#include "io.hh"

// Compl
#include "compl/compl_constr1.hh"
#include "compl/compl_constr2.hh"
#include "compl/compl_constr3.hh"
#include "compl/compl_print.hh"
#include "compl/compl_tag.hh"
#include "compl/compl_tests.hh"		// Test methods


using namespace std;		// C++ standard namespace
using namespace nbautils;	// nbautils-project namespace
using namespace cmpl;		// compl-module namespace

/**
*	\brief Datatype that saves data related to input-arguments of the main-function.
*/
struct Args{
	string input_path;	/// Path including filename

	bool all;			/// True, if the input-file contains several HOA automata that should be all processed
	bool no_cout;		/// True, if the output should not be printed to console
	bool out;			/// True, if the output should be additionally stored in a file
	string output_path;	/// The path of the output-file
	bool stats;			/// True, if statistics should be output

	bool constr1;		/// Construction using STS and rankings
	bool constr2;		/// Construction using PS and rankings
	bool constr3;		/// Construction using PS and rankings with optimizations
};




/**
*	\brief Parses commandline-arguments of the main-function and returns an Args-object in which the according variables are set.
*
*	This method uses the third-party library that is included via args.hxx	(see https://github.com/Taywee/args).
*
*	\param argc		The amount of commandline-arguments (Argument Count).
*	\param argv		The pointer to the first element of the array containing the commandline-arguments (Argument Vector).
*
*	\return		An Args-Object in which the variables are set according to the commandline-arguments.
*/
Args parse_args(int argc, char* argv[]){

	string headerMessage = "Tool for the complementation of Büchi-automata with up to ";
	headerMessage += to_string(max_nba_states);
	headerMessage += " states and up to ";
	headerMessage += to_string(max_nba_syms);
	headerMessage += " letters in the alphabet.";

	args::ArgumentParser parser(headerMessage, "");
	args::Positional<string> input_path(parser, "input_file", "Path and name of a file that contains HOA-automata. If not specified, the default value \"Input.hoa\" will be used.");
	args::HelpFlag help(parser, "help", "Displays this help menu.", {'h', "help"});


	args::Flag all(parser, "all",
		"Choose this option if the input-file contains several HOA-automata. By default, only the first automaton in the input-file is read. If this option is chosen, the tool will create and output the complementary automata via the chosen constructions, before it reads the next input-automaton.",
		{'a', "all"});

	args::Flag no_cout(parser, "no_cout",
		"Choose this option if the output should not be printed to console. If an output-file is specified, the output will still be saved in that file.",
		{'n', "nocout"});

	args::ValueFlag<string> output_path(parser, "output_file",
		"Additionally to the console-output, the output of the calculation is stored in an output-file, that is specified by the given path. If this file already exists, it will be overwritten.",
		{'o', "output"});

	args::Flag stats(parser, "stats",
		"Choose this option if only the statistics of the resulting automata for the selected constructions are relevant. The automata will not be output.",
		{'s', "stats"});

	// Complementation constructions
	args::Group constr(parser,
		"The following options are construction-approaches for complementary Büchi-automata.  In case that several construction methods are chosen, the output will be ordered according to: constr1, constr2, constr3.");

	args::Flag constr1(constr,
		"constr1", "Calls the complementation approach that uses the Slice-Transition-System and rankings.",
		{"constr1"});

	args::Flag constr2(constr, "constr2",
		"Calls the complementation approach that uses the Powerset-automaton and rankings.",
		{"constr2"});

	args::Flag constr3(constr, "constr3",
		"Calls the optimized complementation approach that uses the Powerset-automaton and rankings.",
		{"constr3"});



	// Try to parse, check for exit-cases
	try{
		parser.ParseCLI(argc, argv);
	}catch (args::Help&) {
		std::cout << parser;
		exit(0);
	}catch (args::ParseError& e) {
		cerr << e.what() << endl << parser;
		exit(1);
	}catch (args::ValidationError& e) {
		cerr << e.what() << endl << parser;
		exit(1);
	}


	// Create return-value and fill it
	Args res;

	// Save the path of the input-file
	res.input_path = args::get(input_path);

	// Output file
	string outfile = args::get(output_path);
	ofstream exists2(outfile);
	if(output_path && !exists2){
		cout << endl << "Writing to output-file  \"" << outfile << "\" failed." << endl << endl;
		cout << "Please make sure that the given directory exists and that a filename is specified." << endl;
		exit(1);
	}
	res.out = output_path;
	res.output_path = outfile;


	// Fill in remaining values
	res.all = all;
	res.no_cout = no_cout;
	res.stats = stats;

	res.constr1 = constr1;
	res.constr2 = constr2;
	res.constr3 = constr3;


	return res;
}




/**
*	\brief 	Returns the number of transitions for a given automaton.
*
*	The automaton that is given via the parameter aut is assumed to have less than 2^64 transitions.
*
*	\param aut	The automaton for which the number of transitions should be returned.
*
*	\return 	The number of transitions of the automaton that was given as a parameter.
*/
uint64_t num_trans(auto const& aut){
	uint64_t res = 0;

	for(auto i : aut.states()){
		for(auto x : aut.syms()){
			res = res + aut.succ(i, x).size();
		}
	}

	return res;
}


/**
*	\brief main-function of the compl-module.
*
*	\param argc		The amount of commandline-arguments (Argument Count).
*	\param argv		The pointer to the first element of the array containing the commandline-arguments (Argument Vector).
*
*	\return	0, if the program executed without unexpected problems.
*/
int main(int argc, char* argv[]){

	// Parse arguments
	auto const args = parse_args(argc, argv);

	// If no construction-method is chosen, return 0
	if(!args.constr1 && !args.constr2 && !args.constr3){
		cerr << "No construction method chosen." << endl;
		cerr << "Please choose a construction method, details can be found via '-h' or '--help'." << endl;
		return 0;
	}

	// If this part of the code is reached, args.input_path is defined and an according file does exist
	nbautils::AutStream<nbautils::Aut<string>> autstream(args.input_path);
	Aut<string> aut;

	// Output stream
	ofstream output;
	if(args.out){
		output.open(args.output_path);
	}

	// Check for existence of an automaton in input-file
	if(!autstream.has_next()){
		cerr << "File  \"" << args.input_path << "\"  does not contain an automaton.";
		exit(1);
	}

	// Write meaning of columns for statistics
	if(args.stats){
		if(args.constr1){
			if(!args.no_cout){
				cout << "#states CONSTR1, #transitions CONSTR1";
				if(args.constr2 || args.constr3){ cout << ", "; }	// Delimiter if more stats are added in this output-line
			}
			if(args.out){
				output << "#states CONSTR1, #transitions CONSTR1";
				if(args.constr2 || args.constr3){ output << ", "; }	// Delimiter if more stats are added in this output-line
			}
		}
		if(args.constr2){
			if(!args.no_cout){
				cout << "#states CONSTR2, #transitions CONSTR2";
				if(args.constr3){ cout << ", "; }	// Delimiter if more stats are added in this output-line
			}
			if(args.out){
				output << "#states CONSTR2, #transitions CONSTR2";
				if(args.constr3){ output << ", "; }	// Delimiter if more stats are added in this output-line
			}
		}
		if(args.constr3){
			if(!args.no_cout){cout << "#states CONSTR3, #transitions CONSTR3";}
			if(args.out){output << "#states CONSTR3, #transitions CONSTR3";}
		}

		if(!args.no_cout){cout << endl;}
		if(args.out){output << endl;}
	}

	// Execution according to arguments
	do {
		aut = autstream.parse_next();
        aut.make_colored();

		// Check whether the input-automaton has too many states
		if(aut.states().size() > max_nba_states){
			cerr << "The input-automaton has too many states. Please make sure that the automaton has at most " << max_nba_states << " states." << endl;
			return 0;
		}

		// Check whether the input-automaton has too many letters in the alphabet
		if(aut.syms().size() > max_nba_syms){
			cerr << "The input-automaton has too many letters. Please make sure that the alphabet has at most " << max_nba_syms << " letters." << endl;
			return 0;
		}


		if(args.constr1){
			auto const AL = al_construction(aut, get_adjmat(aut));

			// Statistics-output
			if(args.stats){
				if(!args.no_cout){
					cout << AL.num_states() << ", " << num_trans(AL);
					if(args.constr2 || args.constr3){ cout << ", "; }		// Delimiter if more stats are added in this output-line
				}
				if(args.out){
					output  << AL.num_states() << ", " << num_trans(AL);
					if(args.constr2 || args.constr3){ output << ", "; }		// Delimiter if more stats are added in this output-line
				}
			}
			else{
				// Automaton-output on console
				if(!args.no_cout){ print_aut(AL);}

				// Output to custom file
				if(args.out){ print_aut(AL, output); }
			}
		}

		if(args.constr2){
			auto const Acomp = compl_construction(aut, get_adjmat(aut));

			// Statistics-output
			if(args.stats){
				if(!args.no_cout){
					cout << Acomp.num_states() << ", " << num_trans(Acomp);
					if(args.constr3){ cout << ", "; }		// Delimiter if more stats are added in this output-line
				}
				if(args.out){
					output << Acomp.num_states() << ", " << num_trans(Acomp);
					if(args.constr3){ output << ", "; }		// Delimiter if more stats are added in this output-line
				}
			}
			else{
				// Automaton-output on console
				if(!args.no_cout){ print_aut(Acomp); }

				// Output to custom file
				if(args.out){ print_aut(Acomp, output); }
			}
		}

		if(args.constr3){
			auto const Aopt = compl_construction_opt(aut, get_adjmat(aut));

			// Statistics-output
			if(args.stats){
				if(!args.no_cout){ cout << Aopt.num_states() << ", " << num_trans(Aopt); }
				if(args.out){ output << Aopt.num_states() << ", " << num_trans(Aopt);}
			}
			else{
				// Automaton-output on console
				if(!args.no_cout){ print_aut(Aopt); }

				// Output to custom file
				if(args.out){ print_aut(Aopt, output); }
			}
		}

		if(args.stats){
			if(!args.no_cout){ cout << endl; }
			if(args.out){ output << endl;}
		}

	} while (args.all && autstream.has_next());

	output.close();		// Close stream of output-file

	return 0;
}
