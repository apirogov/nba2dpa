#include <fstream>
#include <iostream>
#include <string>
using namespace std;

// #include <spdlog/spdlog.h>
#include "cpphoafparser/consumer/hoa_consumer_print.hh"
#include "cpphoafparser/parser/hoa_parser.hh"
#include "io.hh"
#include "ps.hh"
#include "scc.hh"
#include "types.hh"
#include "debug.hh"
#include <args.hxx>
using namespace cpphoafparser;

auto parse_args(int argc, char *argv[]) {
  args::ArgumentParser parser("This is a test program.", "");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::Positional<string> input(parser, "INPUTFILE",
                                 "File containing the NBA");
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help) {
    std::cout << parser;
    exit(0);
  } catch (args::ParseError e) {
    cerr << e.what() << endl << parser;
    exit(1);
  } catch (args::ValidationError e) {
    cerr << e.what() << endl << parser;
    exit(1);
  }

  if (!input) {
    return string("");
  }

  return args::get(input);
}

int main(int argc, char *argv[]) {
  auto file = parse_args(argc, argv);

  auto pnba = nbautils::parse_ba(file);
  if (!pnba.second.success) {
    cerr << "no valid NBA parsed!" << endl;
	exit(1);
  }

	printMeta(pnba.second);
	auto& ba = *pnba.first;
    cout << "Input BA has " << ba.adj.size() << " states" << endl;
	auto scci = nbautils::get_scc_info(ba, true);
    printSCCI(*scci);
	cout << "trimming BA.." << endl;
	cout << "removed " << trim_ba(ba, *scci) << " states" << endl;
	printBA(ba, *scci);
    printSCCI(*scci);


  auto ps = nbautils::powerset_product(ba);
  auto pscci = nbautils::get_scc_info(*ps, true);
  // printPS(ps, pscci, true);
  // printPS(ps);
  cout << "PS size: " << ps->adj.size() << " with " <<
	  pscci->sccrep.size() << " SCCs" << endl;
  // printSCCI(pscci);

  cout << "trimming PS.." << endl;
	cout << "removed " << trim_ba(*ps, *pscci) << " states" << endl;
  // printPS(ps, pscci, true);
}
