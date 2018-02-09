#include <fstream>
#include <iostream>
#include <string>
using namespace std;

#include "types.hh"
#include "io.hh"
#include "scc.hh"
#include "ps.hh"
#include "debug.hh"

#include <args.hxx>

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
  if (!pnba) {
    cerr << "no valid NBA parsed!" << endl;
	exit(1);
  }

}
