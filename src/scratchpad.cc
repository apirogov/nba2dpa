#include <iostream>
#include <string>
using namespace std;

#include <spdlog/spdlog.h>
#include <args.hxx>
#include "dpa.hh"
#include "nba.hh"

auto parse_args(int argc, char *argv[]) {
  args::ArgumentParser parser("This is a test program.", "");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::Positional<string> input(parser, "INPUTFILE", "File containing the NBA");
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
    cerr << "Please provide an input file!" << endl;
    exit(1);
  }

  return args::get(input);
}

int main(int argc, char *argv[]) {
  static auto logger = spdlog::stdout_color_mt("times_logger");

  auto file = parse_args(argc, argv);

  logger->debug("Reading automaton from {}", file);

  auto aut = nbautils::spot_nba_from_file(file);
  // in principle here we could do some trimming using spots default stuff

  auto nba = nbautils::NBA(aut);

  auto dpa = nbautils::DPA(nba);

  // auto sccs = spot::scc_info(aut);

  // print NBA again
  for (auto s : nba.states) {
    cout << s;
    if (nba.acc.find(s) != nba.acc.end()) cout << "*";
    cout << endl;
    for (auto i = 0; i < nba.num_syms; i++) {
      if (nba.adj[s][i].empty()) continue;
      cout << "\t" << i << " -> ";
      for (auto t : nba.adj[s][i]) cout << t << ", ";
      cout << endl;
    }
  }

  // auto ps = vector<nbautils::NBA::state_t>{0, 5};
  // for (auto q : nba.powersucc(ps, 0)) cout << q << endl;
}
