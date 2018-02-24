#include <iostream>
#include <string>
using namespace std;

#include <spdlog/spdlog.h>
namespace spd = spdlog;

#include <args.hxx>

#include "io.hh"
#include "scc.hh"


using namespace nbautils;

struct Args {
  using uptr = std::unique_ptr<Args>;

  string file;
  int verbose;

  int nodeinfo;
  int edgeinfo;
};

Args::uptr parse_args(int argc, char *argv[]) {
  args::ArgumentParser parser("nba2dot - visualize nondeterministic Büchi automata", "");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});

  // exactly one input automaton
  args::Positional<string> input(parser, "INPUTFILE", "file containing the NBA (if none given, uses stdin)");

  // logging level -v, -vv, etc.
  args::CounterFlag verbose(parser, "verbose", "Show verbose information", {'v', "verbose"});
  args::ValueFlag<int> nodeinfo(parser, "nodeinfo", "Node info level", {'n', "node-info"});
  args::ValueFlag<int> edgeinfo(parser, "edgeinfo", "Edge info level", {'e', "edge-info"});
  // args::Flag trim(parser, "trim", "Remove dead states from NBA", {'d', "trim"});

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


  auto args = make_unique<Args>(Args());
  if (input) args->file = args::get(input);
  args->verbose = args::get(verbose);
  args->nodeinfo = args::get(nodeinfo);
  args->edgeinfo = args::get(edgeinfo);
  // args->trim = trim;

  return move(args);
}

void dot_header(int fontsize=12, string fontname="Verdana", bool compound=true) {
  cout << "digraph G {" << endl;
  cout << "graph [fontsize=" << fontsize << ", fontname=\"" << fontname << "\", compound=" << (compound ?"true":"false") << "]" << endl;
  cout << "node [style=\"filled\",fillcolor=\"white\",fontsize=" << fontsize << ", fontname=\"" << fontname << "\"]" << endl;

}
void dot_footer() { cout << "}" << endl; }

void dot_state(state_t s, BA const& aut, set<state_t> const& unreach,
       set<state_t> const& dead, set<state_t> const& accsinks, int nodeinfo=0) {
  vector<string> ps;
  auto initial = set<state_t>(cbegin(aut.get_init()), cend(aut.get_init()));

  if (nodeinfo==0)
    ps.push_back("style=\"invis\"");
  if (nodeinfo==1)
    ps.push_back("width=\"0.2\",fixedwidth=true,label=\"\"");
  if (nodeinfo<=2)
    ps.push_back("shape=\"circle\"");
  if (nodeinfo==3)
    ps.push_back("label=\""+ aut.print_tag(s) +"\"");
  if (aut.has_accs(s))
    ps.push_back("peripheries=2");
  if (aut.has_accs(s) || contains(initial, s))
    ps.push_back("penwidth=2");
  if (contains(initial, s))
    ps.push_back("color=\"blue\"");
  bool useless = contains(dead,s) || contains(unreach,s);
  if (useless)
    ps.push_back("fillcolor=\"red\"");
  if (!useless && contains(accsinks, s))
    ps.push_back("fillcolor=\"green\"");

  cout << s << "[";
  for (int i=0; i<(int)ps.size(); i++)
    cout << ps[i] << (i<(int)ps.size()-1 ? "," : "");
  cout << "];" << endl;
}

int main(int argc, char *argv[]) {
  // initialize stuff (args + logging)

  // auto console = spd::stdout_color_mt("log");
  auto log = spd::stdout_logger_mt("log");
  spd::set_pattern("[%Y-%m-%d %H:%M:%S %z] [%l] %v");

  auto args = parse_args(argc, argv);
  if (!args->verbose)
    spd::set_level(spd::level::warn);
  else if (args->verbose == 1)
    spd::set_level(spd::level::info);
  else
    spd::set_level(spd::level::debug);

  // now parse input automaton
  auto auts = nbautils::parse_hoa_ba(args->file, log);

  if (auts.empty()) {
    log->error("Parsing NBAs from {} failed!", args->file.empty() ? "stdin" : args->file);
    exit(1);
  }

  if (auts.size() > 1)
    log->warn("More than one automaton supplied. Processing only first one!");

  auto& aut = *auts.front();

  SCCInfo::uptr auti = get_scc_info(aut, true);
  auto deadscc = get_dead_sccs(aut, *auti);

  set<state_t> dead;
  for (auto s : aut.states()) {
    if (map_has_key(auti->scc, s)) {
        if (contains(deadscc,auti->scc.at(s)))
          dead.emplace(s);
      }
  }

  // detect accepting sinks (acc states with self loop for each sym)
  set<state_t> accsinks;
  auto tmp = get_accepting_sinks(aut);
  accsinks = set<state_t>(cbegin(tmp),cend(tmp));

  dot_header();

  auto states = aut.states();
  auto it = partition(begin(states), end(states), [&](state_t const& a){ return !contains(auti->unreachable, a); });
  states.erase(it, end(states));
  sort(begin(states), end(states), [&](state_t const& a, state_t const& b){
      return auti->scc.at(a) > auti->scc.at(b);
      });

  if (!auti->unreachable.empty()) {
    cout << "subgraph cluster_unreach {" << endl;
    cout << "label = <UNREACHABLE<BR/>(#st=" << auti->unreachable.size() << ")>;" << endl;
    cout << "style = \"filled\"; fillcolor = \"gray\";" << endl;
    for (auto s : auti->unreachable) {
      dot_state(s, aut, auti->unreachable, dead, accsinks, args->nodeinfo);
    }
    cout << "}" << endl;
  }

  int curscc = -1;
  for (auto s : states) {
    int sscc = auti->scc.at(s);
    if (sscc != curscc) {
      if (curscc != -1)
        cout << "}" << endl;
      curscc = sscc;
      cout << "subgraph cluster_scc" << sscc << " {" << endl;
      cout << "label = <<B>SCC " << sscc << "<BR/>(#st=" << auti->sccsz.at(sscc) << ")</B>>;" << endl;
      if (contains(auti->accepting, sscc)) {
        cout << "color = \"green\";" << endl;
      }
      if (contains(auti->rejecting, sscc)) {
        cout << "color = \"red\";" << endl;
      }
    }

    if (args->nodeinfo || auti->sccrep.at(curscc)==s)
      dot_state(s, aut, auti->unreachable, dead, accsinks, args->nodeinfo);
  }
  cout << "}" << endl;

  //node edges
  if (args->nodeinfo) {
    for (auto p : aut.states()) {
      for (auto q : aut.succ(p)) {
        string lbl;
        if (args->edgeinfo) {
          for (sym_t s : aut.outsyms(p)) {
            auto pssuc = aut.succ(p,s);
            if (find(begin(pssuc), end(pssuc), q) != end(pssuc)) {
              lbl += (lbl.empty() ? "" : ",");
              if (args->edgeinfo == 1)
                lbl += to_string(s);
              else if (args->edgeinfo == 2)
                lbl += sym_to_edgelabel(s, aut.get_aps(), true);
            }
          }
        }
        cout << p << " -> " << q << "[label=\""<< lbl << "\"]" << ";" << endl;
      }
    }
  }

  if (!args->nodeinfo) {
    for (auto it : auti->sccrep) {
      auto sccnum = it.first;
      auto sccrep = it.second;
      auto sucs = succ_sccs(aut, *auti, sccnum);
      for (auto sucscc : sucs) {
        auto sucrep = auti->sccrep.at(sucscc);
        cout << sccrep << " -> " << sucrep
          << "[ltail=\"cluster_scc"<<sccnum
          <<"\",lhead=\"cluster_scc"<<sucscc<<"\"];" << endl;
      }
    }
  }

  dot_footer();
}
