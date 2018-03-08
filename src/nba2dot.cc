#include <iostream>
#include <string>
using namespace std;

#include <spdlog/spdlog.h>
namespace spd = spdlog;

#include <args.hxx>

#include "io.hh"
#include "common/scc.hh"
#include "common/algo.hh"


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

void dot_state(state_t s, SWA<string> const& aut, vector<state_t> const& unreach,
       set<state_t> const& dead, vector<state_t> const& accsinks, int nodeinfo=0) {
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
  auto auts = nbautils::parse_hoa(args->file, log);

  if (auts.empty()) {
    log->error("Parsing NBAs from {} failed!", args->file.empty() ? "stdin" : args->file);
    exit(1);
  }

  if (auts.size() > 1)
    log->warn("More than one automaton supplied. Processing only first one!");

  auto& aut = auts.front();

  function<vector<state_t>(state_t)> sucs = [&aut](state_t v){ return aut->succ(v); };
  function<vector<state_t>(state_t,sym_t)> const xsucs = [&aut](state_t v,sym_t s){ return aut->succ(v,s); };
  function<vector<sym_t>(state_t)> const outsyms = [&aut](state_t v){ return aut->outsyms(v); };
  function<bool(state_t)> const ac = [&aut](state_t v){ return aut->has_accs(v); };

  auto states = aut->states();
  auto const unreach = unreachable_states(states, aut->get_init().front(), sucs);


  auto const accsinks = ba_get_accepting_sinks(states, aut->num_syms(), ac, outsyms, xsucs);

  auto const scci = get_sccs(states, sucs, const_true);
  auto const sucsccs = [&](unsigned num){ return succ_sccs(scci, num, sucs); };

  auto const trivial = trivial_sccs(scci, sucs);
  auto const bascl = ba_classify_sccs(scci, ac);

  auto const dead = ba_get_dead_sccs(scci.sccs.size(), bascl.rejecting, trivial, sucsccs);


  dot_header();

  auto it = partition(begin(states), end(states), [&](state_t const& a){ return !contains(unreach, a); });
  states.erase(it, end(states));
  sort(begin(states), end(states),
      [&](state_t const& a, state_t const& b){
        return scci.scc_of.at(a) > scci.scc_of.at(b);
      });

  if (!unreach.empty()) {
    cout << "subgraph cluster_unreach {" << endl;
    cout << "label = <UNREACHABLE<BR/>(#st=" << unreach.size() << ")>;" << endl;
    cout << "style = \"filled\"; fillcolor = \"gray\";" << endl;
    for (auto s :unreach) {
      dot_state(s, *aut, unreach, dead, accsinks, args->nodeinfo);
    }
    cout << "}" << endl;
  }

  int curscc = -1;
  for (auto s : states) {
    int sscc = scci.scc_of.at(s);
    if (sscc != curscc) {
      if (curscc != -1)
        cout << "}" << endl;
      curscc = sscc;
      cout << "subgraph cluster_scc" << sscc << " {" << endl;
      cout << "label = <<B>SCC " << sscc << "<BR/>(#st=" << scci.sccs.at(sscc).size() << ")</B>>;" << endl;
      if (contains(bascl.accepting, sscc)) {
        cout << "color = \"green\";" << endl;
      }
      if (contains(bascl.rejecting, sscc)) {
        cout << "color = \"red\";" << endl;
      }
    }

    if (args->nodeinfo || scci.sccs.at(curscc).front()==s)
      dot_state(s, *aut, unreach, dead, accsinks, args->nodeinfo);
  }
  cout << "}" << endl;

  //node edges
  if (args->nodeinfo) {
    for (auto p : aut->states()) {
      for (auto q : aut->succ(p)) {
        string lbl;
        if (args->edgeinfo) {
          for (sym_t s : aut->outsyms(p)) {
            auto pssuc = aut->succ(p,s);
            if (find(begin(pssuc), end(pssuc), q) != end(pssuc)) {
              lbl += (lbl.empty() ? "" : ",");
              if (args->edgeinfo == 1)
                lbl += to_string(s);
              else if (args->edgeinfo == 2)
                lbl += sym_to_edgelabel(s, aut->get_aps(), true);
            }
          }
        }
        cout << p << " -> " << q << "[label=\""<< lbl << "\"]" << ";" << endl;
      }
    }
  }

  if (!args->nodeinfo) {
    unsigned i=0;
    for (auto it : scci.sccs) {
      auto sccrep = it.front();
      auto sucsccs = succ_sccs(scci, i, sucs);
      for (auto sucscc : sucsccs) {
        auto sucrep = scci.sccs.at(sucscc).front();
        cout << sccrep << " -> " << sucrep
          << "[ltail=\"cluster_scc"<<i
          <<"\",lhead=\"cluster_scc"<<sucscc<<"\"];" << endl;
      }

      ++i;
    }
  }

  dot_footer();
}
