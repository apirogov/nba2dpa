#include <iostream>
#include <string>
#include <cassert>
using namespace std;

#include <spdlog/spdlog.h>
namespace spd = spdlog;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcatch-value"
#include <args.hxx>
#pragma GCC diagnostic pop

#include "metrics/bench.hh"
#include "metrics/memusage.h"

#include "aut.hh"
#include "io.hh"
#include "graph.hh"
#include "common/scc.hh"
#include "ps.hh"
#include "preproc.hh"
#include "detstate.hh"
#include "det.hh"

using namespace nbautils;

struct Args {
  string file;

  int verbose;
  bool stats;
  bool nooutput;

  bool trim;
  bool asinks;
  bool mindfa;

  int mergemode;
  bool weaksat;
  bool puretrees;

  bool psets;
  bool context;

  bool seprej;
  bool sepacc;
  bool cyclicbrk;
  bool sepmix;
  bool optdet;

  //not included: incremental w. locks, using external union of SCC det., simulation stuff
};

Args parse_args(int argc, char *argv[]) {
  args::ArgumentParser parser("nbadet - determinize nondeterministic Büchi automata", "");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});

  args::Positional<string> input(parser, "INPUTFILE",
      "file containing the NBA(s) (if none given, uses <stdin>)");

  // introspection and debugging
  args::CounterFlag verbose(parser, "verbose", "Show verbose information",
      {'v', "verbose"}); // logging level -v, -vv, etc.
  args::CounterFlag stats(parser, "stats", "Output stats about structure",
      {'s', "output-stats"});
  args::Flag nooutput(parser, "nooutput", "Do not print resulting automaton",
      {'x', "no-output"});

  // preprocessing on NBA (known, simple stuff)
  args::Flag trim(parser, "trim", "Kill dead states from NBA.",
      {'k', "trim"});
  args::Flag asinks(parser, "asinks", "Detect and use accepting (pseudo)sinks.",
      {'j', "acc-sinks"});

  // postprocessing
  args::Flag mindfa(parser, "mindfa", "First minimize number of priorities, "
      "then minimize number of states using Hopcroft",
      {'m', "minimize-dfa"});

  // type of update for active ranks
  args::ValueFlag<int> mergemode(parser, "update-mode", "Type of update "
      "(Muller/Schupp, Safra, Maximal merge)",
      {'u', "update-mode"});
  args::Flag puretrees(parser, "pure", "Accepting leaf normal form (acc. states in leaves only)",
      {'l', "pure-trees"});
  args::Flag weaksat(parser, "weak-saturation", "Allow weak saturation",
      {'w', "weak-saturation"});

  // used to weed out redundant SCCs in det. automaton
  args::Flag psets(parser, "powersets", "Use powerset SCCs to guide determinization",
      {'t', "use-powersets"});

  // additional calculations on NBA to optimize construction
  args::Flag context(parser, "context", "Calculate context for separation refinement",
      {'c', "use-context"});

  // enabled optimizations for Safra/Level update
  args::Flag seprej(parser, "seprej", "Separate states in non-accepting SCCs",
      {'n', "sep-rej"});
  args::Flag sepacc(parser, "sepacc", "Separate states in accepting SCCs",
      {'a', "sep-acc"});
  args::Flag cyclicbrk(parser, "cyclicbrk", "Separate states in accepting SCCs, cycle through SCCs",
      {'b', "cyclic-breakpoint"});
  args::Flag sepmix(parser, "sepmix", "Separate states in different SCCs",
      {'e', "sep-mix"});
  args::Flag optdet(parser, "optdet", "Optimize deterministic SCCs by not expanding trees",
      {'d', "opt-det"});


  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help&) {
    std::cout << parser;
    exit(0);
  } catch (args::ParseError& e) {
    cerr << e.what() << endl << parser;
    exit(1);
  } catch (args::ValidationError& e) {
    cerr << e.what() << endl << parser;
    exit(1);
  }

  if (context && !(seprej || sepacc)) {
    spd::get("log")->error("-c without at least one of -a or -n is useless!");
    exit(1);
  }

  if (cyclicbrk && !sepacc) {
    spd::get("log")->error("-b without -a is useless!");
    exit(1);
  }

  // if (cyclicbrk && context) {
  //   spd::get("log")->error("-b and -c don't work together (yet)!");
  //   exit(1);
  // }

  if (optdet && !sepmix) {
    spd::get("log")->error("-d without -e does not work!");
    exit(1);
  }

  if (mergemode && args::get(mergemode) >= static_cast<int>(UpdateMode::num)) {
    spd::get("log")->error("Invalid update mode provided: {}", args::get(mergemode));
    exit(1);
  }

  //fill args
  Args args;
  if (input)
    args.file = args::get(input);

  args.verbose = args::get(verbose);
  args.stats = stats;
  args.nooutput = nooutput;

  args.trim = trim;
  args.asinks = asinks;
  args.mindfa = mindfa;

  args.psets = psets;
  args.context = context;

  args.mergemode = args::get(mergemode);
  args.weaksat = weaksat;
  args.puretrees = puretrees;

  args.seprej = seprej;
  args.sepacc = sepacc;
  args.cyclicbrk = cyclicbrk;
  args.sepmix = sepmix;
  args.optdet = optdet;

  return args;
}

//fill DetConf flags from args
DetConf detconf_from_args(Args const& args) {
  DetConf dc;

  dc.debug = args.verbose > 2;

  dc.update = static_cast<UpdateMode>(args.mergemode);
  dc.weaksat = args.weaksat;
  dc.puretrees = args.puretrees;

  dc.sep_rej = args.seprej;
  dc.sep_acc = args.sepacc;
  dc.sep_acc_cyc = args.cyclicbrk;
  dc.sep_mix = args.sepmix;
  dc.opt_det = args.optdet;

  return dc;
}

//given NBA and detconf without sets, return the corresponding configured sets
DetConfSets get_detconfsets(auto const& aut, DetConf const& dc,
                    shared_ptr<spdlog::logger> log = nullptr) {
  auto const aut_suc = aut_succ(aut);
  auto const scci = get_sccs(aut.states() | ranges::to_vector, aut_suc);
  auto const sccDet = ba_scc_classify_det(aut, scci);
  auto const sccAcc = ba_scc_classify_acc(aut, scci);
  assert(sccAcc.size() == scci.sccs.size());
  int a=0;
  int n=0;
  int m=0;
  for (auto const& it : sccAcc) {
    if (it.second == -1)     n++;
    else if (it.second == 0) m++;
    else if (it.second == 1) a++;
  }
  if (log)
    log->info("SCCs: {} total = {} A + {} N + {} M, of which {} D",
              scci.sccs.size(), a, n, m, sccDet.size());

  return calc_detconfsets(dc, scci, sccAcc, sccDet);
}

//given args and NBA, prepare corresponding determinization config structure
DetConf assemble_detconf(Args const& args, auto const& aut,
                    shared_ptr<spdlog::logger> log = nullptr) {
  auto dc = detconf_from_args(args);
  dc.aut_states = to_bitset<nba_bitset>(aut.states());
  dc.aut_acc    = to_bitset<nba_bitset>( aut.states() | ranges::view::remove_if(
                    [&](state_t s){ return !aut.state_buchi_accepting(s); }));
  //get adj matrix for accelerated powerset calculation
  dc.aut_mat = get_adjmat(aut);

  //get accepting sinks
  dc.aut_asinks = 0;
  if (args.asinks)
    dc.aut_asinks = to_bitset<nba_bitset>(ba_get_acc_sinks(aut, log));

  //calculate 2^AxA context structure and its sccs
  if (args.context)
    dc.ctx = get_context(aut, dc.aut_mat, dc.aut_asinks, log);

  //precompute lots of sets used in construction
  dc.sets = get_detconfsets(aut, dc, log);
  return dc;
}

//output stats about active priorities, numstates, SCCs, different trees overall and per powerset, etc.
void print_stats(auto const& pa) {
  cerr << "#states: " << pa.num_states() << endl;
  cerr << "TODO: stats" << endl;
  /*
    if (args.stats) {
      unordered_map<bitset<256>, int> numsets;
      int mx=0;
      for (auto const st : pa.states()) {
        if (!pa.tag.hasi(st))
          continue;

        Aut<Level> const psh = pa.tag.geti(st).powerset;
        numsets[psh]++;
        if (mx < numsets[psh])
          mx = numsets[psh];
      }
      cerr << pa->num_states() << " states, " << numsets.size() << " psets, each at most " << mx << " times" << endl;
    }
    */
}

PA process_nba(Args const &args, auto& aut, std::shared_ptr<spdlog::logger> log) {
    // -- preprocessing --

    //first trim (unmark trivial states that are accepting, remove useless+unreach SCCs)
    // just in case... usually input is already trim
    if (args.trim)
      ba_trim(aut, log);
    // aut->normalize(); //don't do this, otherwise relationship not clear anymore

    auto const dc = assemble_detconf(args, aut, log);
    if (dc.debug)
      cerr << dc << endl;

    //calculate 2^A and its sccs
    auto const pscon = bench(log,"powerset_construction",
                             WRAP(powerset_construction(aut, dc.aut_mat, dc.aut_asinks)));
    auto const pscon_scci = get_sccs(pscon.states() | ranges::to_vector, aut_succ(pscon));
    log->info("#states in 2^A: {}, #SCCs in 2^A: {}", pscon.num_states(), pscon_scci.sccs.size());
    // print_aut(pscon);

    // -- end of preprocessing --

    //TODO: determinize (optionally using topo)
    auto pa = determinize(aut, dc);
    /*
    PA::uptr pa;
    if (!args.topo)
      pa = bench(log, "determinize", WRAP(determinize(*lc)));
    else
      pa = bench(log, "determinize_topo", WRAP(determinize(*lc, *ps, *psi)));
    */

    if (args.stats) { //show stats before postprocessing
      print_stats(pa);
    }

    // -- begin postprocessing --

    pa.make_complete();
    pa.make_colored();

    //TODO: minimize priorities
    //TODO: minimize states

    /*
    if (args.minpri) {
      bench(log, "minimize number of priorities", WRAP(minimize_priorities(*pa)));
      log->info("resulting number of priorities: {}", pa->get_accsets().size());
      bench(log, "minimize number of states", WRAP(minimize_pa(*pa)));
      log->info("number of states after minimization: {}", pa->num_states());
    }
    */

    // -- end of postprocessing --

    //sanity checks
    assert(pa.get_name() == aut.get_name());
    assert(pa.get_aps() == aut.get_aps());
    // assert(pa.is_complete());
    // assert(pa.is_colored());
    // assert(pa.is_deterministic());

    if (args.stats) { //show stats after postprocessing
      print_stats(pa);
    }

    return pa;
}

int main(int argc, char *argv[]) {
  // initialize stuff (args + logging):

  // auto console = spd::stdout_color_mt("log");
  auto const log = spd::stderr_logger_mt("log");
  spd::set_pattern("[%Y-%m-%d %H:%M:%S %z] [%l] %v");

  auto const args = parse_args(argc, argv);
  if (!args.verbose)
    spd::set_level(spd::level::warn);
  else if (args.verbose == 1)
    spd::set_level(spd::level::info);
  else
    spd::set_level(spd::level::debug);

  auto const totalstarttime = get_time();

  // now parse input automata:
  auto auts = nbautils::AutStream<Aut<string>>(args.file, log);
  while (auts.has_next()) {
    auto aut = auts.parse_next();

    log->info("NBA name: \"{}\", #states: {}, #APs: {}",
              aut.get_name(), aut.num_states(), aut.get_aps().size());

    // sanity of the input
    if (!aut.is_buchi()) {
      log->error("This is not an NBA!");
      exit(1);
    }
    if (aut.num_states() > max_nba_states) {
      log->error("NBA is way too large, I refuse.");
      exit(1);
    }
    if (aut.get_aps().size() > max_nba_syms) {
      log->error("Alphabet is way too large, I refuse.");
      exit(1);
    }

    // NBA -> DPA
    PA const pa = bench(log,"process_nba",
                          WRAP(process_nba(args, aut, log)));

    if (!args.nooutput)
      print_aut(pa);
  }

  log->info("total time: {:.3f} seconds", get_secs_since(totalstarttime));
  log->info("total used memory: {:.3f} MB", (double)getPeakRSS() / (1024 * 1024));
}

//old code for split NBA det.
  /*
vector<SWA<string>::uptr> split_nba(SWA<string> const& aut, bool each_acc_separated) {
  vector<SWA<string>::uptr> ret;

  auto const aut_st = aut.states();
  auto const aut_sucs = swa_succ(aut);
  auto const aut_acc = swa_ba_acc(aut);

  auto const aut_scc = get_sccs(aut_st, aut_sucs);
  auto const aut_cl = ba_classify_sccs(*aut_scc, aut_acc);
  unsigned i=0;
  for (auto const& scc : aut_scc->sccs) {
    if (contains(aut_cl->rejecting, i)) {
      ++i;
      continue;
    }

    //get accepting states of scc
    vector<state_t> sccacc = scc;
    auto it = partition(begin(sccacc), end(sccacc), aut_acc);
    sccacc.erase(it, end(sccacc));
    sort(begin(sccacc), end(sccacc));

    // cerr << seq_to_str(sccacc) << endl;

    //unmark all other accepting states in copy
    auto cur = make_unique<SWA<string>>(aut);

    // cerr << "copied" << endl;

    for (auto const s : aut_st) {
      // if (curaut->has_accs(s))
        if (!contains(sccacc, s))
          cur->set_accs(s, {});
    }

    // cerr << "unmarked" << endl;

    //throw everything else away
    trim_ba(*cur);
    // curaut->normalize();

    if (!each_acc_separated) {
      ret.push_back(move(cur));
    } else {
      for (int i=0; i<(int)sccacc.size(); ++i) {
        auto subcur = make_unique<SWA<string>>(*cur);
        for (int j=0; j<(int)sccacc.size(); ++j) {
          if (i!=j)
            subcur->set_accs(sccacc[j], {});
        }
        // print_hoa(*subcur);
        ret.push_back(move(subcur));
      }
    }

    // cerr << "trimmed" << endl;

    ++i;
  }

  return ret;
}

  for (auto &aut : auts) {
      log->info("splitting NBA");
      auto subauts = split_nba(*aut, args->split > 1);
      log->info("number of subNBAs: {}", subauts.size());
      SWA<PAProdState,naive_unordered_bimap>::uptr res = nullptr;

      for (auto const& subaut : subauts) {
        log->info("number of states in subA: {}", subaut->num_states());
        auto subpa = determinize_nba(*args, *subaut, log);
        // print_hoa(*subpa);

        log->info("calculating union...");
        if (res == nullptr) {
          auto tmp = empty_pa<bool,naive_unordered_bimap>(subpa->get_aps());
          res = pa_union(*tmp, *subpa);
        } else {
          auto tmp = pa_union(*res, *subpa);
          res = move(tmp);
        }

        log->info("completed union. current union size: {}", res->num_states());
        log->info("optimizing union...");
        minimize_priorities(*res);
        minimize_pa(*res);
        // auto equiv = pa_equiv_states(*res);
        // res->quotient(equiv);
        // res->normalize();
        log->info("optimized union. current union size: {}", res->num_states());

        //check language containment to abort earlier
      }
  }
    */
