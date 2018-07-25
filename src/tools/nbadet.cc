#include <iostream>
#include <string>
#include <cassert>
using namespace std;

#include <spdlog/spdlog.h>
namespace spd = spdlog;
#include <args.hxx>

#include "metrics/bench.hh"
#include "metrics/memusage.h"

#include "aut.hh"
#include "io.hh"
#include "graph.hh"
#include "common/scc.hh"
#include "ps.hh"
#include "pa.hh"
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
  bool dsim;
  bool mindfa;

  int mergemode;
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
  args::Flag dsim(parser, "dsim", "Use direct simulation for preprocessing and optimization.",
      {'i', "dir-sim"});

  // postprocessing
  args::Flag mindfa(parser, "mindfa", "First minimize number of priorities, "
      "then minimize number of states using Hopcroft",
      {'m', "minimize-dfa"});

  // type of update for active ranks
  args::ValueFlag<int> mergemode(parser, "N", "Type of update "
      "(0=Muller/Schupp, 1=Safra, 2=Maximal merge)",
      {'u', "update-mode"});
  args::Flag puretrees(parser, "pure", "Accepting leaf normal form (acc. states in leaves only)",
      {'l', "pure-trees"});

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

  if (mergemode && args::get(mergemode) >= static_cast<int>(UpdateMode::num)) {
    spd::get("log")->error("Invalid update mode provided: {}", args::get(mergemode));
    exit(1);
  }

  //fill args
  Args args;
  if (input) {
    string filename = args::get(input);
    ifstream exists(filename);
    if (!exists) {
      spd::get("log")->error("File does not exist: {}", filename);
    }
    args.file = filename;
  }

  args.verbose = args::get(verbose);
  args.stats = stats;
  args.nooutput = nooutput;

  args.trim = trim;
  args.asinks = asinks;
  args.dsim = dsim;
  args.mindfa = mindfa;

  args.psets = psets;
  args.context = context;

  args.mergemode = args::get(mergemode);
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
  auto const scci = get_sccs(aut.states(), aut_suc);
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

  unordered_map<nba_bitset, int> numsets;
  map<pri_t, int> numpri;
  int mx=0;
  for (auto const st : pa.states()) {
    if (!pa.tag.hasi(st))
      continue;

    auto const psh = pa.tag.geti(st).powerset;
    numsets[psh]++;
    if (mx < numsets[psh])
      mx = numsets[psh];

    for (auto const sym : pa.state_outsyms(st)) {
      for (auto const es : pa.succ_edges(st, sym)) {
        numpri[es.second]++;
      }
    }
  }

  cerr << "#states: " << pa.num_states();
  cerr << " #psets: " << numsets.size();
  cerr << " maxseen: " << mx << endl;
  cerr << "prio:\t#: " << endl;
  for (auto const it : numpri) {
    cerr << it.first << "\t" << it.second << endl;
  }
}

PA process_nba(Args const &args, auto& aut, std::shared_ptr<spdlog::logger> log) {
    // -- preprocessing --

    //first trim (unmark trivial states that are accepting, remove useless+unreach SCCs)
    // just in case... usually input is already trim
    if (args.trim)
      ba_trim(aut, log);
    // aut->normalize(); //don't do this, otherwise relationship not clear anymore

    if (args.dsim) {
      auto const simret = ba_direct_sim(aut);
      aut = simret.first;
      //implications. TODO: use those to optimize det.:
      // - if two states in same reachable set and one dominates, keep dominating
      // for (auto const it : simret.second)
      //   cerr << it.first << " <= " << seq_to_str(it.second) << endl;
    }

    auto const dc = assemble_detconf(args, aut, log);
    if (args.verbose >= 2)
      cerr << dc << endl;

    //calculate 2^A and its sccs
    auto const pscon = bench(log,"powerset_construction",
                             WRAP(powerset_construction(aut, dc.aut_mat, dc.aut_asinks)));
    auto const pscon_scci = get_sccs(pscon.states(), aut_succ(pscon));
    log->info("#states in 2^A: {}, #SCCs in 2^A: {}", pscon.num_states(), pscon_scci.sccs.size());
    // print_aut(pscon);


    // -- end of preprocessing --

    //determinize (optionally using psets)
    PA pa;
    if (!args.psets)
      pa = bench(log, "determinize", WRAP(determinize(aut, dc)));
    else
      pa = bench(log, "determinize_with_psets", WRAP(determinize(aut, dc, pscon, pscon_scci)));

    if (args.stats) { //show stats before postprocessing
      print_stats(pa);
    }

    // -- begin postprocessing --
    if (args.mindfa) {
      auto optlog = args.verbose>1 ? log : nullptr;
      log->info("#priorities before: {}", pa.pris().size());
      log->info("#states before: {}", pa.states().size());

      pa.make_colored();
      bench(log, "minimize number of priorities", WRAP(minimize_priorities(pa, optlog)));
      log->info("#priorities after: {}", pa.pris().size());

      pa.make_complete();
      bench(log, "minimize number of states", WRAP(minimize_pa(pa, optlog)));
      log->info("#states after: {}", pa.num_states());
    }

    // -- end of postprocessing --

    //sanity checks
    assert(pa.get_name() == aut.get_name());
    assert(pa.get_aps() == aut.get_aps());
    assert(pa.is_deterministic());

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

    log->info("NBA name: \"{}\", #states: {}, #APs: {} #Syms: {}",
              aut.get_name(), aut.num_states(), aut.get_aps().size(), aut.num_syms());

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
