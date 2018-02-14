#include <fstream>
#include <iostream>
#include <string>
using namespace std;

#include <spdlog/spdlog.h>
namespace spd = spdlog;

#include <args.hxx>

#include "det.hh"
#include "io.hh"
#include "ps.hh"
#include "scc.hh"
#include "types.hh"

#include "bench.hh"
#include "memusage.h"

#include "debug.hh"

using namespace nbautils;

struct Args {
  typedef std::unique_ptr<Args> uptr;

  string file;
  int verbose;
  bool trim;

  LevelUpdateMode lvupdate;
  bool seprej;
  bool sepacc;
  bool context;

  bool topo;
};

Args::uptr parse_args(int argc, char *argv[]) {
  args::ArgumentParser parser("nbadet - determinize nondeterministic Büchi automata", "");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});

  // exactly one input automaton
  args::Positional<string> input(parser, "INPUTFILE",
                                 "file containing the NBA (if none given, uses stdin)");

  // logging level -v, -vv, etc.
  args::CounterFlag verbose(parser, "verbose", "Show verbose information",
                            {'v', "verbose"});

  // preprocessing on NBA
  args::Flag trim(parser, "trim", "Remove dead states from NBA", {'d', "trim"});

  // calculations on NBA to optimize construction
  // args::Flag sim(parser, "sim", "Calculate simulation relation on NBA to inform
  // determinization", {'b', "use-sim"});
  args::Flag context(parser, "context", "Calculate context for separation refinement",
                     {'c', "use-context"});

  // enabled optimizations
  args::Flag seprej(parser, "seprej", "Separate states in non-accepting SCCs",
                    {'n', "separate-rej"});
  args::Flag sepacc(parser, "sepacc", "Separate states in accepting SCCs",
                    {'a', "separate-acc"});

  // type of update
  args::ValueFlag<int> update(parser, "level-update", "Type of update",
                              {'u', "level-update"});

  // construction methods

  // used to weed out redundant SCCs in det. automaton
  args::Flag topo(parser, "topo", "Use powerset SCCs to guide determinization",
                  {'t', "topological"});

  // iterated product construction based determinization
  // args::Flag split(parser, "split", "Determinize all NBA SCCs separately, then
  // combine", {'s', "split"}); args::Flag rabin(parser, "rabin", "Return a smaller Rabin
  // automaton", {'r', "split"});

  // postprocessing
  // args::Flag minpri(parser, "minpri", "Minimize number of priorities", {'m',
  // "minimize-priorities"});

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

  if (context && !(seprej || sepacc)) {
    spd::get("log")->error("-c without at least one of -a or -n is useless!");
    exit(1);
  }

  if (update && args::get(update) >= LevelUpdateMode::num) {
    spd::get("log")->error("Invalid update mode provided: {}", args::get(update));
    exit(1);
  }

  auto args = make_unique<Args>(Args());
  if (input) args->file = args::get(input);
  args->verbose = args::get(verbose);
  args->trim = trim;

  args->lvupdate = static_cast<LevelUpdateMode>(args::get(update));
  args->sepacc = sepacc;
  args->seprej = seprej;
  args->context = context;

  args->topo = topo;

  return move(args);
}

LevelConfig::uptr levelconfig_from_args(Args const &args) {
  auto lc = make_unique<LevelConfig>(LevelConfig());
  lc->update = args.lvupdate;
  lc->sep_rej = args.seprej;
  lc->sep_acc = args.sepacc;
  return move(lc);
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
  BA::uptr aut = nbautils::parse_ba(args->file);
  if (!aut) {
    log->error("failed parsing NBA from {}!", args->file);
    exit(1);
  }
  log->info("number of states in A: {}", aut->adj.size());

  // sanity check size of the input
  if (aut->adj.size() > max_nba_states) {
    spd::get("log")->error("NBA is way too large, I refuse.");
    exit(1);
  }
  if (aut->meta.aps.size() > max_nba_syms) {
    spd::get("log")->error("Alphabet is way too large, I refuse.");
    exit(1);
  }

  // now we start measuring time
  //---------------------------
  auto starttime = get_time();

  // calculate SCCs of NBA if needed
  SCCInfo::uptr auti = nullptr;
  if (args->trim || args->sepacc || args->seprej) {
    auti = get_scc_info(*aut, true);
    log->info("number of SCCs in A: {}", auti->sccrep.size());
  }

  if (args->trim) {
    log->info("trimming original NBA...");
    auto numtrimmed = trim_ba(*aut, *auti);
    log->info("removed {} useless states", numtrimmed);
    log->info("number of states in trimmed A: {}", aut->adj.size());
    log->info("number of SCCs in trimmed A: {}", auti->sccrep.size());
  }

  // calculate 2^A to guide and optimize determinization
  BAPS::uptr ps = nullptr;
  SCCInfo::uptr psi = nullptr;
  if (args->topo) {
    ps = powerset_construction(*aut);
    psi = get_scc_info(*ps, false);
    log->info("number of states in 2^A: {}", ps->adj.size());
    log->info("number of SCCs in 2^A: {}", psi->sccrep.size());
    // printSCCI(*ctxi);
    // printPS(*ctx, *ctxi, true);
  }

  // calculate A x 2^A as context information
  BAPP::uptr ctx = nullptr;
  SCCInfo::uptr ctxi = nullptr;
  if (args->context) {
    ctx = powerset_product(*aut);
    ctxi = get_scc_info(*ctx, true);
    log->info("number of states in Ax2^A: {}", ctx->adj.size());
    log->info("number of SCCs in Ax2^A: {}", ctxi->sccrep.size());
  }

  // configure level update:
  auto lc = levelconfig_from_args(*args);
  // set reference to underlying NBA
  lc->aut = aut.get();
  lc->auti = auti.get();
  // set reference to context NBA
  lc->ctx = ctx.get();
  lc->ctxi = ctxi.get();

  PA::uptr pa = determinize(*lc);
  log->info("number of states in resulting automaton: {}", pa->num_states());

  // TODO: apply postprocessing
  // TODO: output in HOA format

  //---------------------------
  log->info("completed after {:.3f} seconds", get_secs_since(starttime));
  log->info("used {:.3f} MB of memory", (double)getPeakRSS() / (1024 * 1024));
}
