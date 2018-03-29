#include <fstream>
#include <iostream>
#include <string>
#include <cassert>
using namespace std;

#include <spdlog/spdlog.h>
namespace spd = spdlog;

#include <args.hxx>

#include "io.hh"
#include "common/algo.hh"
#include "ba.hh"
#include "ps.hh"
#include "pa.hh"
#include "det.hh"

#include "dev/bench.hh"
#include "dev/memusage.h"

using namespace nbautils;

struct Args {
  using uptr = std::unique_ptr<Args>;

  string file;
  int verbose;

  bool trim;
  int split;

  bool detaccsinks;

  LevelUpdateMode lvupdate;
  bool seprej;
  bool sepacc;
  bool cyclicbrk;
  bool pure;
  bool context;

  bool topo;
  bool optguarded;

  bool minpri;
  bool mindfa;

  bool nooutput;
};

Args::uptr parse_args(int argc, char *argv[]) {
  args::ArgumentParser parser("nbadet - determinize nondeterministic Büchi automata", "");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});

  // exactly one input automaton
  args::Positional<string> input(parser, "INPUTFILE", "file containing the NBA (if none given, uses stdin)");

  // logging level -v, -vv, etc.
  args::CounterFlag verbose(parser, "verbose", "Show verbose information", {'v', "verbose"});

  args::Flag nooutput(parser, "nooutput", "Do not print resulting automaton", {'x', "no-output"});

  // preprocessing on NBA (known, simple stuff)
  args::Flag trim(parser, "trim", "Remove dead states from NBA", {'d', "trim"});
  args::Flag detaccsinks(parser, "detaccsinks", "Detect accepting sinks", {'s', "detect-acc-sinks"});

  // additional calculations on NBA to optimize construction
  args::Flag context(parser, "context", "Calculate context for separation refinement", {'c', "use-context"});

  // enabled optimizations for Safra/Level update
  args::Flag seprej(parser, "seprej", "Separate states in non-accepting SCCs", {'n', "separate-rej"});
  args::Flag sepacc(parser, "sepacc", "Separate states in accepting SCCs", {'a', "separate-acc"});
  args::Flag cyclicbrk(parser, "cyclicbrk", "Separate states in accepting SCCs, cycle through SCCs", {'b', "cyclic-breakpoint"});
  args::Flag pure(parser, "pure", "Accepting leaf normal form (acc. states in leaves only)", {'l', "pure"});

  // type of update
  args::ValueFlag<int> update(parser, "level-update", "Type of update", {'u', "level-update"});

  // construction methods

  // used to weed out redundant SCCs in det. automaton
  args::Flag topo(parser, "topo", "Use powerset SCCs to guide determinization", {'t', "topological"});
  args::Flag optguarded(parser, "optguarded", "Optimize states guarded by higher-priority states", {'g', "opt-guarded"});

  // iterated product construction based determinization
  args::CounterFlag split(parser, "split", "Determinize all NBA SCCs separately, then combine", {'i', "split"});

  // postprocessing
  args::Flag minpri(parser, "minpri", "Minimize number of priorities", {'p', "minimize-priorities"});
  args::Flag mindfa(parser, "mindfa", "Minimize number of states using Hopcroft", {'m', "minimize-dfa"});

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

  if (update && args::get(update) >= static_cast<int>(LevelUpdateMode::num)) {
    spd::get("log")->error("Invalid update mode provided: {}", args::get(update));
    exit(1);
  }

  auto args = make_unique<Args>(Args());
  if (input) args->file = args::get(input);
  args->verbose = args::get(verbose);
  args->trim = trim;
  args->split = args::get(split);

  args->lvupdate = static_cast<LevelUpdateMode>(args::get(update));
  args->detaccsinks = detaccsinks;
  args->sepacc = sepacc || cyclicbrk;
  args->cyclicbrk = cyclicbrk;
  args->seprej = seprej;
  args->context = context;
  args->pure = pure;

  args->topo = topo;
  args->optguarded = optguarded;
  args->minpri = minpri;
  args->mindfa = mindfa;

  args->nooutput = nooutput;

  return move(args);
}

LevelConfig::uptr levelconfig_from_args(Args const &args) {
  auto lc = make_unique<LevelConfig>(LevelConfig());
  lc->update = args.lvupdate;
  lc->optguarded = args.optguarded;
  lc->pure = args.pure;
  lc->sep_rej = args.seprej;
  lc->sep_acc = args.sepacc;
  lc->sep_acc_cyc = args.cyclicbrk;
  return move(lc);
}

auto determinize_nba(Args const &args, SWA<string>& aut, std::shared_ptr<spdlog::logger> log) {
    auto aut_st = aut.states();
    succ_fun<state_t> const aut_sucs = [&aut](state_t v){ return aut.succ(v); };
    function<bool(state_t)> const aut_acc = [&aut](state_t v){ return aut.has_accs(v); };
    succ_sym_fun<state_t, sym_t> const aut_xsucs = [&aut](state_t v,sym_t s){ return aut.succ(v,s); };
    outsym_fun<state_t,sym_t> const aut_osyms = [&aut](state_t v){ return aut.outsyms(v); };

    // calculate SCCs of NBA if needed
    SCCDat<state_t>::uptr aut_scc = nullptr;
    BaSccClassification::uptr aut_cl = nullptr;
    if (args.sepacc || args.seprej) {
      aut_scc = get_sccs(aut_st, aut_sucs);
      aut_cl = ba_classify_sccs(*aut_scc, aut_acc);
      log->info("number of SCCs in A: {}", aut_scc->sccs.size());
    }

    // detect accepting sinks (acc states with self loop for each sym)
    vector<small_state_t> accsinks;
    if (args.detaccsinks) {
      accsinks = to_small_state_t(get_accepting_sinks(aut_st, aut.num_syms(), aut_acc, aut_osyms, aut_xsucs));
      log->info("found {} accepting sinks", accsinks.size());
    }

    // calculate 2^A to guide and optimize determinization
    PS::uptr ps = nullptr;
    SCCDat<state_t>::uptr psi = nullptr;
    if (args.topo) {
      ps = bench(log, "powerset_construction", WRAP(powerset_construction(aut, accsinks)));
      succ_fun<state_t> pssucs = [&ps](state_t v){ return ps->succ(v); };
      psi = bench(log, "get_scc_info",          WRAP(get_sccs(ps->states(), pssucs)));
      log->info("number of states in 2^A: {}", ps->num_states());
      log->info("number of SCCs in 2^A: {}", psi->sccs.size());
      assert(is_deterministic(*ps));
    }

    // calculate A x 2^A as context information
    PP::uptr ctx = nullptr;
    SCCDat<state_t>::uptr ctx_scc = nullptr;
    BaSccClassification::uptr ctx_cl = nullptr;
    if (args.context) {
      ctx = bench(log, "powerset_product", WRAP(powerset_product(aut)));
      succ_fun<state_t> ctxsucs = [&ctx](state_t v){ return ctx->succ(v); };
      function<bool(state_t)> ctx_acc = [&ctx](state_t v){ return ctx->has_accs(v); };
      auto ctxst = ctx->states();

      ctx_scc = bench(log, "get_scc_info", WRAP(get_sccs(ctxst, ctxsucs)));
      ctx_cl = bench(log, "classify_sccs", WRAP(ba_classify_sccs(*ctx_scc, ctx_acc)));
      log->info("number of states in Ax2^A: {}", ctx->num_states());
      log->info("number of SCCs in Ax2^A: {}", ctx_scc->sccs.size());
    }

    // configure level update:
    auto lc = levelconfig_from_args(args);
    //set the accepting sink states
    lc->accsinks = accsinks;
    // set reference to underlying NBA
    lc->aut = &aut;
    lc->aut_scc = aut_scc.get();
    lc->aut_cl = aut_cl.get();

    // set reference to context NBA
    lc->ctx = ctx.get();
    lc->ctx_scc = ctx_scc.get();
    lc->ctx_cl = ctx_cl.get();

    PA::uptr pa;
    if (!args.topo)
      pa = bench(log, "determinize", WRAP(determinize(*lc)));
    else
      pa = bench(log, "determinize_topo", WRAP(determinize(*lc, *ps, *psi)));

    make_complete(*pa);

    log->info("number of states in parity automaton: {}", pa->num_states());

    //sanity checks
    assert(pa->get_name() == aut.get_name());
    assert(pa->get_aps() == aut.get_aps());
    assert(pa->get_init().size() == 1);
    assert(is_deterministic(*pa));
    assert(is_colored(*pa));
    assert(is_complete(*pa));

    if (args.minpri) {
      bench(log, "minimize number of priorities", WRAP(minimize_priorities(*pa)));
      log->info("resulting number of priorities: {}", pa->get_accsets().size());
    }

    if (args.mindfa) {
      bench(log, "minimize number of states", WRAP(minimize_pa(*pa)));
      log->info("number of states after minimization: {}", pa->num_states());
    }

    return move(pa);
}

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

int main(int argc, char *argv[]) {
  // initialize stuff (args + logging)

  // auto console = spd::stdout_color_mt("log");
  auto log = spd::stderr_logger_mt("log");
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

  auto totalstarttime = get_time();
  for (auto &aut : auts) {
    if (aut->acond != Acceptance::BUCHI) {
      log->warn("skipping non-Büchi automaton with name \"{}\"", aut->get_name());
      continue;
    }

    log->info("processing NBA with name \"{}\"", aut->get_name());
    log->info("number of states in A: {}", aut->num_states());
    log->info("number of APs in A: {}", aut->get_aps().size());

    // sanity check size of the input
    if (aut->num_states() > max_nba_states) {
      log->error("NBA is way too large, I refuse.");
      exit(1);
    }
    if (aut->get_aps().size() > max_nba_syms) {
      log->error("Alphabet is way too large, I refuse.");
      exit(1);
    }

    // now we start measuring time
    //---------------------------
    auto starttime = get_time();

    if (args->trim) { // just in case... usually input is already trim
      auto const numtrimmed = trim_ba(*aut);
      // aut->normalize();
      log->info("removed {} useless states", numtrimmed);
    }

    if (args->split==0) { //determinize in one piece
      auto pa = determinize_nba(*args, *aut, log);

      log->info("completed automaton in {:.3f} seconds", get_secs_since(starttime));
      if (!args->nooutput)
        print_hoa(*pa);

    } else {
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

      log->info("completed automaton in {:.3f} seconds", get_secs_since(starttime));
      if (!args->nooutput)
        print_hoa(*res);
    }

  }

  log->info("total time: {:.3f} seconds", get_secs_since(totalstarttime));
  log->info("total used memory: {:.3f} MB", (double)getPeakRSS() / (1024 * 1024));
}
