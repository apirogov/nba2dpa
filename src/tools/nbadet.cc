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

  if (optdet && !sepmix) {
    spd::get("log")->error("-d without -e does not work!");
    exit(1);
  }

  /*
  if (update && args::get(update) >= static_cast<int>(LevelUpdateMode::num)) {
    spd::get("log")->error("Invalid update mode provided: {}", args::get(update));
    exit(1);
  }
  */

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

/*
auto determinize_nba(Args const &args, SWA<string>& aut, std::shared_ptr<spdlog::logger> log) {
    auto aut_st = aut.states();
    succ_fun<state_t> const aut_sucs = [&aut](state_t v){ return aut.succ(v); };
    function<bool(state_t)> const aut_acc = [&aut](state_t v){ return aut.has_accs(v); };
    succ_sym_fun<state_t, sym_t> const aut_xsucs = [&aut](state_t v,sym_t s){ return aut.succ(v,s); };
    outsym_fun<state_t,sym_t> const aut_osyms = [&aut](state_t v){ return aut.outsyms(v); };

    // calculate SCCs of NBA if needed
    SCCDat<state_t>::uptr aut_scc = nullptr;
    BaSccClassification::uptr aut_cl = nullptr;
    if (args.sepacc || args.seprej || args.stats) {
      aut_scc = get_sccs(aut_st, aut_sucs);
      aut_cl = ba_classify_sccs(*aut_scc, aut_acc);
      log->info("number of SCCs in A: {}", aut_scc->sccs.size());
    }

    // detect accepting sinks (acc states with self loop for each sym)
    vector<small_state_t> accsinks;
    if (args.detaccsinks || args.stats) {
      accsinks = to_small_state_t(get_accepting_sinks(aut_st, aut.num_syms(), aut_acc, aut_osyms, aut_xsucs));
      log->info("found {} accepting sinks", accsinks.size());
    }

    if (args.stats) {
      int msccs = aut_scc->sccs.size() - aut_cl->accepting.size() - aut_cl->rejecting.size();
      cerr << "ASCCs: " << aut_cl->accepting.size() << "\tNSCCs: " << aut_cl->rejecting.size() << "\tMSCCs: " << msccs << "\tAsinks: " << accsinks.size() << endl;
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
*/
auto process_nba(Args const &args, auto& aut, std::shared_ptr<spdlog::logger> log) {
    // -- preprocessing --

    //first trim (unmark trivial states that are accepting, remove useless+unreach SCCs)
    // just in case... usually input is already trim
    if (args.trim)
      ba_trim(aut, log);
    // aut->normalize(); //don't do this, otherwise relationship not clear anymore

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

    //next classify all SCCs
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
    log->info("SCCs: {} total = {} A + {} N + {} M, of which {} D",
              scci.sccs.size(), a, n, m, sccDet.size());

    nba_bitset remain = 0;
    for (auto const& it : sccAcc) {
      nba_bitset const tmp = to_bitset<nba_bitset>(scci.sccs.at(it.first));

      if (dc.sep_rej && it.second == -1) { //if we separate NSCCs
        dc.nscc_states |= tmp;

      } else if (dc.sep_acc && it.second == 1) { //if we separate ASCCs
        dc.ascc_states |= tmp;
        if (dc.sep_acc_cyc)
          dc.asccs_states.push_back(tmp);

      } else if (it.second == 0  //MSCC (and others, when we don't separate them)
          || (!dc.sep_rej && it.second == -1)
          || (!dc.sep_acc && it.second ==  1)) {
        if (dc.sep_mix) { //if we handle (M)SCCs all separately
          //if we optimize deterministic (M)SCCs that are not handled otherwise, sep.
          if (dc.opt_det && contains(sccDet, it.first))
            dc.dscc_states |= tmp;
          else
            dc.msccs_states.push_back(tmp);
        } else {
          remain |= tmp;
        }
      } else {
        log->error("Something went wrong! SCC has invalid acceptance class!");
        exit(1);
      }
    }
    if (dc.sep_acc && !dc.sep_acc_cyc)
      dc.asccs_states.push_back(dc.ascc_states);
    if (!dc.sep_mix)
      dc.msccs_states.push_back(remain);

    // TODO: verify that the interplay is correct
    cerr << dc << endl;
    DetState const st = DetState(dc, 1);
    cerr <<  st << endl;
    auto const sucst = st.succ(dc, 0);
    cerr <<  sucst.first << endl;

    //calculate 2^A and its sccs
    auto const pscon = bench(log,"powerset_construction",
                             WRAP(powerset_construction(aut, dc.aut_mat, dc.aut_asinks)));
    auto const pscon_scci = get_sccs(pscon.states() | ranges::to_vector, aut_succ(pscon));
    log->info("#states in 2^A: {}, #SCCs in 2^A: {}", pscon.num_states(), pscon_scci.sccs.size());
    // print_aut(pscon);

    // -- end of preprocessing --


    // for (auto const scc : sccAcc) {
    //   cerr << scc.first << " (" << seq_to_str(scci.sccs.at(scc.first))
    //       << ") -> " << scc.second << " " << contains(sccDet,scc.first) << endl;
    // }

    // -- begin postprocessing --
    aut.make_complete();
    aut.make_colored();

    assert(aut.is_complete());
    assert(aut.is_colored());
    assert(aut.is_deterministic());

    //TODO: minimize priorities
    //TODO: minimize states
    // -- end of postprocessing --

    return aut;
}

int main(int argc, char *argv[]) {
  // initialize stuff (args + logging)

  // auto console = spd::stdout_color_mt("log");
  auto log = spd::stderr_logger_mt("log");
  spd::set_pattern("[%Y-%m-%d %H:%M:%S %z] [%l] %v");

  auto args = parse_args(argc, argv);
  if (!args.verbose)
    spd::set_level(spd::level::warn);
  else if (args.verbose == 1)
    spd::set_level(spd::level::info);
  else
    spd::set_level(spd::level::debug);

  auto totalstarttime = get_time();

  // now parse input automaton
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
    auto const pa = bench(log,"process_nba",
                          WRAP(process_nba(args, aut, log)));

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

    if (!args.nooutput)
      print_aut(pa);
  }

  /*
  for (auto &aut : auts) {


    //---------------------------

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

    }

  }
    */

  log->info("total time: {:.3f} seconds", get_secs_since(totalstarttime));
  log->info("total used memory: {:.3f} MB", (double)getPeakRSS() / (1024 * 1024));
}
