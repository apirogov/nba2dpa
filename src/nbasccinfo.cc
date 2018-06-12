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


int main(int argc, char *argv[]) {
  // initialize stuff (args + logging)

  // auto console = spd::stdout_color_mt("log");
  auto log = spd::stderr_logger_mt("log");
  spd::set_pattern("[%Y-%m-%d %H:%M:%S %z] [%l] %v");

  // now parse input automaton
  auto auts = nbautils::parse_hoa("", log);

  if (auts.empty()) {
    log->error("Parsing NBAs from stdin failed!");
    exit(1);
  }

  for (auto &aut : auts) {
    if (aut->acond != Acceptance::BUCHI) {
      log->warn("skipping non-Büchi automaton with name \"{}\"", aut->get_name());
      continue;
    }

    int num_states = aut->num_states();
    // auto const num_useless = trim_ba(*aut);

    auto aut_st = aut->states();
    succ_fun<state_t> const aut_sucs = [&aut](state_t v){ return aut->succ(v); };
    function<bool(state_t)> const aut_acc = [&aut](state_t v){ return aut->has_accs(v); };
    succ_sym_fun<state_t, sym_t> const aut_xsucs = [&aut](state_t v,sym_t s){ return aut->succ(v,s); };
    outsym_fun<state_t,sym_t> const aut_osyms = [&aut](state_t v){ return aut->outsyms(v); };

    // detect accepting sinks (acc states with self loop for each sym)
    vector<small_state_t> accsinks;
    accsinks = to_small_state_t(get_accepting_sinks(aut_st, aut->num_syms(), aut_acc, aut_osyms, aut_xsucs));
    int num_sinks = accsinks.size();

    // calculate SCCs of NBA if needed
    SCCDat<state_t>::uptr aut_scc = nullptr;
    BaSccClassification::uptr aut_cl = nullptr;
    aut_scc = get_sccs(aut_st, aut_sucs);
    aut_cl = ba_classify_sccs(*aut_scc, aut_acc);
    int num_sccs = aut_scc->sccs.size();
    int num_asccs = aut_cl->accepting.size();
    int num_nsccs = aut_cl->rejecting.size();
    int num_msccs = num_sccs - num_asccs - num_nsccs;

    // calculate 2^A (with accsinks opt)
    PS::uptr ps = nullptr;
    SCCDat<state_t>::uptr psi = nullptr;
    ps = powerset_construction(*aut, accsinks);
    succ_fun<state_t> pssucs = [&ps](state_t v){ return ps->succ(v); };
    psi = get_sccs(ps->states(), pssucs);
    int num_psets = ps->num_states();
    int num_psetsccs = psi->sccs.size();

    //TODO: impression of size distribution of SCCs

    cout << "St:" << num_states
         // << "\tUS:" << num_useless
         << "\tAS:" << num_sinks
         << "\tSCC:" << num_sccs
         << "\tNSCC:" << num_nsccs
         << "\tASCC:" << num_asccs
         << "\tMSCC:" << num_msccs
         << "\tPS:" << num_psets
         << "\tPSCC:" << num_psetsccs
         << endl;
  }
}
