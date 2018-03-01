#pragma once

#include <string>
#include <iostream>
#include "ba.hh"
#include "det.hh"

#include <spdlog/spdlog.h>

namespace nbautils {

std::vector<BA::uptr> parse_hoa_ba(std::string const& filename, std::shared_ptr<spdlog::logger> log=nullptr);

std::string sym_to_edgelabel(sym_t s, std::vector<std::string> const& aps, bool as_aps=false);

//output automaton in HOA format
template<Acceptance A, typename T>
void print_hoa(SWA<A,T> const& aut, ostream &out = cout) {
  out << "HOA: v1" << endl;
  out << "name: \"" << aut.get_name() << "\"" << endl;
  out << "States: " << aut.num_states() << endl;
  out << "Start: " << seq_to_str(aut.get_init(), "&") << endl;
  out << "AP: " << aut.get_aps().size();
  for (auto const& ap : aut.get_aps())
    out << " \"" << ap << "\"";
  out << endl;

  if (A==Acceptance::BUCHI) {
    out << "acc-name: " << "Buchi" << endl << "Acceptance: 1 Inf(0)" << endl;
  } else if (A==Acceptance::PARITY) {
    vector<acc_t> accs = aut.get_accsets();
    int pris = accs.back() + 1;
    out << "acc-name: " << "parity min even " << pris << endl;
    out << "Acceptance: " << pris << " ";
    // for (auto i = 0; i<(int)aut.accsets.size(); i++) {
    //     auto a = accs[i];
    for (auto i = 0; i<pris; i++) {
        auto a = i;
        out << (a%2==0 ? "Inf" : "Fin") << "(" << a << ")";
        if (i!=pris-1) {
          // out << " ";
          out << (a%2==0 ? "|" : "&") << "(";
        }
    }
    for (auto i = 0; i<pris-1; i++)
      out << ")";
    out << endl;
  } else {
    out << "acc-name: none" << endl << "Acceptance: 0 f" << endl;
  }

  out << "properties: trans-labels explicit-labels state-acc";
  if (is_deterministic(aut))
    out << " deterministic";
  if (is_colored(aut))
    out << " colored";
  // if (is_complete(aut))
  //   out << " complete";
  out << endl;

  out << "--BODY--" << endl;
  for (auto const& p : aut.states()) {

    //State: x "tagstr" {accs}
    out << "State: " << p;
    if (aut.tag->hasi(p)) {
      out << " \"" << aut.print_tag(p) << "\"";
    }
    if (aut.has_accs(p)) {
      auto ac = aut.get_accs(p);
      out << " {";
      for (auto i=0; i<(int)ac.size(); i++) {
        out << ac[i];
        if (i!=(int)ac.size()-1)
          out << " ";
      }
      out << "}";
    }
    out << endl;

    //list edges
    for (auto s : aut.outsyms(p)) {
      // cout << (int)s << endl;
      string const lbl = sym_to_edgelabel(s, aut.get_aps());
      for (auto q : aut.succ(p,s))
        out << "[" << lbl << "] " << q << endl;
    }
  }

  out << "--END--" << endl;
}

}  // namespace nbautils
