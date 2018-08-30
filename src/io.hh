#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <memory>

#include <spdlog/spdlog.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "cpphoafparser/parser/hoa_parser.hh"
#include "cpphoafparser/consumer/hoa_consumer_null.hh" //HOAConsumerNull dummy
#pragma GCC diagnostic pop

// #include <boost/optional.hpp>

#include "aut.hh"
#include "pa.hh"

namespace nbautils {

using namespace std;
using namespace cpphoafparser;

std::string sym_to_edgelabel(sym_t s, std::vector<std::string> const& aps, bool as_aps=false);
bool eval_expr(BooleanExpression<AtomLabel>::ptr expr, sym_t val);

// add function to retrieve something after parsing an automaton
// i.e. can be used to read automaton structure or just compute something
template<typename T>
class MyConsumer : public HOAConsumer {
public:
  T retrieve();
};

// e.g. this dummy consumer extracts the number of states only
template<>
class MyConsumer<int> : public HOAConsumerNull {
  std::shared_ptr<spdlog::logger> log;
  int st = 0;
public:
  MyConsumer(std::shared_ptr<spdlog::logger> l) : log(l) {};
  virtual void setNumberOfStates(unsigned int num) override { st = num; }
  int retrieve() { return st; };
};

// our "parser client" constructing automata from parser events
template<>
class MyConsumer<Aut<string>> : public HOAConsumerNull {
 public:
  AcceptanceRepositoryStandard acrep;
  std::unique_ptr<ImplicitEdgeHelper> helper;

  std::shared_ptr<spdlog::logger> log;

  Aut<string> aut;
  bool fixed_sbatba = false;

  MyConsumer(std::shared_ptr<spdlog::logger> l) : log(l) {};
  Aut<string> retrieve() { return aut; }

  virtual bool parserResolvesAliases() override { return true; }
  virtual void setName(const std::string &name) override { aut.set_name(name); }
  virtual void setAPs(const std::vector<std::string> &aps) override { aut.set_aps(aps); }

  virtual void addStartStates(const int_list &stConj) override {
    if (stConj.size() != 1)
      throw std::runtime_error("There must be exactly one initial state!");

    if (!aut.has_state(stConj[0]))
      aut.add_state(stConj[0]);
    aut.set_init(stConj[0]);
  }

  virtual void provideAcceptanceName(const std::string &name,
                                     const std::vector<IntOrString> &extraInfo) override {
    ignore = extraInfo;
    if (name != acrep.ACC_BUCHI && name != acrep.ACC_PARITY)
      throw std::runtime_error("Automaton does not have Büchi or Parity acceptance!");

    if (name == acrep.ACC_PARITY) {
      PAType pat = PAType::MIN_EVEN;
      if (extraInfo.at(0).getString() != "min")
        pat = opposite_polarity(pat);
      if (extraInfo.at(1).getString() != "even")
        pat = opposite_parity(pat);
      aut.set_patype(pat);
    }
  }

  virtual void notifyBodyStart() override {
    helper = make_unique<ImplicitEdgeHelper>(aut.get_aps().size());
  }

  virtual void addState(unsigned int id, std::shared_ptr<std::string> info,
                        label_expr::ptr labelExpr,
                        std::shared_ptr<int_list> accSig) override {
    ignore = labelExpr;
    if (accSig && accSig->size()>1)
      throw std::runtime_error("State can have only one priority mark!");

    helper->startOfState(id);

    if (!aut.has_state(id)) {
      aut.add_state(id);
    }

    if (accSig) {
      if (fixed_sbatba && !aut.is_sba())
        throw std::runtime_error("Priorities should be either at states or at edges!");

      if (!fixed_sbatba) {
        aut.set_sba(true);
        fixed_sbatba = true;
      }

      aut.set_pri(id,(*accSig)[0]);
    }

    if (info) {
      aut.tag.put(*info, id);
    } else {
      aut.tag.put(to_string(id), id); //if untagged, we tag with the original id number
    }
  }

  virtual void addEdgeImplicit(unsigned int sId, const int_list &conjSucs,
                               std::shared_ptr<int_list> accSig) override {
    if (!aut.num_syms())
      throw std::runtime_error("Edge list, but no atomic propositions!");
    if (accSig && accSig->size()>1)
      throw std::runtime_error("Transition can have only one priority mark!");
    if (accSig && !fixed_sbatba)
      fixed_sbatba = true;
    if (fixed_sbatba && accSig && aut.is_sba())
      throw std::runtime_error("Priorities should be either at states or at edges!");

    auto const pri = accSig ? accSig->at(0) : -1;

    //the helper tells us which edge this is
    auto const sym = helper->nextImplicitEdge();
    for (auto const trg : conjSucs) {
      if (!aut.has_state(trg))
        aut.add_state(trg);
      aut.add_edge(sId,sym,trg,pri);
    }
  }

  virtual void addEdgeWithLabel(unsigned int sId, label_expr::ptr labelExpr,
                                const int_list &conjSucs,
                                std::shared_ptr<int_list> accSig) override {
    if (!aut.num_syms())
      throw std::runtime_error("Edge list, but no atomic propositions!");
    if (accSig && accSig->size()>1)
      throw std::runtime_error("Transition can have only one priority mark!");
    if (accSig && !fixed_sbatba)
      fixed_sbatba = true;
    if (fixed_sbatba && accSig && aut.is_sba())
      throw std::runtime_error("Priorities should be either at states or at edges!");

    auto const pri = accSig ? accSig->at(0) : -1;

    // check boolean expression against all combinations and add successors accordingly
    // store successors in set for now - keeps them unique and sorted
    for (sym_t const sym : aut.syms()) {
      if (eval_expr(labelExpr, sym) && !conjSucs.empty())
        for (auto const trg : conjSucs) {
          if (!aut.has_state(trg))
            aut.add_state(trg);
          aut.add_edge(sId,sym,trg,pri);
        }
    }
  }

  virtual void notifyEndOfState(unsigned int stateId) override {
    ignore = stateId;
    helper->endOfState();
  }

  // finalize - copy successor edge sets into vectors
  virtual void notifyEnd() override {
    if (!fixed_sbatba) { //automaton without acceptance sets can be seen as statebased
      aut.set_sba(true);
    }
    aut.tag_to_str = default_printer<string>();
  }
  //   for (auto &it : edges) {
  //     auto state = it.first;
  //     for (auto const &sym : it.second) { //add edges, if present
  //       auts.back()->set_succs(state, sym.first,
  //                              vector<state_t>(cbegin(sym.second),cend(sym.second)));
  //     }
  //   }
  // }

  // virtual void notifyAbort() override { }

  virtual void notifyWarning(const std::string &warning) override {
    if (log)
      log->warn(warning);
    else
      std::cerr << "Warning: " << warning << std::endl;
  }
};

template<typename T>
class AutStream {
  bool garbage = false;
  bool is_file = false;
  ifstream fin;
  istream* in;
  std::shared_ptr<spdlog::logger> log;

public:
  AutStream(istream &instream, std::shared_ptr<spdlog::logger> logger=nullptr) {
    log = logger;
    in = instream;
  }

  AutStream(string const &filename, std::shared_ptr<spdlog::logger> logger=nullptr) {
    log = logger;
    is_file = !filename.empty();
    if (is_file)
      fin = ifstream(filename.c_str(), ifstream::in);
    auto& instream = is_file ? fin : cin;
    in = &instream;
  }

  bool has_next() const { return in->good() && !garbage; }

  // return next automaton. can fail and will return "default value"
  T parse_next() {
    try {
      auto mc = make_shared<MyConsumer<T>>(log);

      HOAParser::parse(*in, mc);
      //hack required, as parser consumes a token from next automaton too, for some reason Oo
      if (in->good()) {
        //if it looks like there comes another automaton, put eaten token back
        //I have no idea why neither unget nor putback works generically on both kinds...
        if (is_file) {
          in->unget(); in->unget(); in->unget(); in->unget(); in->unget(); //unget "HOA: "
        } else {
          //putback "HOA: "
          in->putback(' '); in->putback(':');
          in->putback('A'); in->putback('O'); in->putback('H');
        }
      }

      return mc->retrieve();
    }

    catch (char const* e) { //aborted and did not deliver next. not fatal
      if (log)
        log->error(e);
      else
        std::cerr << e << std::endl;
      return {};
    }

    catch (std::exception &e) { //syntax error stuff. rather fatal. we give up.
      if (log)
        log->error(e.what());
      else
        std::cerr << e.what() << std::endl;
      garbage = true; //input stream is broken, we better abort
      return {};
    }

  }
};


//output automaton in HOA format with parity min even acceptance
template<typename T>
void print_aut(Aut<T> const& aut, ostream &out = cout) {
  assert(aut.get_patype() == PAType::MIN_EVEN);
  bool sba = aut.is_sba();

  out << "HOA: v1" << endl;
  out << "name: \"" << aut.get_name() << "\"" << endl;
  out << "States: " << aut.num_states() << endl;
  out << "Start: " << aut.get_init() << endl;
  out << "AP: " << aut.get_aps().size();
  for (auto const& ap : aut.get_aps())
    out << " \"" << ap << "\"";
  out << endl;

  if (aut.pris().size() > 0) {
    int pris = aut.pris().back() + 1;
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

  out << "properties: trans-labels explicit-labels";
  if (sba)
    out << " state-acc";
  if (aut.is_deterministic())
    out << " deterministic";
  if (aut.is_colored())
    out << " colored";
  if (aut.is_complete())
    out << " complete";
  out << endl;

  out << "--BODY--" << endl;
  for (auto const& p : aut.states()) {

    //State: x "tagstr" {accs}
    out << "State: " << p;
    if (aut.tag.hasi(p)) {
      out << " \"";
      aut.print_state_tag(out,p);
      out << "\"";
    }
    if (sba && aut.has_pri(p)) {
      out << " {" << aut.get_pri(p) << "}";
    }
    out << endl;

    //list edges
    for (auto s : aut.state_outsyms(p)) {
      // cout << (int)s << endl;
      string const lbl = sym_to_edgelabel(s, aut.get_aps());
      for (auto e : aut.succ_edges(p,s)) {
        out << "[" << lbl << "] " << e.first;
        if (!sba && e.second >= 0)
          out << " {" << e.second << "}";
        out << endl;
      }
    }
  }

  out << "--END--" << endl;
}

}  // namespace nbautils
