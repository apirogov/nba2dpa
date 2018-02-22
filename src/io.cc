#include "io.hh"
#include <fstream>
#include <memory>
#include <string>
#include <utility>
#include <tuple>
#include "types.hh"
#include "det.hh"

#include <spdlog/spdlog.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "cpphoafparser/consumer/hoa_consumer_print.hh"
#include "cpphoafparser/parser/hoa_parser.hh"
#include "cpphoafparser/util/implicit_edge_helper.hh"
#pragma GCC diagnostic pop

using namespace std;
using namespace cpphoafparser;
using namespace nbautils;

// helper:
// take a boolean expression (from parser) and a symbol (bits represent ap)
// return whether the represented valuation satisfies the formula
bool eval_expr(BooleanExpression<AtomLabel>::ptr expr, sym_t val) {
  unsigned apidx = 0;
  switch (expr->getType()) {
    case BooleanExpression<AtomLabel>::OperatorType::EXP_TRUE:
      return true;
    case BooleanExpression<AtomLabel>::OperatorType::EXP_FALSE:
      return false;
    case BooleanExpression<AtomLabel>::OperatorType::EXP_ATOM:
      apidx = expr->getAtom().getAPIndex();
      return (val >> apidx) & 1;
    case BooleanExpression<AtomLabel>::OperatorType::EXP_NOT:
      return !eval_expr(expr->getLeft(), val);
    case BooleanExpression<AtomLabel>::OperatorType::EXP_AND:
      return eval_expr(expr->getLeft(), val) && eval_expr(expr->getRight(), val);
    case BooleanExpression<AtomLabel>::OperatorType::EXP_OR:
      return eval_expr(expr->getLeft(), val) || eval_expr(expr->getRight(), val);
  }
  //can not happen
  throw std::runtime_error("There must be an unhandled BooleanExpression eval case! FIXME");
  return false;
}

// our "parser client" constructing automata from parser events
class NBAConsumer : public HOAConsumer {
 public:
  AcceptanceRepositoryStandard acrep;
  std::unique_ptr<ImplicitEdgeHelper> helper;

  std::shared_ptr<spdlog::logger> log;

  vector<BA::uptr> auts;
  std::map<state_t, std::map<sym_t, std::set<state_t>>> adj;


  NBAConsumer(std::shared_ptr<spdlog::logger> l) : log(l) {};

  virtual bool parserResolvesAliases() override { return true; }

  virtual void notifyHeaderStart(const std::string &version) override {
    ignore = version;
    auts.push_back(move(make_unique<BA>()));
  }

  virtual void setNumberOfStates(unsigned int numStates) override { ignore = numStates; }

  virtual void addStartStates(const int_list &stateConjunction) override {
    if (stateConjunction.size() != 1)
      throw std::runtime_error("There must be exactly one initial state!");
    auts.back()->set_init(stateConjunction);
  }

  virtual void addAlias(const std::string &name, label_expr::ptr labelExpr) override { ignore = name; ignore = labelExpr; }

  virtual void setAPs(const std::vector<std::string> &aps) override {
    auts.back()->set_aps(aps);
  }

  virtual void setAcceptanceCondition(unsigned int numSets,
                                      acceptance_expr::ptr accExpr) override {
    if (numSets != 1)
      throw std::runtime_error("There must be exactly one accepting set!");
    ignore = accExpr;
  }

  virtual void provideAcceptanceName(const std::string &name,
                                     const std::vector<IntOrString> &extraInfo) override {
    ignore = extraInfo;
    if (name != acrep.ACC_BUCHI)
      throw std::runtime_error("Automaton does not have Büchi acceptance!");
  }

  virtual void setName(const std::string &name) override { auts.back()->set_name(name); }

  virtual void setTool(const std::string &name,
                       std::shared_ptr<std::string> version) override {
    ignore = name;
    ignore = version;
  }

  virtual void addProperties(const std::vector<std::string> &properties) override { ignore = properties; }

  virtual void addMiscHeader(const std::string &name,
                             const std::vector<IntOrString> &content) override { ignore = name; ignore = content; }

  virtual void notifyBodyStart() override {
    helper = make_unique<ImplicitEdgeHelper>(auts.back()->get_aps().size());
  }

  virtual void addState(unsigned int id, std::shared_ptr<std::string> info,
                        label_expr::ptr labelExpr,
                        std::shared_ptr<int_list> accSignature) override {
    helper->startOfState(id);

    if (!auts.back()->has_state(id))
      auts.back()->add_state(id);
    auts.back()->tag->put(id, id); //we mark the original id number in the tag
    if (accSignature) auts.back()->set_accs(id,{0});
    ignore = info;
    ignore = labelExpr;
  }

  virtual void addEdgeImplicit(unsigned int stateId, const int_list &conjSuccessors,
                               std::shared_ptr<int_list> accSignature) override {
    if (accSignature)
      throw std::runtime_error("State-based NBA can not have edge acceptance!");
    if (!auts.back()->num_syms())
      throw std::runtime_error("Edge list, but no atomic propositions!");

    auto sym = helper->nextImplicitEdge(); //the helper tells us which edge this is
    if (!conjSuccessors.empty())
        copy(cbegin(conjSuccessors), cend(conjSuccessors),
             inserter(adj[stateId][sym], end(adj[stateId][sym])));
  }

  virtual void addEdgeWithLabel(unsigned int stateId, label_expr::ptr labelExpr,
                                const int_list &conjSuccessors,
                                std::shared_ptr<int_list> accSignature) override {
    if (accSignature)
      throw std::runtime_error("State-based NBA can not have edge acceptance!");
    if (!auts.back()->num_syms())
      throw std::runtime_error("Edge list, but no atomic propositions!");

    // check boolean expression against all combinations and add successors accordingly
    // store successors in set for now - keeps them unique and sorted
    for (auto sym = 0; sym < (int)auts.back()->num_syms(); sym++) {
      if (eval_expr(labelExpr, sym) && !conjSuccessors.empty())
        copy(cbegin(conjSuccessors), cend(conjSuccessors),
             inserter(adj[stateId][sym], end(adj[stateId][sym])));
    }
  }

  virtual void notifyEndOfState(unsigned int stateId) override {
    ignore = stateId;
    helper->endOfState();
  }

  // finalize - copy successor edge sets into vectors
  virtual void notifyEnd() override {
    for (auto &it : adj) {
      auto state = it.first;
      for (auto const &sym : it.second) { //add edges, if present
        auts.back()->adj[state][sym.first] = {};
        copy(cbegin(sym.second), cend(sym.second),
             back_inserter(auts.back()->adj[state][sym.first]));
      }
    }
    adj = {}; //reset edges storage for next automaton
  }

  virtual void notifyAbort() override {
    auts.pop_back(); //drop current automaton
    adj = {}; //reset edges
  }

  virtual void notifyWarning(const std::string &warning) override {
    if (log) log->warn(warning);
    // std::cerr << "Warning: " << warning << std::endl;
  }
};

namespace nbautils {

// Try to read an arbitrary number of NBAs from given stream. If anything goes wrong, return nothing.
vector<BA::uptr> parse_hoa_ba(istream& instream, std::shared_ptr<spdlog::logger> log) {
  try {
    auto nc = make_shared<NBAConsumer>(log);
    while (instream.good()) {
      HOAParser::parse(instream, nc);
      //hack required, as parser consumes a token from next automaton too, for some reason Oo
      if (instream.good()) { //if it looks like there comes another automaton, put eaten token back
        instream.putback(' '); instream.putback(':');
        instream.putback('A'); instream.putback('O'); instream.putback('H');
      }
    }
    return move(nc->auts);

  } catch (std::exception &e) {
    if (log)
      log->error(e.what());
    else
      std::cerr << e.what() << std::endl;
    return {};
  }
}

//empty string as input file -> read stdin
vector<BA::uptr> parse_hoa_ba(string const &filename, std::shared_ptr<spdlog::logger> log) {
  bool givenfile = !filename.empty();
  ifstream fin;
  if (givenfile) fin = ifstream(filename.c_str(), ifstream::in);
  auto& instream = givenfile ? fin : cin;
  return parse_hoa_ba(instream, log);
}

std::string sym_to_edgelabel(sym_t s, std::vector<std::string> const& aps,bool as_aps) {
  stringstream elbl;
  auto tmp=s;
  for (size_t b=0; b<aps.size(); b++) {
    if (!(tmp&1))
      elbl << "!";
    if (as_aps)
      elbl << aps[b];
    else
      elbl << b;
    if (b<aps.size()-1)
      elbl << "&";
    tmp = tmp >> 1;
  }
  return elbl.str();
}

//output automaton in HOA format
void print_hoa_pa(PA const& aut, ostream &out) {
  out << "HOA: v1" << endl;
  out << "name: \"" << aut.get_name() << "\"" << endl;
  out << "States: " << aut.num_states() << endl;
  out << "Start: ";
  auto const& init = aut.get_init();
  for (size_t i=0; i<init.size(); i++)
    out << init.at(i) << (i<init.size()-1 ? "&" : "");
  out << endl;
  out << "AP: " << aut.get_aps().size();
  for (auto const& ap : aut.get_aps())
    out << " \"" << ap << "\"";
  out << endl;

  vector<acc_t> accs(cbegin(aut.get_accsets()),cend(aut.get_accsets()));
  int pris = accs.back() + 1;
  out << "acc-name: " << "parity min even " << pris << endl;
  out << "Acceptance: " << pris << " ";
  // for (auto i = 0; i<(int)aut.accsets.size(); i++) {
  //     auto a = accs[i];
  for (auto i = 0; i<pris; i++) {
      auto a = i;
      out << (a%2==0 ? "Inf" : "Fin") << "(" << a << ")";
      if (i!=pris-1) {
        out << " ";
        out << (a%2==0 ? "|" : "&") << " (";
      }
  }
  for (auto i = 0; i<pris-1; i++)
    out << ")";
  out << endl;

  out << "properties: trans-labels explicit-labels state-acc colored deterministic" << endl;

  out << "--BODY--" << endl;
  for (auto &it : aut.adj) {
    auto p = it.first;

    //State: x "tagstr" {accs}
    out << "State: " << p;
    if (aut.tag->hasi(p)) {
      out << " \"" << aut.tag->geti(p) << "\"";
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
