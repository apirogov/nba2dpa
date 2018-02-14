#include "io.hh"
#include <fstream>
#include <string>
#include <utility>
#include <tuple>
#include "types.hh"

#include <spdlog/spdlog.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "cpphoafparser/consumer/hoa_consumer_print.hh"
#include "cpphoafparser/parser/hoa_parser.hh"
#pragma GCC diagnostic pop

namespace spd = spdlog;

using namespace std;
using namespace cpphoafparser;
using namespace nbautils;

class NBAConsumer : public HOAConsumer {
 private:
  bool eval_expr(label_expr::ptr expr, sym_t val) {
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

 public:
  bool success = false;
  std::map<state_t, std::map<sym_t, std::set<state_t>>> adj;
  BA::uptr aut = make_unique<BA>(BA());

  virtual bool parserResolvesAliases() override { return true; }

  virtual void notifyHeaderStart(const std::string &version) override { ignore = version; }

  virtual void setNumberOfStates(unsigned int numStates) override { ignore = numStates; }

  virtual void addStartStates(const int_list &stateConjunction) override {
    if (stateConjunction.size() != 1)
      throw std::runtime_error("There must be exactly one initial state!");
    aut->init = stateConjunction[0];
  }

  virtual void addAlias(const std::string &name, label_expr::ptr labelExpr) override { ignore = name; ignore = labelExpr; }

  virtual void setAPs(const std::vector<std::string> &aps) override {
    aut->meta.aps = aps;
  }

  virtual void setAcceptanceCondition(unsigned int numSets,
                                      acceptance_expr::ptr accExpr) override {
    if (numSets != 1)
      throw std::runtime_error("There must be exactly one accepting set!");
    ignore = accExpr;
    // TODO: check büchi acceptance!
    // out << "Acceptance: " << numberOfSets << " " << *accExpr << std::endl;
  }

  // TODO: check büchi acceptance!
  virtual void provideAcceptanceName(const std::string &name,
                                     const std::vector<IntOrString> &extraInfo) override {
    ignore = name;
    ignore = extraInfo;
  }

  virtual void setName(const std::string &name) override { aut->meta.name = name; }

  virtual void setTool(const std::string &name,
                       std::shared_ptr<std::string> version) override {
    ignore = name;
    ignore = version;
  }

  virtual void addProperties(const std::vector<std::string> &properties) override { ignore = properties; }

  virtual void addMiscHeader(const std::string &name,
                             const std::vector<IntOrString> &content) override { ignore = name; ignore = content; }

  virtual void notifyBodyStart() override {}

  virtual void addState(unsigned int id, std::shared_ptr<std::string> info,
                        label_expr::ptr labelExpr,
                        std::shared_ptr<int_list> accSignature) override {
    adj[id] = {};
    if (accSignature) aut->acc[id] = Unit();
    ignore = info;
    ignore = labelExpr;
  }

  virtual void addEdgeImplicit(unsigned int stateId, const int_list &conjSuccessors,
                               std::shared_ptr<int_list> accSignature) override {
    /*
    if (accSignature)
      throw std::runtime_error("State-based NBA can not have edge acceptance!");
    if (!aut->meta.num_syms)
      throw std::runtime_error("Edge list, but no atomic propositions!");
    if (conjSuccessors.size() != (size_t)aut->meta.num_syms)
      throw std::runtime_error("Edge list length not equal 2^AP!");
    */
    ignore = stateId;
    ignore = conjSuccessors;
    ignore = accSignature;
    throw std::runtime_error("Implicit edges are not supported!");
  }

  virtual void addEdgeWithLabel(unsigned int stateId, label_expr::ptr labelExpr,
                                const int_list &conjSuccessors,
                                std::shared_ptr<int_list> accSignature) override {
    if (accSignature)
      throw std::runtime_error("State-based NBA can not have edge acceptance!");
    if (!aut->num_syms())
      throw std::runtime_error("Edge list, but no atomic propositions!");

    // store successors in set for now - keeps them unique and sorted
    for (auto sym = 0; sym < (int)aut->num_syms(); sym++) {
      if (eval_expr(labelExpr, sym))
        copy(begin(conjSuccessors), end(conjSuccessors),
             inserter(adj[stateId][sym], end(adj[stateId][sym])));
    }
  }

  virtual void notifyEndOfState(unsigned int stateId) override { ignore = stateId; }

  // finalize - copy successor sets into vectors
  virtual void notifyEnd() override {
    for (auto &it : adj) {
      auto state = it.first;
      aut->adj[state] = {};
      for (auto const &sym : it.second) {
        aut->adj[state][sym.first] = {};
        copy(begin(sym.second), end(sym.second),
             back_inserter(aut->adj[state][sym.first]));
      }
    }
    success = true;
  }

  virtual void notifyAbort() override {}

  virtual void notifyWarning(const std::string &warning) override {
    spd::get("log")->warn(warning);
    // std::cerr << "Warning: " << warning << std::endl;
  }
};

namespace nbautils {

BA::uptr parse_ba(string const &filename) {
  spd::get("log")->debug("parse_ba({}) ", filename);

  bool givenfile = !filename.empty();
  ifstream fin;
  if (givenfile) fin = ifstream(filename.c_str(), ifstream::in);
  try {
    HOAConsumer::ptr consumer(new NBAConsumer());
    auto nc = std::static_pointer_cast<NBAConsumer>(consumer);
    HOAParser::parse(givenfile ? fin : cin, consumer);
    return move(nc->aut);
  } catch (std::exception &e) {
    spd::get("log")->error(e.what());
    // std::cerr << e.what() << std::endl;
    return nullptr;
  }
}
}  // namespace nbautils
