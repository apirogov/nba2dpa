#include "aut.hh"
#include "io.hh"

#include <fstream>
#include <memory>
#include <string>

#include <spdlog/spdlog.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "cpphoafparser/parser/hoa_parser.hh"
// #include "cpphoafparser/consumer/hoa_consumer_print.hh" //demo: prints back hoa
#include "cpphoafparser/util/implicit_edge_helper.hh"
#pragma GCC diagnostic pop

using namespace std;
using namespace cpphoafparser;
using namespace nbautils;

namespace nbautils {

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
      return  eval_expr(expr->getLeft(), val) && eval_expr(expr->getRight(), val);
    case BooleanExpression<AtomLabel>::OperatorType::EXP_OR:
      return  eval_expr(expr->getLeft(), val) || eval_expr(expr->getRight(), val);
  }
  //can not happen
  throw std::runtime_error("There must be an unhandled BooleanExpression eval case! FIXME");
  return false;
}

//decode conjunction of APs from a symbol number
std::string sym_to_edgelabel(sym_t s, std::vector<std::string> const& aps,bool as_aps) {
  if (aps.empty()) //no aps and edge -> can always be taken
    return "t";

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

}
