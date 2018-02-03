#include <vector>

#include "dpa.hh"
#include "level.hh"
#include "nba.hh"
#include "relorder.hh"

using namespace std;
using namespace nbautils;

namespace nbautils {

DPA::DPA(NBA const &nba) : nba_(make_shared<NBA>(nba)) {}

// Level DPA::succ(sym_t const &x) {
// 	return Level();
// }

}  // namespace nbautils
