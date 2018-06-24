#include "aut.hh"
#include "detstate.hh"
#include "common/util.hh"
#include <iostream>
#include <algorithm>
#include <bitset>
#include <iostream>
#include <sstream>
#include <vector>
using namespace std;
using namespace nbautils;

namespace nbautils {

//precompute all sets often used in successor calculation
DetConfSets calc_detconfsets(DetConf const& dc, SCCDat const& scci,
    BASccAClass const& sccAcc, set<unsigned> const& sccDet) {
  DetConfSets ret;

  nba_bitset remain = 0;
  for (auto const& it : sccAcc) {
    nba_bitset const tmp = to_bitset<nba_bitset>(scci.sccs.at(it.first));

    if (dc.sep_rej && it.second == -1) { //if we separate NSCCs
      ret.nscc_states |= tmp;

    } else if (dc.sep_acc && it.second == 1) { //if we separate ASCCs
      ret.ascc_states |= tmp;
      if (dc.sep_acc_cyc)
        ret.asccs_states.push_back(tmp);

    } else if (it.second == 0  //MSCC (and others, when we don't separate them)
        || (!dc.sep_rej && it.second == -1)
        || (!dc.sep_acc && it.second ==  1)) {
      if (dc.sep_mix) { //if we handle (M)SCCs all separately
        //if we optimize deterministic (M)SCCs that are not handled otherwise, sep.
        if (dc.opt_det && contains(sccDet, it.first))
          ret.dscc_states |= tmp;
        else
          ret.msccs_states.push_back(tmp);
      } else {
        remain |= tmp;
      }
    } else {
      cerr << "Something went wrong! SCC has invalid acceptance class!" << endl;
      exit(1);
    }
  }
  if (dc.sep_acc && !dc.sep_acc_cyc)
    ret.asccs_states.push_back(ret.ascc_states);
  if (!dc.sep_mix)
    ret.msccs_states.push_back(remain);

  return ret;
}

std::ostream& operator<<(std::ostream& os, DetConf const& dc) {
  os << "DetConf {" << endl;
  os << "aut_mat: " << !dc.aut_mat.empty() << endl;
  os << "aut_states: " << pretty_bitset(dc.aut_states) << endl;
  os << "aut_acc: "    << pretty_bitset(dc.aut_acc) << endl;
  os << "aut_asinks: " << pretty_bitset(dc.aut_asinks) << endl;
  os << "ctx: "        << !dc.ctx.empty() << endl;

  os << "nscc_states: " <<  pretty_bitset(dc.sets.nscc_states) << endl;
  os << "ascc_states: " <<  pretty_bitset(dc.sets.ascc_states) << endl;
  for (auto const& as : dc.sets.asccs_states)
    os << "\tascc: " <<  pretty_bitset(as) << endl;
  os << "dscc_states: " <<  pretty_bitset(dc.sets.dscc_states) << endl;
  os << "msccs:";
  for (auto const& ms : dc.sets.msccs_states)
    os << "\t" <<  pretty_bitset(ms) << endl;

  os << "update: " << (int)dc.update << endl;
  os << "options: ";
  if (dc.sep_rej)     os << "seprej ";
  if (dc.sep_acc)     os << "sepacc ";
  if (dc.sep_acc_cyc) os << "sepacccyc ";
  if (dc.sep_mix)     os << "sepmix ";
  if (dc.opt_det)     os << "optdet ";
  if (dc.weaksat)     os << "weaksat ";
  if (dc.puretrees)   os << "puretrees ";
  os << endl;
  os << "}" << endl;
  return os;
}

std::ostream& operator<<(std::ostream& os, DetState const& s) {
  os << "N: " << pretty_bitset(s.nsccs)
    << "\t(AC: " << pretty_bitset(s.asccs)
    << ", AB: " << pretty_bitset(s.asccs_buf) << "):" << s.asccs_pri
    << "\tD: (";
  vector<string> tmp;
  for (auto const& dscc : s.dsccs)
    tmp.push_back(pretty_bitset(dscc.first) + ":" + to_string(dscc.second));
  os << seq_to_str(tmp, ", ") << ")";
  tmp.clear();
  os << "\tM: (";
  for (auto const& mscc : s.msccs) {
    vector<string> tmp2;
    for (auto const& ms : mscc) {
      tmp2.push_back(pretty_bitset(ms.first) + ":" + to_string(ms.second));
    }
    tmp.push_back(seq_to_str(tmp2, ", "));
  }
  os << seq_to_str(tmp, " | ") << ")";
  return os;
}

string DetState::to_string() const {
  stringstream ss;
  ss << *this;
  return ss.str();
}

//componentwise equality
bool DetState::operator==(DetState const& o) const {
  return powerset == o.powerset
    &&      nsccs == o.nsccs
    &&  asccs_buf == o.asccs_buf
    &&      asccs == o.asccs
    &&  asccs_pri == o.asccs_pri
    &&      dsccs == o.dsccs
    &&      msccs == o.msccs
    ;
}

//returns Parent and Left-border relationship for a list of ranks
//(left border = left sibling for the last one popped in the inner while loop)
pair<vector<int>,vector<int>> unflatten(vector<pair<nba_bitset, pri_t>> const& rank) {
  int const n = rank.size();
  pair<vector<int>,vector<int>> ret = make_pair(vector<int>(n),vector<int>(n));
  stack<int> s;
  for (int i=n-1; i>=0; i--) {
    while (!s.empty() && rank[i].second < rank[s.top()].second) {
      ret.second[s.top()] = i;
      s.pop();
    }
    ret.first[i] = s.empty() ? -1 : s.top();
    s.push(i);
  }
  while (!s.empty()) {
    ret.second[s.top()] = -1;
    s.pop();
  }
  return ret;
}

DetState::DetState() {}

//put everything into corresponding setsn
DetState::DetState(DetConf const& dc, nba_bitset const& qs) {
  powerset = qs;
  nsccs = qs & dc.sets.nscc_states;
  asccs_buf = qs & dc.sets.ascc_states;

  int cur_fresh = 1;
  nba_bitset tmp = qs & dc.sets.dscc_states;
  if (tmp != 0)
    dsccs.push_back(make_pair(tmp, cur_fresh++));

  int const numm = dc.sets.msccs_states.size();
  msccs.resize(numm); //allocate as many as MSCCs
  for (auto const i : ranges::view::ints(0, numm)) {
    tmp = qs & dc.sets.msccs_states.at(i);
    if (tmp != 0)
      msccs.at(i).push_back(make_pair(tmp, cur_fresh++));
  }
}

// ----------------------------------------------------------------------------

//apply successor set function on each set separately, inplace
void successorize_all(DetConf const& dc, DetState& s, sym_t const x) {
  auto const psucc = [&dc,x](auto const& bset){
    return powersucc(dc.aut_mat, bset, x, dc.aut_asinks); };

  s.powerset  = psucc(s.powerset);
  s.nsccs     = psucc(s.nsccs);
  s.asccs_buf = psucc(s.asccs_buf);
  s.asccs     = psucc(s.asccs);
  for (auto const i : ranges::view::ints(0, (int)s.dsccs.size()))
    s.dsccs[i].first = psucc(s.dsccs[i].first);
  for (auto const i : ranges::view::ints(0, (int)s.msccs.size()))
    for (auto const j : ranges::view::ints(0, (int)s.msccs[i].size()))
      s.msccs[i][j].first = psucc(s.msccs[i][j].first);
}

//remove wrong located states, return their collection
//the sets are provided separately as they might be modified for context
nba_bitset extract_switchers(DetConf const& dc, DetConfSets const& sts, DetState& s) {
  nba_bitset switchers;

  switchers |= s.nsccs     & (~sts.nscc_states & dc.aut_states);
  switchers |= s.asccs     & (~sts.ascc_states & dc.aut_states);
  switchers |= s.asccs_buf & (~sts.ascc_states & dc.aut_states);
  s.nsccs     &= sts.nscc_states;
  s.asccs     &= sts.ascc_states;
  s.asccs_buf &= sts.ascc_states;

  for (auto const i : ranges::view::ints(0, (int)s.dsccs.size())) {
    switchers |= s.dsccs[i].first & (~sts.dscc_states & dc.aut_states);
    s.dsccs[i].first &= sts.dscc_states;
  }
  for (auto const i : ranges::view::ints(0, (int)s.msccs.size()))
    for (auto const j : ranges::view::ints(0, (int)s.msccs[i].size())) {
      switchers |= s.msccs[i][j].first & (~sts.msccs_states[i] & dc.aut_states);
      s.msccs[i][j].first &= sts.msccs_states[i];
    }
  return switchers;
}

//given a safra forest tuple, split accepting states into fresh children with fresh ranks
void expand_row(DetConf const& dc, vector<pair<nba_bitset, pri_t>>& row, pri_t& cur_fresh) {
  vector<pair<nba_bitset, pri_t>> nrow;
  nrow.reserve(2*row.size());
  for (auto& it : row) {
    nrow.push_back(make_pair(it.first &  dc.aut_acc, cur_fresh++));
    nrow.push_back(make_pair(it.first & ~dc.aut_acc, it.second));
  }
  nrow.shrink_to_fit();
  swap(row, nrow);
}

//split acc successors into extra nodes for tree-organized sets (MSCCs)
void expand_trees(DetConf const& dc, DetState &s, pri_t& cur_fresh) {
  for (auto& mscc : s.msccs)
    expand_row(dc, mscc, cur_fresh);
}

//given a ranked tuple, keep leftmost occurence of each state
void left_normalize_row(vector<pair<nba_bitset, pri_t>>& row) {
  nba_bitset seen = 0;
  for (auto& s : row) {
    s.first &= ~seen;
    seen |= s.first;
  }
}

//keep leftmost occurence of each state
void left_normalize(DetState &s) {
  s.asccs_buf &= ~s.asccs; //keep in buffer only ones not already reached in active
  left_normalize_row(s.dsccs);
  for (auto& mscc : s.msccs)
    left_normalize_row(mscc);
}

//integrate states that needed to switch SCCs into corresponding buckets
//create new priorities if necessary
//the sets are provided separately as they might be modified for context
void integrate_switchers(DetConf const& dc, DetConfSets const& sts, DetState& s,
    nba_bitset const& switchers, pri_t& cur_fresh) {
  s.nsccs     |= switchers & sts.nscc_states;
  s.asccs_buf |= switchers & sts.ascc_states;

  nba_bitset const detswitchers = switchers & sts.dscc_states;
  if (detswitchers != 0)
    s.dsccs.push_back(make_pair(detswitchers, cur_fresh++));

  // if (dc.sep_mix) { //if MSCCs separately, add new tree root for switchers
    for (auto const i : ranges::view::ints(0, (int)s.msccs.size())) {
      nba_bitset const mswitchers = switchers & sts.msccs_states[i];
      if (mswitchers != 0)
        s.msccs[i].push_back(make_pair(mswitchers, cur_fresh++));
    }
    /*
  } else {
    //if MSCC root was not present (due to cleanup e.g.)
    nba_bitset const mswitchers = switchers & sts.msccs_states.at(0);
    if (s.msccs.at(0).empty()) {
      if (mswitchers != 0)
        s.msccs.at(0).push_back(make_pair(mswitchers, cur_fresh++));
    } else { //just merge into existing root
      s.msccs.at(0).back().first |= mswitchers;
    }
  }
  */
}


// remove unnecessary empty sets (in dscc and mscc)
void cleanup_empty(DetState &s) {
  auto const is_empty = [](auto const& it){ return it.first == 0; };
  s.dsccs.erase(std::remove_if(begin(s.dsccs), end(s.dsccs), is_empty), end(s.dsccs));
  for (auto& mscc : s.msccs)
    mscc.erase(std::remove_if(begin(mscc), end(mscc), is_empty), end(mscc));
}

// normalize priorities (0..<=n consecutively)
void normalize_prios(DetState &s) {
  //collect
  vector<pri_t> used_pris;
  used_pris.push_back(s.asccs_pri);
  for (auto const& it : s.dsccs)
    used_pris.push_back(it.second);
  for (auto const& mscc : s.msccs)
    for (auto const& it : mscc)
      used_pris.push_back(it.second);

  //get new numbering
  ranges::sort(used_pris);
  map<pri_t, pri_t> f;
  for (auto const i : ranges::view::ints(0, (int)used_pris.size())) {
    f[used_pris[i]] = i;
  }

  //apply
  auto const update_pri = [&f](pri_t& old){ old = f[old]; }; //change prio inplace
  update_pri(s.asccs_pri);
  for (auto& it : s.dsccs)
    update_pri(it.second);
  for (auto& mscc : s.msccs)
    for (auto& it : mscc)
      update_pri(it.second);
}

//convert a rank + event to corresponding fired prio
inline pri_t rank_to_prio(pri_t r, bool good) {
  return 2*(r+1)-(good ? 0 : 1);
}

inline pair<pri_t, bool> prio_to_event(pri_t pri) {
  return make_pair((pri+1)/2 - 1, pri%2==0);

}

//perform breakpoints, detect saturation/death etc
pri_t perform_actions(DetConf const& dc, DetConfSets const& sts,
    DetState const& old, DetState &s, pri_t& cur_fresh) {

  pri_t fired = 2*max_nba_states+1;
  auto const fire = [&fired,&cur_fresh](pri_t& p, bool good){
    fired = min(fired, rank_to_prio(p, good));
    if (!good) //kill rank if it was bad
      p = cur_fresh++;
  };

  //ASCC handling -- fire priority and also swap/cycle set
  int offset = -1; //last/current active ASCC
  if (dc.sep_acc_cyc) {
    for (auto const i : ranges::view::ints(0, (int)sts.asccs_states.size())) {
      if ((old.asccs & sts.asccs_states.at(i))!=0) { //last active SCC found?
        offset = i;
        break;
      }
    }

    //remove states that are in asccs but left the current ascc back into ascc buffer
    if (offset > -1) {
      s.asccs_buf |= s.asccs & ~sts.asccs_states.at(offset);
      s.asccs     &= s.asccs &  sts.asccs_states.at(offset);
    }
  }

  if (s.asccs != 0) {
    if (dc.debug)
      cerr << "ASCC alive" << endl;

    fire(s.asccs_pri, true);
  } else { //breakpoint - either all or the current ASCC died out
    if (dc.debug)
      cerr << "ASCC breaks" << endl;
    fire(s.asccs_pri, false);

    if (!dc.sep_acc_cyc) { //regular breakpoint (swap sets)
      swap(s.asccs, s.asccs_buf);

    } else { //SCC rotating (cyclic) breakpoint
      //TODO: does this harmonize with context pseudo-set?

      if (dc.debug)
        cerr << "offset: " << offset << endl;
      //move next in order from buffer to active
      for (auto const i : ranges::view::ints(0, (int)sts.asccs_states.size())) {
        nba_bitset const cand = s.asccs_buf & sts.asccs_states.at((1+offset+i) % sts.asccs_states.size());
        if (cand != 0) { //next non-empty successor SCC found
          s.asccs_buf &= ~cand; //remove from buffer
          s.asccs = cand; //add to active
          break;
        }
      }
    }
  }

  if (dc.debug)
    cerr << "scanning DSCCs";

  //DSCC handling -- get oldest active
  for (auto& it : s.dsccs) {
    if (it.first == 0) { //set died
      fire(it.second, false);
    } else if ((it.first & dc.aut_acc) != 0) { //set has acc. state
      fire(it.second, true);
    }
  }

  if (dc.debug)
    cerr << " scanning MSCCs" << endl;

  //MSCC handling -- get oldest active, update if MS/Safra
  for (auto& mscc : s.msccs) {
    if (mscc.size() == 0)
      continue; //nothing to do

    //calculate tree stuff
    vector<int> p;
    vector<int> l;
    tie(p,l) = unflatten(mscc);
    if (dc.debug) {
      cerr << "computed tree" << endl;
      cerr << "P: " << seq_to_str(p) << endl;
      cerr << "L: " << seq_to_str(l) << endl;
    }

    // check emptiness, saturation = empty & union of children not empty
    // using the fact that children come before parents we can just run left to right
    vector<bool> node_empty(mscc.size(), true);
    vector<bool> node_saturated(mscc.size(), false);
    vector<int> rightmost_ne_child(mscc.size(), -1); //for müller schupp update

    for (auto const i : ranges::view::ints(0, (int)mscc.size())) {
      if (dc.debug)
        cerr << i << ": ";

      auto const hostempty = mscc[i].first == 0;
      if ((!hostempty || !node_empty[i]) && p[i]!=-1) //me or a child non-empty -> parent non-empty
        node_empty[p[i]] = false;
      if (hostempty && !node_empty[i])  //me empty, children non-empty -> saturated
        node_saturated[i] = true;

      if (hostempty) { //some good or bad event
        if (!node_saturated[i]) { //dead node, kill rank
          fire(mscc[i].second, false);
          if (dc.debug)
            cerr << "dead ";
        } else { //saturated
          fire(mscc[i].second, true);
          if (dc.debug)
            cerr << "strd ";

          if (dc.update == UpdateMode::MUELLERSCHUPP) {
            auto const rc=rightmost_ne_child[i];
            // cerr << "rc: " << rc << " ";
            swap(mscc[i].first, mscc[rc].first);
            // cerr << "merged child" << endl;

          } else if (dc.update == UpdateMode::SAFRA) {
            //collect states of subtree
            nba_bitset subtree = 0;
            for (auto j=l[i]+1; j<i; j++) {
              subtree |= mscc[j].first;
              mscc[j].first = 0;
            }
            if (!dc.puretrees) {
              mscc[i].first = subtree;
            } else { //keep acc/non-acc separated
              mscc[i-1].first = subtree &  dc.aut_acc;
              mscc[i].first   = subtree & ~dc.aut_acc;
            }
          }
        }
      } else {
        if (dc.debug)
          cerr << "neut ";
      }

      //track rightmost nonempty child for the parent
      if (mscc[i].first!=0 && p[i]!= -1)
        rightmost_ne_child[p[i]] = i;

    }
    if (dc.debug)
      cerr << endl;

    //TODO: full collapse for DSCCs and MSCCs, treat all DSCCs separately
  /*
  } else {
    unsigned int active_rk = 2*realn;
    // unsigned int active_ix = 2*realn;

    // find most important rank
    for (auto i=0; i<2*realn; i++) {
      auto const hostempty = suclvl.tups[i].empty();
        if (debug)
          cerr << i << " with rank " << *tupranks[i] << endl;
      if (hostempty && *tupranks[i] < active_rk) {
        if (debug)
          cerr << "empty " << i << " with better rank " << *tupranks[i] << endl;
        // active_ix = i;
        active_rk = *tupranks[i];
      }
    }
    if (debug)
      cerr << "active rank: " << active_rk << endl;

    // merge regions
    if (active_rk < (unsigned int)2*realn) {
      for (auto i=0; i<2*realn-1; i++) {
        if (*tupranks[i] > active_rk && *tupranks[i+1] >= active_rk) { //push stuff through region
          if (*tupranks[i] < *tupranks[i+1]) //push smaller rank forward
            swap(tupranks[i], tupranks[i+1]);
          rord.kill(tupranks[i]); //kill bigger rank
          //collect states
          copy(begin(suclvl.tups[i]), end(suclvl.tups[i]), back_inserter(suclvl.tups[i+1]));
          suclvl.tups[i].clear();
        } else if (*tupranks[i] == active_rk || (*tupranks[i]>active_rk && *tupranks[i+1] < active_rk)) { //right border of region
            if (suclvl.tups[i].empty()) {
              suclvl.prio = min(suclvl.prio, rank_to_prio(*tupranks[i], false));
              rord.kill(tupranks[i]);
            } else {
              suclvl.prio = min(suclvl.prio, rank_to_prio(*tupranks[i], true));
              sort(begin(suclvl.tups[i]), end(suclvl.tups[i]));
            }
        }
      }
      sort(begin(suclvl.tups[2*realn-1]), end(suclvl.tups[2*realn-1]));
      if (active_rk == *tupranks[2*realn-1])
        suclvl.prio = min(suclvl.prio, rank_to_prio(*tupranks[2*realn-1], true));
    }
  }
  */
  }

  return fired;
}

pair<DetState, pri_t> DetState::succ(DetConf const& dc, sym_t x) const {
  bool const& debug = dc.debug;
  if (debug) {
    cerr << "begin " << (int)x << " succ of: " << *this << endl;
  }

  DetState ret(*this); //clone current state
  pri_t cur_fresh = 2*max_nba_states+1; //some for sure unused rank

  successorize_all(dc, ret, x); //calculate successors component-wise

  if (ret.powerset == 0) { //empty successor?
    return {};
  }
  if ((ret.powerset & dc.aut_asinks) != 0) { //check for acc sink reach for early completion
    DetState sink(dc, dc.aut_asinks);
    return make_pair(sink, 0); //good priority fired, sink reached
  }

  //modify sets applying current context, if one provided
  unique_ptr<DetConfSets> modsets = nullptr;
  if (!dc.ctx.empty()) {
    modsets = make_unique<DetConfSets>(dc.sets);
    auto const& curctx = dc.ctx.at(ret.powerset);
    //add relative N/A states
    if (dc.sep_rej)
      modsets->nscc_states = dc.sets.nscc_states | curctx.second;
    if (dc.sep_acc) {
      modsets->ascc_states = dc.sets.ascc_states | curctx.first;
      if (dc.sep_acc_cyc)
        modsets->asccs_states.push_back(~dc.sets.ascc_states & curctx.first);
      else
        modsets->asccs_states.front() |= curctx.first;
    }
    if (debug) {
      cerr << "RA: " << pretty_bitset(modsets->ascc_states);
      cerr << " RN: " << pretty_bitset(modsets->nscc_states) << endl;
    }

    //remove from others
    nba_bitset const tmp = ~(modsets->ascc_states | modsets->nscc_states) & dc.aut_states;
    modsets->dscc_states &= tmp;
    for (auto& mscc : modsets->msccs_states)
      mscc &= tmp;
  }
  auto const& cursets = !modsets ? dc.sets : *modsets; //these will be used for switchers

  left_normalize(ret);
  if (debug) {
    cerr << "normsucc: " << ret << endl;
  }
  expand_trees(dc, ret, cur_fresh);
  if (debug) {
    cerr << "expand: " << ret << endl;
  }
  nba_bitset const switchers = extract_switchers(dc, cursets, ret);
  if (debug) {
    cerr << "-switchers: " << pretty_bitset(switchers) << endl;
    cerr << "half-step: " << ret << endl;
  }

  // half-transition done. now check saturation stuff, get best active, kill ranks...
  pri_t const active_pri = perform_actions(dc, cursets, *this, ret, cur_fresh);
  if (debug) {
    cerr << "merges: " << ret << endl;
  }

  integrate_switchers(dc, cursets, ret, switchers, cur_fresh);
  if (debug) {
    cerr << "+switchers: " << ret << endl;
  }

  // finalize by cleaning empty + fixing priorities
  cleanup_empty(ret);
  normalize_prios(ret);
  if (debug) {
    cerr << "clean-up: " << ret << " -> " << active_pri << endl;
  }
  if (debug) {
    cerr << "----" << endl;
  }
  return make_pair(ret, active_pri);
}

/*
DetState DetState::succ(DetState const& dc, sym_t x) const {

  // ----------------------------------------------------------------

  //calculate tree structure
  vector<int> p;
  vector<int> l;
  tie(p,l) = unflatten(suclvl.tupo);

  if (debug) {
    cerr << "calculated parent + left sibling:" << endl;
    for (auto i=0; i<(int)p.size(); i++) {
      cerr << "p[" << i <<"] =" << p[i] << ", ";
      cerr << "l[" << i <<"] =" << l[i] << endl;
    }
  }

  // ----------------------------------------------------------------
*/

}  // namespace nbautils
