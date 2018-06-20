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

std::ostream& operator<<(std::ostream& os, DetConf const& dc) {
  os << "DetConf {" << endl;
  os << "aut_mat: " << !dc.aut_mat.empty() << endl;
  os << "aut_states: " << pretty_bitset(dc.aut_states) << endl;
  os << "aut_acc: "    << pretty_bitset(dc.aut_acc) << endl;
  os << "aut_asinks: " << pretty_bitset(dc.aut_asinks) << endl;
  os << "ctx: "        << !dc.ctx.empty() << endl;

  os << "nscc_states: " <<  pretty_bitset(dc.nscc_states) << endl;
  os << "ascc_states: " <<  pretty_bitset(dc.ascc_states) << endl;
  for (auto const& as : dc.asccs_states)
    os << "\tascc: " <<  pretty_bitset(as) << endl;
  os << "dscc_states: " <<  pretty_bitset(dc.dscc_states) << endl;
  os << "msccs:";
  for (auto const& ms : dc.msccs_states)
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

pair<vector<int>,vector<int>> unflatten(vector<pri_t> const& rank) {
  int const n = rank.size();
  pair<vector<int>,vector<int>> ret = make_pair(vector<int>(n),vector<int>(n));
  stack<int> s;
  ret.first[n-1] = n-1;
  ret.second[n-1] = -1;
  s.push(n-1);
  for (int i=n-2; i>=0; i--) {
    while (rank[i] < rank[s.top()]) {
      ret.second[s.top()] = i;
      s.pop();
    }
    ret.first[i] = s.top();
    s.push(i);
  }
  while (!s.empty()) {
    ret.second[s.top()] = -1;
    s.pop();
  }

  return ret;
}

DetState::DetState() {}

//put everything into one set... in later steps it will differentiate
DetState::DetState(DetConf const& dc, nba_bitset const& qs) {
  msccs.resize(dc.msccs_states.size());
  msccs.at(0) = {make_pair(qs, 1)};
  powerset = qs;
}

//convert a rank + event to corresponding fired prio
inline pri_t rank_to_prio(pri_t r, bool good) {
  return 2*(r+1)-(good ? 0 : 1);
}

//apply successor set function on each set separately, inplace
void successorize_all(DetConf const& dc, DetState& s, sym_t const x) {
  auto const psucc = [&dc,x](auto const& bset){
    return powersucc(dc.aut_mat, bset, x, dc.aut_asinks); };

  s.powerset  = psucc(s.powerset);
  s.nsccs     = psucc(s.nsccs);
  s.asccs_buf = psucc(s.asccs_buf);
  s.asccs     = psucc(s.asccs);
  for (auto i : ranges::view::ints(0, (int)s.dsccs.size()))
    s.dsccs[i].first = psucc(s.dsccs[i].first);
  for (auto i : ranges::view::ints(0, (int)s.msccs.size()))
    for (auto j : ranges::view::ints(0, (int)s.msccs[i].size()))
      s.msccs[i][j].first = psucc(s.msccs[i][j].first);
}

//remove wrong located states, return their collection
nba_bitset extract_switchers(DetConf const& dc, DetState& s) {
  nba_bitset switchers;
  switchers |= s.nsccs     & (~dc.nscc_states & dc.aut_states);
  switchers |= s.asccs     & (~dc.ascc_states & dc.aut_states);
  switchers |= s.asccs_buf & (~dc.ascc_states & dc.aut_states);
  s.nsccs     &= dc.nscc_states;
  s.asccs     &= dc.ascc_states;
  s.asccs_buf &= dc.ascc_states;
  for (auto i : ranges::view::ints(0, (int)s.dsccs.size())) {
    switchers |= s.dsccs[i].first & (~dc.dscc_states & dc.aut_states);
    s.dsccs[i].first &= dc.dscc_states;
  }
  for (auto i : ranges::view::ints(0, (int)s.msccs.size()))
    for (auto j : ranges::view::ints(0, (int)s.msccs[i].size())) {
      switchers |= s.msccs[i][j].first & (~dc.msccs_states[i] & dc.aut_states);
      s.msccs[i][j].first &= dc.msccs_states[i];
    }
  return switchers;
}

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
void integrate_switchers(DetConf const& dc, DetState& s,
    nba_bitset const& switchers, pri_t& cur_fresh) {
  s.nsccs     |= switchers & dc.nscc_states;
  s.asccs_buf |= switchers & dc.ascc_states;

  nba_bitset const detswitchers = switchers & dc.dscc_states;
  if (detswitchers != 0)
    s.dsccs.push_back(make_pair(detswitchers, cur_fresh++));

  if (dc.sep_mix) { //if MSCCs separately, add new tree root for switchers
    for (auto i : ranges::view::ints(0, (int)s.msccs.size())) {
      nba_bitset const mswitchers = switchers & dc.msccs_states[i];
      if (mswitchers != 0)
        s.msccs[i].push_back(make_pair(mswitchers, cur_fresh++));
    }
  } else {
    //if MSCC root was not present (due to cleanup e.g.)
    nba_bitset const mswitchers = switchers & dc.msccs_states[0];
    if (s.msccs.empty() && mswitchers != 0)
      s.msccs.push_back({make_pair(mswitchers, cur_fresh++)});
    else //just merge into existing root
      s.msccs[0].back().first |= mswitchers;
  }
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
  for (auto i : ranges::view::ints(0, (int)used_pris.size())) {
    f[used_pris[i]] = i;
  }
  auto const update_pri = [&f](pri_t& old){ old = f[old]; }; //change prio inplace

  //apply
  update_pri(s.asccs_pri);
  for (auto& it : s.dsccs)
    update_pri(it.second);
  for (auto& mscc : s.msccs)
    for (auto& it : mscc)
      update_pri(it.second);
}

pair<DetState, pri_t> DetState::succ(DetConf const& dc, sym_t x) const {
  bool const& debug = dc.debug;
  // bool const& debug = true;
  if (debug) {
    cerr << "begin succ of: " << *this << endl;
  }

  DetState ret(*this); //clone current state
  pri_t cur_fresh = 2*max_nba_states+1; //some for sure unused rank

  successorize_all(dc, ret, x);

  //check for acc sink reach
  if ((ret.powerset & dc.aut_asinks) != 0) {
    DetState sink(dc, dc.aut_asinks);
    return make_pair(sink, 0); //good priority fired, sink reached
  }

  // -- in any order --
  left_normalize(ret);
  // cerr << pretty_bitset(s.powerset) << endl;
  nba_bitset const switchers = extract_switchers(dc, ret);
  expand_trees(dc, ret, cur_fresh);
  // -- end of in any order --
  cerr << "extracted switchers: " << pretty_bitset(switchers) << endl;
  cerr << "half-step: " << ret << endl;

  // half-transition done. now check saturation stuff, get best active, kill ranks...
  pri_t const active_pri = 0; //= get_vip_active(dc, ret) ; //some default worst case priority

  //TODO: perform breakpoints, detect saturation/death etc

  // finalize by reintegrating states that switched and fixing priorities
  integrate_switchers(dc, ret, switchers, cur_fresh);
  cerr << "integrated switchers: " << ret << endl;
  cleanup_empty(ret);
  cerr << "cleaned up: " << ret << endl;
  normalize_prios(ret);
  return make_pair(ret, active_pri);
}

/*
DetState DetState::succ(DetState const& dc, sym_t x) const {
  //successor helper
  auto xsucc = [&lvc,&x](auto const& ps){ return powersucc(*lvc.aut, ps, x); };
  auto sucpset = xsucc(states()); //successor powerset (without structure)

  //check for accepting sinks (then we're happily done)
  if (!set_intersect_empty(sucpset,lvc.accsinks)) {
    if (debug) {
      cerr << "reached accepting sink. returning accepting sink level!" << endl;
    }
    auto ret = Level(lvc,lvc.accsinks);
    ret.prio = 0; //must accept
    return ret;
  }

  //define helpers to detect whether state is in (relative) acc/rej SCC
  auto is_xscc = [&lvc,&sucpset](auto const& as, auto const &cs, state_t s){
    //state is in xscc of original automaton?
    bool xscc = contains(as, lvc.aut_scc->scc_of.at(s));
    //state is in relative xscc of current context (if given)?
    bool cxscc = false;
    if (lvc.ctx) {
      auto stt = make_pair(sucpset, s);
      if (lvc.ctx->tag->has(stt)) {
        auto st = lvc.ctx->tag->get(stt);
        cxscc = contains(cs, lvc.ctx_scc->scc_of.at(st));
      }
    }
    return xscc || cxscc;
  };
  auto is_nscc = [&lvc,&is_xscc](state_t s){ return is_xscc(lvc.aut_cl->rejecting, lvc.ctx_cl->rejecting, s); };
  auto is_ascc = [&lvc,&is_xscc](state_t s){ return is_xscc(lvc.aut_cl->accepting, lvc.ctx_cl->accepting, s); };

  //helper that says whether a state is accepting in original NBA
  auto is_acc = [&lvc](state_t s){ return lvc.aut->has_accs(s); };

  // define left-reduce helper
  set<Level::state_t> used_sucs; //track already used successors
  //left-reduce: tmp = suctups[i] \ used_sucs; used_sucs = used_sucs U tmp
  auto leftreduce = [&](vector<Level::state_t> const& v){
      vector<Level::state_t> tmp;
      set_difference(begin(v), end(v),
                     begin(used_sucs), end(used_sucs),
                     back_inserter(tmp));
      copy(begin(tmp), end(tmp), inserter(used_sucs, end(used_sucs)));
      return tmp;
  };

  // ----------------------------------------------------------------
  // calculate successor sets
  vector<vector<Level::state_t>> suctups;
  transform(begin(tups),end(tups),back_inserter(suctups), xsucc);

  int oldaccscc = -1; //no known previous scc
  if (lvc.sep_acc_cyc) {
    //for cyclic ASCC probing we need to know the current one, if any
    if (!tups.back().empty())
      oldaccscc = lvc.aut_scc->scc_of.at(tups.back().front());
  }

  //split out accepting scc successors for breakpoint construction sets
  vector<vector<Level::state_t>> sascc;
  if (lvc.sep_acc) {
    copy(end(suctups)-2, end(suctups), back_inserter(sascc));
    suctups.erase(end(suctups)-2, end(suctups));
  }

  if (debug) {
    tmplv.tupo = tupo; tmplv.tups = suctups;
    copy(begin(sascc), end(sascc), back_inserter(tmplv.tups));
    cerr << "after powersucc: " << tmplv.to_string() << endl;
  }

  // ----------------------------------------------------------------
  // first take care of breakpoint sets

  // what is in good accepting scc stays there, rest goes back into root
  if (lvc.sep_acc) {
    if (suctups.empty()) { //ensure that there is a safra root (edge case!)
      suctups.push_back({});
    }

    for (int i=1; i>=0; i--) {
      //get unused successors
      vector<Level::state_t> tmp = leftreduce(sascc[i]);

      //separate the ones in ASCCs and the ones not
      auto sep = stable_partition(begin(tmp), end(tmp), is_ascc);
      sascc[i] = vector<Level::state_t>(begin(tmp), sep);
      auto out = vector<Level::state_t>(sep, end(tmp));

      //merge the ones not into safra root and allow to use them
      tmp.clear();
      set_union(begin(suctups.back()),end(suctups.back()), begin(out), end(out), back_inserter(tmp));
      swap(suctups.back(), tmp);

      tmp.clear();
      set_difference(begin(used_sucs),end(used_sucs), begin(out), end(out), back_inserter(tmp));
      used_sucs = set<state_t>(begin(tmp), end(tmp));
    }

    if (debug) {
      tmplv.tupo = tupo; tmplv.tups = suctups;
      copy(begin(sascc), end(sascc), back_inserter(tmplv.tups));
      cerr << "reduced ASCCs:" << tmplv.to_string() << endl;
    }
  }

  // ----------------------------------------------------------------
  // perform partial safra tree update

  Level suclvl;

  auto oldtupo = tupo;
  if (lvc.sep_acc && tupo.size()==1) { //edge case if safra tree was dead
    oldtupo = {1, 0};
  }

  int n = tups.size(); //number of sets in tuple
  int realn = suctups.size(); //number of sets in actual safra tree

  //allocate twice as many tuples for next
  suclvl.tups = vector<vector<state_t>>(2*realn, vector<state_t>());
  suclvl.tupo = vector<ord_t>(2*realn, 0);

  // prioritize left, split F/not F, forward token to right
  for (auto i=0; i<realn; i++) {
    vector<Level::state_t> tmp = leftreduce(suctups[i]);

    //separate F / not F successors: (l,r) = (Δ_F, Δ_notF)
    auto sep = stable_partition(begin(tmp), end(tmp), is_acc);
    copy(begin(tmp), sep,      back_inserter(suclvl.tups[2*i]));
    copy(sep,        end(tmp), back_inserter(suclvl.tups[2*i+1]));

    //pass token -> old one to right, new one to left
    suclvl.tupo[2*i] = oldtupo.size()+i;
    suclvl.tupo[2*i+1] = oldtupo[i];
  }

  if (debug) {
    tmplv = suclvl;
    copy(begin(sascc), end(sascc), back_inserter(tmplv.tups));
    if (lvc.sep_acc) tmplv.tupo.push_back(tupo.back());
    cerr << "left-reduce + split:" << tmplv.to_string() << endl;
  }

  // ----------------------------------------------------------------
  // now take care of keeping NSCC-states in root

  if (lvc.sep_rej) {

    set<Level::state_t> nsccst;
    for (auto &t : suclvl.tups) {
      copy_if(begin(t),end(t),inserter(nsccst,end(nsccst)),is_nscc);
      t.erase(remove_if(t.begin(),t.end(),is_nscc), end(t));
    }

    // relocate NSCC states into root
    vector<Level::state_t> tmp;
    set_union(begin(suclvl.tups.back()),end(suclvl.tups.back()),
              begin(nsccst), end(nsccst), back_inserter(tmp));
    swap(suclvl.tups.back(), tmp);

    if (debug) {
      tmplv = suclvl;
      copy(begin(sascc), end(sascc), back_inserter(tmplv.tups));
      if (lvc.sep_acc) tmplv.tupo.push_back(tupo.back());
      cerr << "after relocating NSCC-states into root:" << tmplv.to_string() << endl;
    }
  }

  // ----------------------------------------------------------------
  // relocate remaining ASCC states in safra nodes into brkpt set
  if (lvc.sep_acc) {
    set<Level::state_t> asccst;
    for (auto &t : suclvl.tups) {
      copy_if(begin(t),end(t),inserter(asccst,end(asccst)),is_ascc);
      t.erase(remove_if(t.begin(),t.end(),is_ascc), end(t));
    }

    vector<Level::state_t> tmp;
    auto& brkset = sascc.front();
    set_union(begin(brkset),end(brkset),
              begin(asccst), end(asccst), back_inserter(tmp));
    swap(brkset, tmp);

    if (debug) {
      tmplv = suclvl;
      copy(begin(sascc), end(sascc), back_inserter(tmplv.tups));
      if (lvc.sep_acc) tmplv.tupo.push_back(tupo.back());
      cerr << "after relocating ASCC-states into extra set:" << tmplv.to_string() << endl;
    }
  }

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

  auto ascc_ord = tupo.back(); //save old ascc ord (if sep_acc)
  if (lvc.sep_acc) // add the token of breakpoint construction back into list to get a token
    suclvl.tupo.push_back(tupo.back());

  //create token manager with enough slots
  RelOrder rord(suclvl.tupo.size());
  //get all tokens
  auto tupranks = rord.from_ranks(suclvl.tupo);

  // ----------------------------------------------------------------

  //default: low impact (high) bad value
  suclvl.prio = rank_to_prio(n, false);

  if (debug) {
    cerr << "def prio: " << suclvl.prio << endl;
  }

  if (lvc.sep_acc) {
    bool ascc_empty = sascc.back().empty();
    if (ascc_empty) { //breakpoint?
          rord.kill(tupranks.back());
          suclvl.prio = min(suclvl.prio, rank_to_prio(ascc_ord, false));

          if (debug)
            cerr << "reached ASCC breakpoint" << endl;

          if (!lvc.sep_acc_cyc) {
            swap(sascc.front(), sascc.back()); //breakpoint (right was empty) -> just swap sets

            if (debug)
              cerr << "swapped ASCC buffer and active set" << endl;

          } else if (!sascc.front().empty()){ //put next scc in order into the set, if any
            //TODO: cyclic update seems to have a bug still FIXME
            //possible problem/fix: need to put _succs_ of next SCC in order (sort on pred SCC)

            // auto tmp = sascc.front();

            auto tmp = *(end(tups)-2);

            stable_sort(begin(tmp), end(tmp),
                 [&](state_t const& a, state_t const& b){
                 return lvc.aut_scc->scc_of.at(a) < lvc.aut_scc->scc_of.at(b); });
            int nextscc = -1;
            //take first bigger
            for (auto st : tmp) {
              int stscc = lvc.aut_scc->scc_of.at(st);
              if (stscc > oldaccscc) {
                nextscc = stscc;
                break;
              }
            }
            //or start cycle anew
            if (nextscc == -1 && !tmp.empty())
              nextscc = lvc.aut_scc->scc_of.at(tmp.front());

            // get the predecessors with the next SCC
            tmp = states();
            auto it = stable_partition(begin(tmp), end(tmp), [&lvc, &nextscc](auto s){
                return lvc.aut_scc->scc_of.at(s) == (unsigned)nextscc; });
            tmp.erase(it, end(tmp));

            if (nextscc > -1) { //move just states of cyclically next scc into right set

              // separate successors of preds and non-successors
              auto const predsuc = xsucc(tmp);
              sascc.back() = set_intersect(sascc.front(), predsuc);
              sascc.front() = set_diff(sascc.front(), predsuc);

              // auto it = stable_partition(begin(tmp), end(tmp), [&lvc, &nextscc](auto s){
              //     return lvc.aut_scc->scc_of.at(s) != (unsigned)nextscc; });
              // copy(it, end(tmp), back_inserter(sascc.back()));
              // tmp.erase(it, end(tmp));
              // sort(begin(tmp), end(tmp));
              // sascc.front() = tmp;

              if (debug)
                cerr << "as ASCC " << oldaccscc << " died, now it's the turn of states from ASCC " << nextscc << endl;
            } else {
              if (debug)
                cerr << "no states in buffering ASCC set" << endl;
            }

          }
          // cerr << seq_to_str(sascc.front()) << " , " << seq_to_str(sascc.back()) << endl;
    } else {
          suclvl.prio = min(suclvl.prio, rank_to_prio(ascc_ord, true));
    }

    if (debug) {
      tmplv = suclvl;
      copy(begin(sascc), end(sascc), back_inserter(tmplv.tups));
      tmplv.tupo = rord.to_ranks(tupranks);
      cerr << "after updating ASCC set:" << tmplv.to_string() << endl;
    }
  }

  // ----------------------------------------------------------------
  // check emptiness, saturation = empty & union of children not empty
  // using the fact that children come before parents we can just run left to right
  vector<bool> node_empty(2*realn, true);
  vector<bool> node_saturated(2*realn, false);
  vector<int> rightmost_ne_child(2*realn, -1); //for müller schupp update

  if (lvc.update != LevelUpdateMode::FULLMERGE) {

  for (auto i=0; i<2*realn; i++) {
    auto const hostempty = suclvl.tups[i].empty();
    if (!hostempty || !node_empty[i]) //me or children non-empty -> parent non-empty
      node_empty[p[i]] = false;
    if (hostempty && !node_empty[i]) //me empty, children non-empty-> saturated
      node_saturated[i] = true;

    if (hostempty) { //some good or bad event
      if (node_saturated[i]) { //saturated node
        if (debug)
          cerr << i << " strd ";

        suclvl.prio = min(suclvl.prio, rank_to_prio(suclvl.tupo[i], true));

        if (lvc.update == LevelUpdateMode::MUELLERSCHUPP) {
          auto const rc=rightmost_ne_child[i];
          move(begin(suclvl.tups[rc]),
               end(suclvl.tups[rc]),
               back_inserter(suclvl.tups[i]));
          suclvl.tups[rc].clear();
          rord.kill(tupranks[rc]);
        } else if (lvc.update == LevelUpdateMode::SAFRA) {
          for (auto j=l[i]+1; j<i; j++) {
            move(begin(suclvl.tups[j]),
                end(suclvl.tups[j]),
                back_inserter(suclvl.tups[i]));
            suclvl.tups[j].clear();
            if (!lvc.pure || j!=i-1)
              rord.kill(tupranks[j]);
          }
          sort(begin(suclvl.tups[i]), end(suclvl.tups[i]));

          //TODO: find rare bug, make it work with full collapse too
          if (lvc.pure) { //push out accepting back into a leaf
            auto it = stable_partition(begin(suclvl.tups[i]), end(suclvl.tups[i]),
                [&is_acc](state_t s){ return !is_acc(s);});
            if (it != begin(suclvl.tups[i]) && it != end(suclvl.tups[i])) { // if node is mixed
              if (debug)
                cerr << i << " purified ";
              //separate accepting back out into child
              suclvl.tups[i-1] = vector<small_state_t>(it, end(suclvl.tups[i]));
              suclvl.tups[i].erase(it, suclvl.tups[i].end());
            } else {
              rord.kill(tupranks[i-1]);
            }
          }
        }
      } else { //dead node, kill rank
        if (debug)
          cerr << i << " dead ";

        suclvl.prio = min(suclvl.prio, rank_to_prio(suclvl.tupo[i], false));
        rord.kill(tupranks[i]);
      }
    }

    //track rightmost nonempty child for the parent
    if (!suclvl.tups[i].empty() && p[i]!=i)
      rightmost_ne_child[p[i]] = i;
  }
  if (debug)
    cerr << endl;

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

  if (debug) {
    tmplv = suclvl;
    copy(begin(sascc), end(sascc), back_inserter(tmplv.tups));
    tmplv.tupo = rord.to_ranks(tupranks);
    cerr << "after completed merges:" << tmplv.to_string() << endl;
  }

  // ----------------------------------------------------------------
  //cleanup remaining empty sets (all empty sets must be killed)
  auto itt = remove_if(begin(suclvl.tups), end(suclvl.tups),
                      [](auto const& s){return s.empty();});
  suclvl.tups.erase(itt, end(suclvl.tups));

  //keep the non-killed ranks, i.e. 0-num_nonemptysets-1
  size_t num_keep = suclvl.tups.size() + (lvc.sep_acc ? 1 : 0);
  rord.normalize();
  auto ito = stable_partition(begin(tupranks), end(tupranks),
                      [&](auto s){return *s<num_keep;});
  tupranks.erase(ito, end(tupranks));

  if (lvc.sep_acc) {
    // insert new accscc set back
    for (auto const& v : sascc)
      suclvl.tups.push_back(v);
  }

  //get new resulting rank order
  suclvl.tupo = rord.to_ranks(tupranks);

  if (debug) {
    cerr << "after cleanup:" << suclvl.to_string() << endl;
  }

  // ----------------------------------------------------------------
  assert(suclvl.states() == sucpset); //consistent with powerset

  // finalize by assigning powerset stamp
  suclvl.powerset = add_powerset_hash(*lvc.aut, suclvl);

  if (debug) {
    cerr << "completed succ" << endl;
  }

  return suclvl;
}
*/

}  // namespace nbautils
