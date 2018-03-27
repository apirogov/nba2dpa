#include <iostream>
#include "io.hh"
#include "common/algo.hh"
#include "common/part_refinement.hh"
#include "pa.hh"

using namespace std;
using namespace nbautils;

// arbitrary sequence (with iterators) to string, intercalated with separator
std::string seq_to_str(vector<string> const& s, std::string const& sep=",") {
  std::stringstream ss;
  auto last = --std::end(s);
  for (auto it = std::cbegin(s); it!=std::cend(s); ++it) {
    ss << *it;
    if (it != last) {
      ss << sep;
    }
  }
  return ss.str();
}

enum class ConstrType { C_FALSE, C_TRUE, C_ATOM, C_AND, C_OR };
struct Constr {
  ConstrType t;
  pair<state_t, state_t> v;
  vector<Constr> d;

  Constr(ConstrType t_) : t(t_) {}
  Constr(state_t a, state_t b) : t(ConstrType::C_ATOM), v(make_pair(min(a,b),max(a,b))) {}
  Constr(ConstrType t_, vector<Constr> const& dat) : t(t_), d(dat) {}

  string to_string() const {
    vector<string> parts;
    for (auto const& it : d)
      parts.push_back(it.to_string());

    if (t==ConstrType::C_TRUE)
      return "1";
    if (t==ConstrType::C_FALSE)
      return "0";
    if (t==ConstrType::C_ATOM)
      return std::to_string(v.first)+","+std::to_string(v.second);
    if (t==ConstrType::C_AND)
      return seq_to_str(parts," & ");
    if (t==ConstrType::C_OR)
      return seq_to_str(parts," | ");
  }
};

// TODO: optimize (prevent recalc.)
// problem: this yields (hopefully) states with finitely witnessed prio preserving, but does not say
// anything about whether states can be merged without harm!!
bool equiv(auto const& aut, auto const& scci, pair<state_t, state_t> start, set<pair<pair<state_t,int>,pair<state_t,int>>>& vis,
    pair<pair<state_t,int>,pair<state_t,int>> cur) {
  state_t a; state_t b; int pa; int pb;
  tie(a,pa) = cur.first;
  tie(b,pb) = cur.second;

  int ca = aut.get_accs(a).front();
  int cb = aut.get_accs(b).front();
  // cerr << "check " << a << "(" << ca << "), " << b << "(" << cb << ") : " << pa << "," << pb << endl;

  //good or bad cycle. if it is priority preserving and potentially another equivalence pair, assume it's good
  if (contains(vis, cur)) {
    //a bad cycle that can not be resolved (i.e. we reached an incompatible pair again or
    //did not preserve the priorities)
    if ((ca!=cb || pa!=pb))
      return false;

    // successful self-reduction (preserved pris and reached a combination of the current pair)
    if (pa == pb && (a == start.first || a==start.second) && (b==start.first || b==start.second))
      return true;
  }

  //if we haven't marked this state, mark to detect cycles
  vis.emplace(cur);

  // if (a==start.first || a==start.second || b==start.first || b==start.second)
  //   if (aut.get_accs(a).front()!=aut.get_accs(b).front())
  //     return false;

  if (scci.scc_of.at(a) != scci.scc_of.at(b)) //only applied within same SCC!
    return false;

  //pri preserved + converged -> witness for future
  if (a==b && pa==pb)
    return true;

  //pri preserved, compatible pair of states other than the current investigated
  //-> might also be a witness pair?
  if (pa==pb && ca==cb && ((a!=start.first && a!=start.second) || (b!=start.first && b!=start.second))) {
    // cerr << "recur" << endl;
    if (equiv(aut, scci, make_pair(a,b), vis, make_pair(make_pair(a,ca),make_pair(b,cb))))
      return true;
  }

  //check for all successor words that there is a witness
  for (int i=0; i<aut.num_syms(); i++) {
    auto suca = aut.succ(a, i).front();
    auto sucb = aut.succ(b, i).front();
    auto sucpa = (int)aut.get_accs(suca).front();
    auto sucpb = (int)aut.get_accs(sucb).front();

    bool ret = equiv(aut, scci, start, vis, make_pair(make_pair(suca, min(pa, sucpa)),
                                            make_pair(sucb, min(pb, sucpb))));
    if (!ret)
      return false;
    // cerr << "-> " << suca << "," << sucb << " " << sucpa << ":" << sucpb << endl;
  }

  return true;
}

//candidates. TODO: need to verify their interdependencies. construct merge graph (a -> b
//if a can be redirected into b <=> b does not visit a other than on converge?)
vector<vector<state_t>> pa_equiv(auto const& aut, auto const& scci) {
  auto sts = aut.states();
  map<state_t, state_t> eq;
  for (auto const st : sts)
    eq[st] = st;

  size_t expl = 0;

  // auto icycp = incomp_cycle_pairs(aut);
  // for (auto p : icycp) {
  //   cout << "incomp " << p.first << "," << p.second << endl;
  // }

  for (int i=0; i<(int)sts.size(); ++i) {
    for (int j=i+1; j<(int)sts.size(); ++j) {
      if (aut.get_accs(sts[i]) != aut.get_accs(sts[j]))
        continue; //different priorities can not be color-equivalent

      // if (scci.scc_of.at(sts[i]) != scci.scc_of.at(sts[j])) //only works within same SCC!
      //   continue;

      // if (contains(icycp, make_pair(sts[i], sts[j])))
      //   continue; //no need to check again

      int pi = aut.get_accs(sts[i]).front();
      int pj = aut.get_accs(sts[j]).front();
      // int p = aut.get_accsets().back();
      // cerr << sts[i] << "," << sts[j] << " " << pi << ":" << pj << endl;
      // cerr << sts[i] << "," << sts[j] << endl;
      set<pair<pair<state_t,int>,pair<state_t,int>>> vis;
      // cerr << "---" << endl;
      if (equiv(aut, scci, make_pair(sts[i],sts[j]), vis , make_pair(make_pair(sts[i], pi), make_pair(sts[j], pj)))) {
        eq[sts[j]] = eq[sts[i]];
        // cerr << "eq " << sts[i] << "," << sts[j] << endl;
      } else {
        // icycp.emplace(make_pair(sts[i], sts[j]));
      }
      expl += vis.size();
    }
  }
  cerr << "explored: " << expl << endl;

  stable_sort(begin(sts), end(sts), [&](state_t a, state_t b){ return eq.at(a)<eq.at(b); });
  return group_by(sts, [&](state_t a, state_t b){ return eq.at(a)==eq.at(b); });
}

// given set in sorted(!!) vector, merge states into one rep. merged mustnot be initial
// after merge graph not normalized
void det_merge_states(auto& aut, vector<state_t> const& others, state_t rep) {
  // cerr << "merging " << seq_to_str(others) << " into " << rep << endl;
  if (others.empty()) return; //nothing to do

  // add redirecting edges into all class members to representative
  for (auto const st : aut.states()) {
      if (st == rep || sorted_contains(others, st))
        continue;

    for (auto const sym : aut.outsyms(st)) {
      auto const suc = aut.succ(st, sym).front();
      if (sorted_contains(others, suc)) {
        aut.set_succs(st, sym, {rep});

        // cerr << st  << " - " << sym << " -> " << suc << " to " << rep << endl;
      }
    }
  }

  //redirect internal edges
  for (auto const sym : aut.outsyms(rep)) {
      auto const suc = aut.succ(rep, sym).front();
      if (!sorted_contains(others, suc))
        continue;

      auto osuc = aut.succ(suc, sym).front();
      if (sorted_contains(others, osuc)) {
        osuc = rep;
      }
      aut.set_succs(rep, sym, {osuc});

      // cerr << rep  << " - " << sym << " -> " << suc << " to " << osuc << endl;
  }

  // finally kill the others!!
  aut.remove_states(others);
}

void det_quotient(auto& aut, vector<vector<state_t>> const& equiv) {
  assert(aut.get_init().size()==1);
  //assert setvec, disjoint, etc

  auto initial = aut.get_init().front();
  bool seenini = false;
  for (auto ecl : equiv) {
    auto rep = ecl.back();
    if (!seenini) {
      auto it = lower_bound(begin(ecl), end(ecl), initial);
      if (it != end(ecl) && *it == initial) {
        ecl.erase(it);
        rep = initial;
        seenini = true;
      } else {
        ecl.pop_back();
      }
    } else {
      ecl.pop_back();
    }

    det_merge_states(aut, ecl, rep);
  }
}

//this is for trying stuff out. do whatever you want here.
int main(int argc, char *argv[]) {
  (void)argv[0];
  (void)argc;
  auto auts = parse_hoa("");

  if (auts.empty())
    cout << "something went wrong!" << endl;

  auto &aut = *auts[0];
  cerr << "autsz: " << aut.num_states() << endl;
  auto aut_st = aut.states();
  succ_fun<state_t> const aut_sucs = [&aut](state_t v){ return aut.succ(v); };
  auto scci = get_sccs(aut_st, aut_sucs);

  // cerr << "is empty? " << pa_is_empty(aut) << endl;
  // cerr << "good set: " << seq_to_str(find_acc_pa_scc(aut)) << endl;
  // vector<state_t> pref;
  // vector<state_t> cyc;
  // tie(pref, cyc) = get_acc_pa_run(aut);
  // cerr << seq_to_str(pref) << " | " << seq_to_str(cyc) << endl;

  // auto aut_prod = pa_union(aut1, aut2);
  // minimize_priorities(*aut_prod);
  // minimize_pa(*aut_prod);
  // aut_prod->tag_to_str = [](auto const& t){ return ""; };
  // print_hoa(*aut_prod);

  auto eqv = pa_equiv(aut, *scci);
  for (auto eq : eqv)
    if (eq.size() > 1)
      cerr << seq_to_str(eq) << endl;
  det_quotient(aut, eqv);
  aut.normalize();
  print_hoa(aut);
}
