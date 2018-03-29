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


// TODO: optimize (prevent recalc.)
// returns whether a~b and whether b ~> a and a ~> b (a ~> b := a can be safely merged into b)
// bool equiv(auto const& aut, auto const& scci, pair<state_t, state_t> start,
pair<bool,pair<bool,bool>> equiv(auto const& aut, auto const& scci, pair<state_t, state_t> start,
    set<pair<pair<state_t,int>,pair<state_t,int>>>& vis,
    pair<pair<state_t,int>,pair<state_t,int>> cur, int k=0) {
  auto const triv_true = make_pair(true, make_pair(true, true));
  auto const triv_false = make_pair(false, make_pair(false, false));

  state_t a; state_t b; int pa; int pb;
  tie(a,pa) = cur.first;
  tie(b,pb) = cur.second;

  int ca = aut.get_accs(a).front();
  int cb = aut.get_accs(b).front();

  // for (int i=0; i<k; ++i)
    // cerr << " ";
  // cerr << "check " << a << "(" << ca << "), " << b << "(" << cb << ") : " << pa << "," << pb << endl;
  if (vis.size() % 1000 == 0)
    cerr << vis.size() << endl;

  //good or bad cycle. if it is priority preserving and potentially another equivalence pair, assume it's good
  if (contains(vis, cur)) {
    // vis.erase(cur);

    // successful self-reduction (preserved pris and reached a combination of the current candidate pair)
    // NOTE: subsumed by next case... just for conceptual clarity
    if (pa == pb && (a == start.first || a==start.second) && (b==start.first || b==start.second))
      return triv_true;

    // this is a candidate state without violations that we've already seen. we can assume
    // here that it is good. it still can be falsified above
    if (pa == pb && ca == cb)
      return triv_true;

    //a bad cycle that can not be resolved (i.e. we reached an incompatible pair again or
    //did not preserve the priorities)
    if ((ca!=cb || pa!=pb))
      return triv_false;
  }

  //if we haven't marked this state, mark to detect cycles
  vis.emplace(cur);

  //pri preserved + converged -> witness for future (simplest good case)
  if (a==b && pa==pb)
    return triv_true;


  // if (a==start.first || a==start.second || b==start.first || b==start.second)
  //   if (aut.get_accs(a).front()!=aut.get_accs(b).front())
  //     return false;

  // if (scci.scc_of.at(a) != scci.scc_of.at(b)) //only applied within same SCC!
  //   return triv_false;

  //check for all successor words that there is a witness
  auto ret = triv_true;
  // path from left start state over right start state, but not converge
  if (a == start.second && b != start.first && b != start.second)
    ret.second.first = false;
  // path from right start state over left start state, but not converge
  if (b == start.first  && a != start.second && a != start.first)
    ret.second.second = false;

  for (int i=0; i<aut.num_syms(); i++) {
    auto suca = aut.succ(a, i).front();
    auto sucb = aut.succ(b, i).front();
    // cerr << a << " " << i << " " << suca << " , " << b << " " << i << " " << sucb << endl;
    auto sucpa = (int)aut.get_accs(suca).front();
    auto sucpb = (int)aut.get_accs(sucb).front();

    //pri preserved, compatible pair of states other than the current investigated
    //-> might also be a witness pair?
    auto tmp = triv_false;
    if (pa==pb && ca==cb && ((suca!=start.first && suca!=start.second) || (sucb!=start.first && sucb!=start.second))) {
      // cerr << "recur" << endl;
      state_t norm_a = suca; //min(a,b);
      state_t norm_b = sucb; //max(a,b);
      tmp = equiv(aut, scci, make_pair(norm_a,norm_b), vis, make_pair(make_pair(norm_a,sucpa),make_pair(norm_b,sucpb)), k+1);
    }

    if (!tmp.first) {
      tmp = equiv(aut, scci, start, vis, make_pair(make_pair(suca, min(pa, sucpa)),
                                         make_pair(sucb, min(pb, sucpb))), k+1);
    }

    if (!tmp.first)
      return triv_false;

    ret.second.first = ret.second.first && tmp.second.first;
    ret.second.second = ret.second.second && tmp.second.second;
    // cerr << "-> " << suca << "," << sucb << " " << sucpa << ":" << sucpb << endl;
  }

  return ret;
}

//calculate converging states and merge orders
//TODO: for each ~ class calc SCCs of ~>, decompose dag paths, merge accordingly
vector<vector<state_t>> pa_equiv(auto& aut, auto const& scci) {
  auto sts = aut.states();
  map<state_t, state_t> eq;
  for (auto const st : sts)
    eq[st] = st;

  size_t expl = 0;

  // auto icycp = incomp_cycle_pairs(aut);
  // for (auto p : icycp) {
  //   cout << "incomp " << p.first << "," << p.second << endl;
  // }

  // int merges=0;
  for (int i=0; i<(int)sts.size(); ++i) {
    for (int j=i+1; j<(int)sts.size(); ++j) {
      // if (aut.get_init().front() == sts[i])
      //   continue; //leave initial state alone (for now)

      //some state merged already
      // if (!aut.has_state(sts[i]))
      //   break;
      // if (!aut.has_state(sts[j]))
      //   continue;

      if (aut.get_accs(sts[i]) != aut.get_accs(sts[j]))
        continue; //different priorities can not be color-equivalent

      // if (merges>=1)
      //   break;

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
      auto tmp = equiv(aut, scci, make_pair(sts[i],sts[j]), vis,
                                  make_pair(make_pair(sts[i], pi), make_pair(sts[j], pj)));
      cerr << sts[i] << "~" << sts[j] << " | " << tmp.second.first<<","<<tmp.second.second<< endl;
      if (tmp.first && tmp.second.first && tmp.second.second) { //TODO: SCC stuff
        eq[sts[j]] = eq[sts[i]];

        // auto aut2(aut);
        // det_quotient(aut2, {{sts[i], sts[j]}});
        // if (!pa_equivalent(aut, aut2))
        //   cerr << "merge failure" << endl;
        // det_quotient(aut, {{sts[i], sts[j]}});
        // merges++;
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

  // cerr << pa_equivalent(aut, aut) << endl;
  // auto aut2(aut);
  // complement_pa(aut2);
  // cerr << pa_equivalent(aut, aut2) << endl;


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
  // det_quotient(aut, eqv);
  aut.normalize();
  print_hoa(aut);
}
