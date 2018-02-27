#pragma once

namespace nbautils {

enum class PAType { MIN_EVEN, MIN_ODD, MAX_EVEN, MAX_ODD };

inline constexpr PAType opposite_parity(PAType pt) {
  switch (pt) {
    case PAType::MIN_EVEN: return PAType::MIN_ODD;
    case PAType::MIN_ODD:  return PAType::MIN_EVEN;
    case PAType::MAX_EVEN: return PAType::MAX_ODD;
    case PAType::MAX_ODD:  return PAType::MAX_EVEN;
  }
}

inline constexpr PAType opposite_polarity(PAType pt) {
  switch (pt) {
    case PAType::MIN_EVEN: return PAType::MAX_EVEN;
    case PAType::MIN_ODD:  return PAType::MAX_ODD;
    case PAType::MAX_EVEN: return PAType::MIN_EVEN;
    case PAType::MAX_ODD:  return PAType::MIN_ODD;
  }
}

inline constexpr bool pa_acc_is_min(PAType a) {
  return a==PAType::MIN_EVEN || a==PAType::MIN_ODD;
}
inline constexpr bool pa_acc_is_even(PAType a) {
  return a==PAType::MIN_EVEN || a==PAType::MAX_EVEN;
}
inline constexpr bool pa_acc_is_max(PAType a) {
  return !pa_acc_is_min(a);
}
inline constexpr bool pa_acc_is_odd(PAType a) {
  return !pa_acc_is_even(a);
}

inline constexpr bool same_parity(PAType a, PAType b) {
  return (pa_acc_is_even(a) == pa_acc_is_even(b));
}
inline constexpr bool same_parity(int a, int b) {
  return (a%2==0)==(b%2==0);
}
inline constexpr bool same_polarity(PAType a, PAType b) {
  return (pa_acc_is_min(a) == pa_acc_is_min(b));
}
inline constexpr bool good_priority(PAType a, int p) {
  return p%2 == (pa_acc_is_even(a) ? 0 : 1);
}
inline auto stronger_priority_f(PAType a) {
  return pa_acc_is_max(a) ? [](int p, int q){ return max(p,q); }
                          : [](int p, int q){ return min(p,q); };
}

// takes source and target priority acceptance conditions and a list of priorities of the
// source type, returns function mapping source pris to semantically equiv. in target cond.
inline function<int(int)> priority_transformer(PAType from, PAType to, int minpri, int maxpri) {
  function<int(int)> idmap = identity;
  if (from==to) //nothing to do
    return identity;

  function<int(int)> parflip = [minpri](int a){ return a+(minpri%2==0 ? 1 : -1); };
  auto adapt_par = same_parity(from,to)   ? idmap : parflip;

  // auto mintmp = min(adapt_par(minpri), adapt_par(maxpri));
  auto maxtmp = max(adapt_par(minpri), adapt_par(maxpri));
  function<int(int)> polflip = [maxtmp](int a){ return -a + (maxtmp + maxtmp%2); };
  auto adapt_pol = same_polarity(from,to) ? idmap : polflip;

  return [=](int v){ return adapt_pol(adapt_par(v)); };
}

}
