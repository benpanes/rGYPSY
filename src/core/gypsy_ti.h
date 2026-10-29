#ifndef GYPSY_TI_H
#define GYPSY_TI_H

#include "gypsy_eq.h"

// ---------------------------------------------------------------------------
// Tree improvement: site index adjusted for genetic worth, following
// section 2.4 of
//
//   Luo, D. and Thomas, B.R. 2021. An analysis of age-age correlations in
//   white spruce and lodgepole pine and how it applies to the growth and
//   yield projection system (GYPSY) in Alberta. Forest Ecology and
//   Management 482: 118865. https://doi.org/10.1016/j.foreco.2020.118865
//
//   AgeAgeCor: age-age correlation between selection age and total age 50,
//              from the adjusted Lambeth equations (their Table 2)
//                r = c + m ln(sel_age / 50)
//   SI_TI:     improved SI_t = ((gw / 100) r + 1) SI_t
//   SimTI (gypsy_sim.h) runs a full simulation from the improved SI
//
// gw is genetic worth for height at the selection age (%). The adjustment
// is always to the total age SI (top height at total age 50); a breast
// height SI is converted to and from it.
// Parameters exist for SW and PL only.
// ---------------------------------------------------------------------------

namespace gypsy {

// Adjusted Lambeth equation parameters by species (NaN: none)
//                          AW   SB   SW     PL
const double TI_C[4]   = {NaN, NaN, 0.993, 1.017};
const double TI_M[4]   = {NaN, NaN, 0.255, 0.445};

inline bool has_ti(int s) { return s == SW || s == PL; }

inline double AgeAgeCor(int s, double sel_age) {
  return TI_C[s] + TI_M[s]*std::log(sel_age/50.0);
}

// gw of 0 or NaN leaves si_t unchanged (the vector version and SimTI warn
// about missing SW and PL values)
inline double SI_TI(int s, double si_t, double gw, double sel_age) {
  if (miss(gw) || gw == 0) return si_t;
  return ((gw/100.0)*AgeAgeCor(s, sel_age) + 1.0)*si_t;
}

// Checks the genetic worth inputs and returns gw with missing SW and PL
// values set to 0, with a warning listing those rows. Rows with a genetic
// worth must be SW or PL with a positive selection age; throws listing the
// rows that aren't.
inline std::vector<double> check_ti(const std::vector<int>& spp, const std::vector<double>& gw,
                                    const std::vector<double>& sel_age, Messages& warn) {
  std::vector<double> out = gw;
  std::vector<int> bad_spp, bad_age, no_gw;
  for (std::size_t i = 0; i < spp.size(); i++) {
    if (miss(gw[i]) && has_ti(spp[i])) { out[i] = 0; no_gw.push_back(static_cast<int>(i)); }
    if (miss(gw[i]) || gw[i] == 0) continue;
    if      (!has_ti(spp[i]))     bad_spp.push_back(static_cast<int>(i));
    else if (!(sel_age[i] > 0))   bad_age.push_back(static_cast<int>(i));
  }
  if (!bad_spp.empty())
    throw Error("Genetic worth can only be given for SW and PL. Check row(s): " + join_rows(bad_spp));
  if (!bad_age.empty())
    throw Error("Genetic worth needs a selection age > 0. Check row(s): " + join_rows(bad_age));
  if (!no_gw.empty())
    warn.push_back("No genetic worth for SW / PL row(s) " + join_rows(no_gw) + "; treated as 0 (no adjustment).");
  return out;
}

// Vector version, one value per row
inline std::vector<double> SI_TI(const std::vector<std::string>& spp, const std::vector<double>& si_t,
                                 const std::vector<double>& gw, const std::vector<double>& sel_age,
                                 Messages& warn) {
  std::size_t n = spp.size();
  for (std::size_t m : {si_t.size(), gw.size(), sel_age.size()})
    if (m != n) throw Error("spp, si_t, gw, and sel_age must all be the same length.");
  std::vector<int> s(n), bad;
  for (std::size_t i = 0; i < n; i++) {
    s[i] = spp_index(spp[i]);
    if (s[i] < 0) bad.push_back(static_cast<int>(i));
  }
  if (!bad.empty()) throw Error("Species must be AW, SB, SW, or PL. Check row(s): " + join_rows(bad));
  std::vector<double> g = check_ti(s, gw, sel_age, warn);

  std::vector<double> out(n);
  for (std::size_t i = 0; i < n; i++) out[i] = SI_TI(s[i], si_t[i], g[i], sel_age[i]);
  return out;
}

} // namespace gypsy

#endif
