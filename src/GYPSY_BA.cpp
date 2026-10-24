#include <Rcpp.h>
#include <cmath>
#include <string>
#include "GYPSY_eq.h"
using namespace Rcpp;

// ---------------------------------------------------------------------------
// Annual basal area increment (m2/ha/yr), equations in GYPSY_eq.h
//
// @param spp     "AW", "SB", "SW" or "PL"
// @param bhage   breast height age
// @param sc      breast height age based species composition (0-1)
// @param den0    initial density (AW at BH age 0, others at total age 0)
// @param si_b    site index at BH age 50
// @param psi     percent stocking index (0 when unknown)
// @param ba      current basal area of the species
// @param sdf_aw, sdf_sb, sdf_sw, sdf_pl
//                competitor SDFs used by SW and PL; pass 0 for species
//                that have not reached breast height
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector GYPSY_BAinc_cpp(CharacterVector spp,
                              NumericVector   bhage,
                              NumericVector   sc,
                              NumericVector   den0,
                              NumericVector   si_b,
                              NumericVector   psi,
                              NumericVector   ba,
                              NumericVector   sdf_aw,
                              NumericVector   sdf_sb,
                              NumericVector   sdf_sw,
                              NumericVector   sdf_pl) {
  using namespace gypsy;
  int n = spp.size();

  NumericVector res(n);
  for (int i = 0; i < n; i++) {
    int s = spp_index(as<std::string>(spp[i]));
    double inc = NaN;
    if (s >= 0) {
      double comp[4] = {nz(sdf_aw[i]), nz(sdf_sb[i]), nz(sdf_sw[i]), nz(sdf_pl[i])};
      inc = BAinc(s, bhage[i], sc[i], den0[i], si_b[i], psi[i], ba[i], comp);
    }
    res[i] = miss(inc) ? NA_REAL : inc;
  }
  return res;
}
