#include <Rcpp.h>
#include <cmath>
#include <string>
#include "GYPSY_eq.h"
using namespace Rcpp;

// ---------------------------------------------------------------------------
// Total volume, merchantable volume and merchantable density (SAS examp18)
//
// @param spp        "AW", "SB", "SW" or "PL"
// @param ba         basal area (m2/ha)
// @param topht      top height (m)
// @param n          density (stems/ha)
// @param bhage      breast height age
// @param stump_ht   stump height (m, <= 0.3)
// @param stump_dob  minimum stump diameter outside bark (cm)
// @param top_dib    minimum top diameter inside bark (cm)
//
// Follows the SAS conditions: total volume needs top height > 1.3 m;
// merchantable values also need n > 0 (through QMD). Merchantable volume is
// rounded to 1e-5 as in the SAS. Values that can't be computed are NA.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List GYPSY_Vol_cpp(CharacterVector spp,
                         NumericVector   ba,
                         NumericVector   topht,
                         NumericVector   n,
                         NumericVector   bhage,
                         NumericVector   stump_ht,
                         NumericVector   stump_dob,
                         NumericVector   top_dib) {
  using namespace gypsy;
  int len = spp.size();

  NumericVector tvol(len, NA_REAL), mvol(len, NA_REAL), mn(len, NA_REAL), qmd(len, NA_REAL);
  for (int i = 0; i < len; i++) {
    int s = spp_index(as<std::string>(spp[i]));
    if (s < 0 || miss(ba[i])) continue;

    double k = (sv(n[i]) > 0) ? std::sqrt(ba[i]/(n[i]/10000.0)) : NaN;
    if (!miss(k)) qmd[i] = k;

    if (sv(topht[i]) > 1.3) {
      double tv = TotVol(s, ba[i], topht[i]);
      tvol[i] = tv;
      if (sv(k) > 0) {
        MerchResult m = Merch(s, tv, k, topht[i], bhage[i], n[i], stump_ht[i], stump_dob[i], top_dib[i]);
        if (!miss(m.mv)) mvol[i] = sas_round(m.mv, 0.00001);
        if (!miss(m.mn)) mn[i]   = m.mn;
      }
    }
  }
  return List::create(
    _["tvol"] = tvol,
    _["mvol"] = mvol,
    _["mn"]   = mn,
    _["qmd"]  = qmd
  );
}
