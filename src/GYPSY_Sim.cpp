#include <Rcpp.h>
#include <cmath>
#include <string>
#include "GYPSY_eq.h"
using namespace Rcpp;

// ---------------------------------------------------------------------------
// Stand-level GYPSY components (equations in GYPSY_eq.h)
//
//   GYPSY_PSI_cpp: percent stocking index (SAS examp4)
//   GYPSY_SDF_cpp: stand density factor and initial density (SAS examp5)
//   GYPSY_Ht_cpp:  top height from SI_t and total age
//   GYPSY_N_cpp:   density at age from SDF
//
// PSI and SDF take one row per stand with species in columns, because each
// species' value depends on the others (AW -> SB -> SW -> PL).
// ---------------------------------------------------------------------------

// ---- Percent Stocking Index -----------------------------------------------
// PSI = percent stocking at total age 50, capped at 100.
// As in the SAS, missing AW, SB and SW PSIs are returned as 0, PL as NA.

// [[Rcpp::export]]
Rcpp::List GYPSY_PSI_cpp(CharacterVector stand,

                         NumericVector   aw_ps,
                         NumericVector   aw_age_t,
                         NumericVector   aw_si_b,

                         NumericVector   sb_ps,
                         NumericVector   sb_age_t,
                         NumericVector   sb_si_b,

                         NumericVector   sw_ps,
                         NumericVector   sw_age_t,
                         NumericVector   sw_si_b,

                         NumericVector   pl_ps,
                         NumericVector   pl_age_t,
                         NumericVector   pl_si_b) {
  using namespace gypsy;
  int n = stand.size();
  const double lg51sq = log51*log51;

  NumericVector aw_res(n), sb_res(n), sw_res(n), pl_res(n);
  for (int i = 0; i < n; i++) {
    double aw_psi = NaN, sb_psi = NaN, sw_psi = NaN, pl_psi = NaN;

    // ---- AW
    if (sv(aw_si_b[i]) > 0 && sv(aw_age_t[i]) > 0 && sv(aw_ps[i]) > 0) {
      const double b1 = -14.9555;
      const double b2 =  0.45162;
      const double b3 =  1.59953;
      double la = std::log(aw_age_t[i] + 1.0), ls = std::log(aw_si_b[i]);
      aw_psi = aw_ps[i]*(1.0 + std::exp(b1 + b2*la*la + b3*ls))/(1.0 + std::exp(b1 + b2*lg51sq + b3*ls));
    }
    if (aw_psi > 100) aw_psi = 100;

    // ---- SB
    if (sv(sb_si_b[i]) > 0 && sv(sb_age_t[i]) > 0 && sv(sb_ps[i]) > 0) {
      const double b1 = -16.7154;
      const double b2 =  0.369843;
      const double b3 =  2.808537;
      double la = std::log(sb_age_t[i] + 1.0), ls = std::log(sb_si_b[i]);
      sb_psi = sb_ps[i]*(1.0 + std::exp(b1 + b2*la*la + b3*ls))/(1.0 + std::exp(b1 + b2*lg51sq + b3*ls));
    }
    if (sb_psi > 100) sb_psi = 100;

    aw_psi = nz(aw_psi);

    // ---- SW (AW PSI enters as z1*(psi_aw/b3)/50; zero when no AW)
    if (sv(sw_si_b[i]) > 0 && sv(sw_age_t[i]) > 0 && sv(sw_ps[i]) > 0) {
      const double b1 = -21.8721;
      const double b2 =  0.395767;
      const double b3 =  4.128582;
      double la = std::log(sw_age_t[i] + 1.0), ls = std::log(sw_si_b[i]);
      double z  = aw_psi/b3/50.0;
      sw_psi = sw_ps[i]*(1.0 + std::exp(b1 + b2*la*la + b3*ls + z))/(1.0 + std::exp(b1 + b2*lg51sq + b3*ls + z));
    }
    if (sw_psi > 100) sw_psi = 100;

    sb_psi = nz(sb_psi);
    sw_psi = nz(sw_psi);

    // ---- PL
    if (sv(pl_si_b[i]) > 0 && sv(pl_age_t[i]) > 0 && sv(pl_ps[i]) > 0) {
      const double b1 = -12.2731;
      const double b2 =  0.23981;
      const double b3 =  2.17610;
      const double b4 =  0.02492;
      const double b5 =  0.08990;
      const double b6 =  0.17058;
      double la = std::log(pl_age_t[i] + 1.0), ls = std::log(pl_si_b[i]);
      double z  = b4*aw_psi + b5*std::sqrt(sb_psi) + b6*std::sqrt(sw_psi);
      pl_psi = pl_ps[i]*(1.0 + std::exp(b1 + b2*la*la + b3*ls + z))/(1.0 + std::exp(b1 + b2*lg51sq + b3*ls + z));
    }
    if (pl_psi > 100) pl_psi = 100;

    aw_res[i] = aw_psi;
    sb_res[i] = sb_psi;
    sw_res[i] = sw_psi;
    pl_res[i] = miss(pl_psi) ? NA_REAL : pl_psi;
  }
  return List::create(
    _["aw_psi"] = aw_res,
    _["sb_psi"] = sb_res,
    _["sw_psi"] = sw_res,
    _["pl_psi"] = pl_res
  );
}

// ---- Stand Density Factor -------------------------------------------------
// SDF = density at age 50 (AW: BH age 50, others total age 50), solved from
// the observed density. den0 = density at age 0 from that SDF (AW at BH age
// 0, others at total age 0), used by the BA increment models.
// As in the SAS, missing AW, SB and SW SDFs are returned as 0, PL as NA.

// [[Rcpp::export]]
Rcpp::List GYPSY_SDF_cpp(CharacterVector stand,

                         NumericVector   aw_age_b,
                         NumericVector   aw_n,
                         NumericVector   aw_si_b,
                         NumericVector   aw_psi,

                         NumericVector   sb_age_t,
                         NumericVector   sb_n,
                         NumericVector   sb_si_b,
                         NumericVector   sb_psi,

                         NumericVector   sw_age_t,
                         NumericVector   sw_n,
                         NumericVector   sw_si_b,
                         NumericVector   sw_psi,

                         NumericVector   pl_age_t,
                         NumericVector   pl_n,
                         NumericVector   pl_si_b,
                         NumericVector   pl_psi,

                         int             iter = 1000,
                         double          tol  = 1e-8) {
  using namespace gypsy;
  int n = stand.size();

  const NumericVector* age[4] = {&aw_age_b, &sb_age_t, &sw_age_t, &pl_age_t};
  const NumericVector* den[4] = {&aw_n,     &sb_n,     &sw_n,     &pl_n};
  const NumericVector* sib[4] = {&aw_si_b,  &sb_si_b,  &sw_si_b,  &pl_si_b};
  const NumericVector* psi[4] = {&aw_psi,   &sb_psi,   &sw_psi,   &pl_psi};
  const double min_den[4] = {20.0, 0.0, 0.0, 0.0};   // SAS: den_aw>20, others den>0

  NumericVector sdf_res[4], den0_res[4];
  for (int s = 0; s < 4; s++) {
    sdf_res[s]  = NumericVector(n);
    den0_res[s] = NumericVector(n);
  }

  for (int i = 0; i < n; i++) {
    double sdf[4] = {NaN, NaN, NaN, NaN};
    for (int s = 0; s < 4; s++) {   // order matters: SW uses AW, PL uses AW, SW, SB
      double a = (*age[s])[i], d = (*den[s])[i], si = (*sib[s])[i], p = (*psi[s])[i];
      double den0 = NaN;
      if (sv(d) > min_den[s] && sv(a) > 0 && sv(si) > 0) {
        bool converged;
        sdf[s] = SolveSDF(s, d, a, si, p, sdf, iter, tol, converged);
        if (!converged)
          Rcpp::warning("Row " + std::to_string(i + 1) + ": SDF did not converge.");
        den0 = Density(s, 0.0, sdf[s], si, p, sdf);
      }
      if (s != PL) sdf[s] = nz(sdf[s]);
      sdf_res[s][i]  = miss(sdf[s]) ? NA_REAL : sdf[s];
      den0_res[s][i] = miss(den0)   ? NA_REAL : den0;
    }
  }
  return List::create(
    _["aw_sdf"]  = sdf_res[AW],  _["sb_sdf"]  = sdf_res[SB],
    _["sw_sdf"]  = sdf_res[SW],  _["pl_sdf"]  = sdf_res[PL],
    _["aw_den0"] = den0_res[AW], _["sb_den0"] = den0_res[SB],
    _["sw_den0"] = den0_res[SW], _["pl_den0"] = den0_res[PL]
  );
}

// ---- Exported top height function -----------------------------------------
// [[Rcpp::export]]
NumericVector GYPSY_Ht_cpp(CharacterVector spp,
                           NumericVector   si_t,
                           NumericVector   age_t) {
  int n = spp.size();

  NumericVector res(n);
  for (int i = 0; i < n; i++) {
    int s = gypsy::spp_index(as<std::string>(spp[i]));
    double ht = (s < 0) ? NA_REAL : gypsy::Ht(s, si_t[i], age_t[i]);
    res[i] = std::isnan(ht) ? NA_REAL : ht;
  }
  return res;
}

// ---- Exported density function --------------------------------------------
// Density at age (AW: BH age, others: total age) from the species' SDF.
// sdf_aw, sdf_sw and sdf_sb are the other species' SDFs (SW uses AW;
// PL uses AW, SW and SB); 0 or NA when absent.

// [[Rcpp::export]]
NumericVector GYPSY_N_cpp(CharacterVector spp,
                          NumericVector   age,
                          NumericVector   sdf,
                          NumericVector   si_b,
                          NumericVector   psi,
                          NumericVector   sdf_aw,
                          NumericVector   sdf_sw,
                          NumericVector   sdf_sb) {
  using namespace gypsy;
  int n = spp.size();

  NumericVector res(n);
  for (int i = 0; i < n; i++) {
    int s = spp_index(as<std::string>(spp[i]));
    double d = NaN;
    if (s >= 0 && sv(sdf[i]) > 0) {
      double comp[4] = {sdf_aw[i], sdf_sb[i], sdf_sw[i], NaN};
      d = Density(s, age[i], sdf[i], si_b[i], psi[i], comp);
    }
    res[i] = miss(d) ? NA_REAL : d;
  }
  return res;
}
