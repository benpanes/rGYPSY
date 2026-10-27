#ifndef GYPSY_STAND_H
#define GYPSY_STAND_H

#include "gypsy_eq.h"

// ---------------------------------------------------------------------------
// Stand-level GYPSY components, one stand at a time
//
//   StandPSI: percent stocking index (SAS examp4)
//   StandSDF: stand density factor and initial density (SAS examp5)
//
// Inputs and outputs are indexed by species (AW, SB, SW, PL), because each
// species' value depends on the others (AW -> SB -> SW -> PL).
// ---------------------------------------------------------------------------

namespace gypsy {

// ---- Percent Stocking Index -------------------------------------------------
// PSI = percent stocking at total age 50, capped at 100.
// As in the SAS, missing AW, SB and SW PSIs are 0, PL is NaN.

inline void StandPSI(const double ps[4], const double age_t[4], const double si_b[4], double psi[4]) {
  const double lg51sq = log51*log51;
  double aw_psi = NaN, sb_psi = NaN, sw_psi = NaN, pl_psi = NaN;

  // ---- AW
  if (sv(si_b[AW]) > 0 && sv(age_t[AW]) > 0 && sv(ps[AW]) > 0) {
    const double b1 = -14.9555;
    const double b2 =  0.45162;
    const double b3 =  1.59953;
    double la = std::log(age_t[AW] + 1.0), ls = std::log(si_b[AW]);
    aw_psi = ps[AW]*(1.0 + std::exp(b1 + b2*la*la + b3*ls))/(1.0 + std::exp(b1 + b2*lg51sq + b3*ls));
  }
  if (aw_psi > 100) aw_psi = 100;

  // ---- SB
  if (sv(si_b[SB]) > 0 && sv(age_t[SB]) > 0 && sv(ps[SB]) > 0) {
    const double b1 = -16.7154;
    const double b2 =  0.369843;
    const double b3 =  2.808537;
    double la = std::log(age_t[SB] + 1.0), ls = std::log(si_b[SB]);
    sb_psi = ps[SB]*(1.0 + std::exp(b1 + b2*la*la + b3*ls))/(1.0 + std::exp(b1 + b2*lg51sq + b3*ls));
  }
  if (sb_psi > 100) sb_psi = 100;

  aw_psi = nz(aw_psi);

  // ---- SW (AW PSI enters as z1*(psi_aw/b3)/50; zero when no AW)
  if (sv(si_b[SW]) > 0 && sv(age_t[SW]) > 0 && sv(ps[SW]) > 0) {
    const double b1 = -21.8721;
    const double b2 =  0.395767;
    const double b3 =  4.128582;
    double la = std::log(age_t[SW] + 1.0), ls = std::log(si_b[SW]);
    double z  = aw_psi/b3/50.0;
    sw_psi = ps[SW]*(1.0 + std::exp(b1 + b2*la*la + b3*ls + z))/(1.0 + std::exp(b1 + b2*lg51sq + b3*ls + z));
  }
  if (sw_psi > 100) sw_psi = 100;

  sb_psi = nz(sb_psi);
  sw_psi = nz(sw_psi);

  // ---- PL
  if (sv(si_b[PL]) > 0 && sv(age_t[PL]) > 0 && sv(ps[PL]) > 0) {
    const double b1 = -12.2731;
    const double b2 =  0.23981;
    const double b3 =  2.17610;
    const double b4 =  0.02492;
    const double b5 =  0.08990;
    const double b6 =  0.17058;
    double la = std::log(age_t[PL] + 1.0), ls = std::log(si_b[PL]);
    double z  = b4*aw_psi + b5*std::sqrt(sb_psi) + b6*std::sqrt(sw_psi);
    pl_psi = ps[PL]*(1.0 + std::exp(b1 + b2*la*la + b3*ls + z))/(1.0 + std::exp(b1 + b2*lg51sq + b3*ls + z));
  }
  if (pl_psi > 100) pl_psi = 100;

  psi[AW] = aw_psi; psi[SB] = sb_psi; psi[SW] = sw_psi; psi[PL] = pl_psi;
}

// ---- Stand Density Factor ---------------------------------------------------
// SDF = density at age 50 (AW: BH age 50, others total age 50), solved from
// the observed density. den0 = density at age 0 from that SDF (AW at BH age
// 0, others at total age 0), used by the BA increment models.
// age = AW breast height age, total age for the others.
// As in the SAS, missing AW, SB and SW SDFs are 0, PL is NaN.

inline void StandSDF(const double age[4], const double den[4], const double si_b[4], const double psi[4],
                     int iter, double tol, double sdf[4], double den0[4], bool converged[4]) {
  const double min_den[4] = {20.0, 0.0, 0.0, 0.0};   // SAS: den_aw>20, others den>0

  for (int s = 0; s < 4; s++) { sdf[s] = NaN; den0[s] = NaN; converged[s] = true; }
  for (int s = 0; s < 4; s++) {   // order matters: SW uses AW, PL uses AW, SW, SB
    if (sv(den[s]) > min_den[s] && sv(age[s]) > 0 && sv(si_b[s]) > 0) {
      sdf[s]  = SolveSDF(s, den[s], age[s], si_b[s], psi[s], sdf, iter, tol, converged[s]);
      den0[s] = Density(s, 0.0, sdf[s], si_b[s], psi[s], sdf);
    }
    if (s != PL) sdf[s] = nz(sdf[s]);
  }
}

} // namespace gypsy

#endif
