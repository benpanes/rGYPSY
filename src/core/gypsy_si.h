#ifndef GYPSY_SI_H
#define GYPSY_SI_H

#include "gypsy_eq.h"

// ---------------------------------------------------------------------------
// Site index and age solvers (fixed-point iteration on the top height
// equations)
//
//   SolveSI:  SI from top height and age (any age type)
//   SolveAge: age from SI and top height; with ht = NaN, also converts
//             between SI types and gives years to breast / stump height
//
// Each species step takes the current si0 guess and returns
//   si0 = (si0_old + si1) / 2  - averaged, becomes the next iteration input
//   si1 = raw si1 before averaging
// Convergence check: abs(si0_averaged - si1_raw) < tol
// ---------------------------------------------------------------------------

namespace gypsy {

// ---- SI from height and age -------------------------------------------------

struct SIStep { double si0, si1, si_bh, y2bh, y2sh; };

// at: age type, 'T' total, 'B' breast height, 'S' stump height
inline SIStep AW_SI(double age, char at, double ht, double si0, double bh, double sh) {
  const double b1 =  9.908888;
  const double b2 = -3.92451;
  const double b3 = -0.32778;
  const double b4 =  0.134376;

  double lsi0 = std::log(si0);
  double k1   = std::exp(b1 + b2*std::sqrt(log51) + b3*lsi0*lsi0 + b4*sq50);
  double k2   = std::pow(si0, b3*lsi0);
  double k3b  = (si0*(1.0+k1)/bh - 1.0) / (std::exp(b1)*std::exp(b4*sq50)*k2);
  double k3s  = (si0*(1.0+k1)/sh - 1.0) / (std::exp(b1)*std::exp(b4*sq50)*k2);
  double Y2BH = std::exp(std::pow(std::log(k3b)/b2, 2)) - 1.0;
  double Y2SH = std::exp(std::pow(std::log(k3s)/b2, 2)) - 1.0;

  double tage;
  if      (at == 'T') tage = age;
  else if (at == 'B') tage = age+Y2BH;
  else                tage = age+Y2SH;

  double x10  = 1.0 + std::exp(b1 + b2*std::sqrt(std::log(tage+1.0)) + b3*lsi0*lsi0 + b4*sq50);
  double x20  = 1.0 + std::exp(b1 + b2*std::sqrt(log51)              + b3*lsi0*lsi0 + b4*sq50);
  double si1  = ht * x10 / x20;
  double si0n = (si0 + si1) / 2.0;

  double lsi1  = std::log(si1);
  double si_bh = si1 * (1.0 + std::exp(b1 + b2*std::sqrt(log51)                   + b3*lsi1*lsi1 + b4*sq50)) /
                       (1.0 + std::exp(b1 + b2*std::sqrt(std::log(50.0+Y2BH+1.0)) + b3*lsi1*lsi1 + b4*sq50));
  return {si0n, si1, si_bh, Y2BH, Y2SH};
}

inline SIStep SB_SI(double age, char at, double ht, double si0, double bh, double sh) {
  const double b1 =  14.56236;
  const double b2 =  -6.04705;
  const double b3 =  -1.53715;
  const double b4 =   0.240174;

  double lsi0 = std::log(si0);
  double k1   = std::exp(b1 + b2*std::sqrt(log51) + b3*lsi0 + b4*sq50);
  double k2   = std::pow(si0, b3);
  double k3b  = (si0*(1.0+k1)/bh - 1.0) / (std::exp(b1)*std::exp(b4*sq50)*k2);
  double k3s  = (si0*(1.0+k1)/sh - 1.0) / (std::exp(b1)*std::exp(b4*sq50)*k2);
  double Y2BH = std::exp(std::pow(std::log(k3b)/b2, 2)) - 1.0;
  double Y2SH = std::exp(std::pow(std::log(k3s)/b2, 2)) - 1.0;

  double tage;
  if      (at == 'T') tage = age;
  else if (at == 'B') tage = age+Y2BH;
  else                tage = age+Y2SH;

  double x10  = 1.0 + std::exp(b1 + b2*std::sqrt(std::log(tage+1.0)) + b3*lsi0 + b4*sq50);
  double x20  = 1.0 + std::exp(b1 + b2*std::sqrt(log51)              + b3*lsi0 + b4*sq50);
  double si1  = ht * x10 / x20;
  double si0n = (si0 + si1) / 2.0;

  double lsi1  = std::log(si1);
  double si_bh = si1 * (1.0 + std::exp(b1 + b2*std::sqrt(log51)                   + b3*lsi1 + b4*sq50)) /
                       (1.0 + std::exp(b1 + b2*std::sqrt(std::log(50.0+Y2BH+1.0)) + b3*lsi1 + b4*sq50));
  return {si0n, si1, si_bh, Y2BH, Y2SH};
}

inline SIStep SW_SI(double age, char at, double ht, double si0, double bh, double sh) {
  const double b1 =  12.14943;
  const double b2 =  -3.77051;
  const double b3 =  -0.28534;
  const double b4 =   0.165483;

  double lsi0 = std::log(si0);
  double k1   = std::exp(b1 + b2*std::sqrt(log2501) + b3*lsi0*lsi0 + b4*sq50);
  double k2   = std::pow(si0, b3*lsi0);
  double k3b  = (si0*(1.0+k1)/bh - 1.0) / (std::exp(b1)*std::exp(b4*sq50)*k2);
  double k3s  = (si0*(1.0+k1)/sh - 1.0) / (std::exp(b1)*std::exp(b4*sq50)*k2);
  double Y2BH = std::sqrt(std::exp(std::pow(std::log(k3b)/b2, 2)) - 1.0);
  double Y2SH = std::sqrt(std::exp(std::pow(std::log(k3s)/b2, 2)) - 1.0);

  double tage;
  if      (at == 'T') tage = age;
  else if (at == 'B') tage = age+Y2BH;
  else                tage = age+Y2SH;

  double x10  = 1.0 + std::exp(b1 + b2*std::sqrt(std::log(tage*tage+1.0)) + b3*lsi0*lsi0 + b4*sq50);
  double x20  = 1.0 + std::exp(b1 + b2*std::sqrt(log2501)                  + b3*lsi0*lsi0 + b4*sq50);
  double si1  = ht * x10 / x20;
  double si0n = (si0 + si1) / 2.0;

  double lsi1  = std::log(si1);
  double y50   = 50.0 + Y2BH;
  double si_bh = si1 * (1.0 + std::exp(b1 + b2*std::sqrt(log2501)                  + b3*lsi1*lsi1 + b4*sq50)) /
                       (1.0 + std::exp(b1 + b2*std::sqrt(std::log(y50*y50 + 1.0)) + b3*lsi1*lsi1 + b4*sq50));
  return {si0n, si1, si_bh, Y2BH, Y2SH};
}

inline SIStep PL_SI(double age, char at, double ht, double si0, double bh, double sh) {
  const double b1 =  12.84571;
  const double b2 =  -5.73936;
  const double b3 =  -0.91312;
  const double b4 =   0.150668;

  double lsi0 = std::log(si0);
  double k1   = std::exp(b1 + b2*std::sqrt(log51) + b3*lsi0 + b4*sq50);
  double k2   = std::pow(si0, b3);
  double k3b  = (si0*(1.0+k1)/bh - 1.0) / (std::exp(b1)*std::exp(b4*sq50)*k2);
  double k3s  = (si0*(1.0+k1)/sh - 1.0) / (std::exp(b1)*std::exp(b4*sq50)*k2);
  double Y2BH = std::exp(std::pow(std::log(k3b)/b2, 2)) - 1.0;
  double Y2SH = std::exp(std::pow(std::log(k3s)/b2, 2)) - 1.0;

  double tage;
  if      (at == 'T') tage = age;
  else if (at == 'B') tage = age+Y2BH;
  else                tage = age+Y2SH;

  double x10  = 1.0 + std::exp(b1 + b2*std::sqrt(std::log(tage+1.0)) + b3*lsi0 + b4*sq50);
  double x20  = 1.0 + std::exp(b1 + b2*std::sqrt(log51)              + b3*lsi0 + b4*sq50);
  double si1  = ht * x10 / x20;
  double si0n = (si0 + si1) / 2.0;

  double lsi1  = std::log(si1);
  double si_bh = si1 * (1.0 + std::exp(b1 + b2*std::sqrt(log51)                   + b3*lsi1 + b4*sq50)) /
                       (1.0 + std::exp(b1 + b2*std::sqrt(std::log(50.0+Y2BH+1.0)) + b3*lsi1 + b4*sq50));
  return {si0n, si1, si_bh, Y2BH, Y2SH};
}

struct SIResult { double si_t, si_bh, y2bh, y2sh, diff; bool converged; };

// SI from top height and age, starting from a guess of 10
inline SIResult SolveSI(int s, double ht, double age, char at,
                        int iter, double tol, double bh, double sh) {
  double si0 = 10.0, si1 = 0.0, si_bh = 0.0, y2bh = 0.0, y2sh = 0.0;
  double diff = 10.0;
  bool converged = false;

  for (int it = 0; it < iter; it++) {
    SIStep res;
    if      (s == AW) res = AW_SI(age, at, ht, si0, bh, sh);
    else if (s == SB) res = SB_SI(age, at, ht, si0, bh, sh);
    else if (s == SW) res = SW_SI(age, at, ht, si0, bh, sh);
    else              res = PL_SI(age, at, ht, si0, bh, sh);

    si0   = res.si0;
    si1   = res.si1;
    si_bh = res.si_bh;
    y2bh  = res.y2bh;
    y2sh  = res.y2sh;

    diff = std::abs(si0 - si1);
    if (diff < tol) {
      converged = true;
      break;
    }
  }
  return {si1, si_bh, y2bh, y2sh, diff, converged};
}

// ---- Age from SI and height -------------------------------------------------

struct AgeStep { double si0, si1, sibh, tage, y2bh, y2sh; };

// st: SI type, 'B' (sibh given, SI_t solved by iteration) or 'T' (si0 is
// SI_t, one pass)
inline AgeStep AW_Age(double sibh, char st, double ht, double si0, double bh, double sh) {
  const double b1 =  9.908888;
  const double b2 = -3.92451;
  const double b3 = -0.32778;
  const double b4 =  0.134376;

  double lsi0 = std::log(si0);
  double k1   = std::exp(b1 + b2*std::sqrt(log51) + b3*lsi0*lsi0 + b4*sq50);
  double k2   = std::pow(si0, b3*lsi0);
  double k3b  = (si0*(1.0+k1)/bh - 1.0) / (std::exp(b1)*std::exp(b4*sq50)*k2);
  double k3s  = (si0*(1.0+k1)/sh - 1.0) / (std::exp(b1)*std::exp(b4*sq50)*k2);
  double Y2BH = std::exp(std::pow(std::log(k3b)/b2, 2)) - 1.0;
  double Y2SH = std::exp(std::pow(std::log(k3s)/b2, 2)) - 1.0;
  double x10  = 1.0 + std::exp(b1 + b2*std::sqrt(std::log(50+Y2BH+1.0)) + b3*lsi0*lsi0 + b4*sq50);
  double x20  = 1.0 + std::exp(b1 + b2*std::sqrt(log51)                 + b3*lsi0*lsi0 + b4*sq50);

  double si1 = sibh*x10/x20;

  double si0n, sibhn;
  if (st == 'B') { si0n = (si0 + si1) / 2.0; sibhn = sibh; }
  else           { si0n = si0;               sibhn = si0*x20/x10; }

  double k3x       = (si0n*(1+k1)/ht-1)/(std::exp(b1)*std::exp(b4*sq50)*k2);
  double logk3x_b2 = std::log(k3x)/b2;
  double tage      = std::exp(logk3x_b2*logk3x_b2)-1;

  return {si0n, si1, sibhn, tage, Y2BH, Y2SH};
}

inline AgeStep SB_Age(double sibh, char st, double ht, double si0, double bh, double sh) {
  const double b1 =  14.56236;
  const double b2 =  -6.04705;
  const double b3 =  -1.53715;
  const double b4 =   0.240174;

  double lsi0 = std::log(si0);
  double k1   = std::exp(b1 + b2*std::sqrt(log51) + b3*lsi0 + b4*sq50);
  double k2   = std::pow(si0, b3);
  double k3b  = (si0*(1.0+k1)/bh - 1.0) / (std::exp(b1)*std::exp(b4*sq50)*k2);
  double k3s  = (si0*(1.0+k1)/sh - 1.0) / (std::exp(b1)*std::exp(b4*sq50)*k2);
  double Y2BH = std::exp(std::pow(std::log(k3b)/b2, 2)) - 1.0;
  double Y2SH = std::exp(std::pow(std::log(k3s)/b2, 2)) - 1.0;
  double x10  = 1.0 + std::exp(b1 + b2*std::sqrt(std::log(50+Y2BH+1.0)) + b3*lsi0 + b4*sq50);
  double x20  = 1.0 + std::exp(b1 + b2*std::sqrt(log51)                 + b3*lsi0 + b4*sq50);

  double si1 = sibh*x10/x20;

  double si0n, sibhn;
  if (st == 'B') { si0n = (si0 + si1) / 2.0; sibhn = sibh; }
  else           { si0n = si0;               sibhn = si0*x20/x10; }

  double k3x       = (si0n*(1+k1)/ht-1)/(std::exp(b1)*std::exp(b4*sq50)*k2);
  double logk3x_b2 = std::log(k3x)/b2;
  double tage      = std::exp(logk3x_b2*logk3x_b2)-1;

  return {si0n, si1, sibhn, tage, Y2BH, Y2SH};
}

inline AgeStep SW_Age(double sibh, char st, double ht, double si0, double bh, double sh) {
  const double b1 =  12.14943;
  const double b2 =  -3.77051;
  const double b3 =  -0.28534;
  const double b4 =   0.165483;

  double lsi0 = std::log(si0);
  double k1   = std::exp(b1 + b2*std::sqrt(log2501) + b3*lsi0*lsi0 + b4*sq50);
  double k2   = std::pow(si0, b3*lsi0);
  double k3b  = (si0*(1.0+k1)/bh - 1.0) / (std::exp(b1)*std::exp(b4*sq50)*k2);
  double k3s  = (si0*(1.0+k1)/sh - 1.0) / (std::exp(b1)*std::exp(b4*sq50)*k2);
  double Y2BH = std::sqrt(std::exp(std::pow(std::log(k3b)/b2, 2)) - 1.0);
  double Y2SH = std::sqrt(std::exp(std::pow(std::log(k3s)/b2, 2)) - 1.0);
  double x10  = 1.0 + std::exp(b1 + b2*std::sqrt(std::log((50+Y2BH)*(50+Y2BH)+1.0)) + b3*lsi0*lsi0 + b4*sq50);
  double x20  = 1.0 + std::exp(b1 + b2*std::sqrt(log2501)                           + b3*lsi0*lsi0 + b4*sq50);

  double si1 = sibh*x10/x20;

  double si0n, sibhn;
  if (st == 'B') { si0n = (si0 + si1) / 2.0; sibhn = sibh; }
  else           { si0n = si0;               sibhn = si0*x20/x10; }

  double k3x       = (si0n*(1+k1)/ht-1)/(std::exp(b1)*std::exp(b4*sq50)*k2);
  double logk3x_b2 = std::log(k3x)/b2;
  double tage      = std::sqrt(std::exp(logk3x_b2*logk3x_b2)-1);

  return {si0n, si1, sibhn, tage, Y2BH, Y2SH};
}

inline AgeStep PL_Age(double sibh, char st, double ht, double si0, double bh, double sh) {
  const double b1 =  12.84571;
  const double b2 =  -5.73936;
  const double b3 =  -0.91312;
  const double b4 =   0.150668;

  double lsi0 = std::log(si0);
  double k1   = std::exp(b1 + b2*std::sqrt(log51) + b3*lsi0 + b4*sq50);
  double k2   = std::pow(si0, b3);
  double k3b  = (si0*(1.0+k1)/bh - 1.0) / (std::exp(b1)*std::exp(b4*sq50)*k2);
  double k3s  = (si0*(1.0+k1)/sh - 1.0) / (std::exp(b1)*std::exp(b4*sq50)*k2);
  double Y2BH = std::exp(std::pow(std::log(k3b)/b2, 2)) - 1.0;
  double Y2SH = std::exp(std::pow(std::log(k3s)/b2, 2)) - 1.0;
  double x10  = 1.0 + std::exp(b1 + b2*std::sqrt(std::log(50+Y2BH+1.0)) + b3*lsi0 + b4*sq50);
  double x20  = 1.0 + std::exp(b1 + b2*std::sqrt(log51)                 + b3*lsi0 + b4*sq50);

  double si1 = sibh*x10/x20;

  double si0n, sibhn;
  if (st == 'B') { si0n = (si0 + si1) / 2.0; sibhn = sibh; }
  else           { si0n = si0;               sibhn = si0*x20/x10; }

  double k3x       = (si0n*(1+k1)/ht-1)/(std::exp(b1)*std::exp(b4*sq50)*k2);
  double logk3x_b2 = std::log(k3x)/b2;
  double tage      = std::exp(logk3x_b2*logk3x_b2)-1;

  return {si0n, si1, sibhn, tage, Y2BH, Y2SH};
}

// tage = total age; stump and breast height age are tage - y2sh, tage - y2bh
struct AgeResult { double tage, y2bh, y2sh, si_t, si_bh, diff; bool converged; };

// Age from SI and top height (ht may be NaN when only the SI conversion or
// years to breast / stump height are needed).
//   st = 'B': si is SI_bh; SI_t is solved by iteration from si
//   st = 'T': si is SI_t; one pass
inline AgeResult SolveAge(int s, double si, double ht, char st,
                          int iter, double tol, double bh, double sh) {
  auto step = [&](double sibh, double si0) {
    if      (s == AW) return AW_Age(sibh, st, ht, si0, bh, sh);
    else if (s == SB) return SB_Age(sibh, st, ht, si0, bh, sh);
    else if (s == SW) return SW_Age(sibh, st, ht, si0, bh, sh);
    else              return PL_Age(sibh, st, ht, si0, bh, sh);
  };

  if (st == 'B') {
    double sibh = si, si0 = si;   // BH index as the starting guess
    double si1 = 0.0, tage = 0.0, y2bh = 0.0, y2sh = 0.0, diff = 10.0;
    bool converged = false;
    for (int it = 0; it < iter; it++) {
      AgeStep res = step(sibh, si0);
      si0  = res.si0;
      si1  = res.si1;
      sibh = res.sibh;
      tage = res.tage;
      y2bh = res.y2bh;
      y2sh = res.y2sh;

      diff = std::abs(si0 - si1);
      if (diff < tol) {
        converged = true;
        break;
      }
    }
    return {tage, y2bh, y2sh, si1, sibh, diff, converged};
  }

  AgeStep res = step(0.0, si);
  return {res.tage, res.y2bh, res.y2sh, si, res.sibh, 0.0, true};
}

} // namespace gypsy

#endif
