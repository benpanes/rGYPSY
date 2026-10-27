#ifndef GYPSY_EQ_H
#define GYPSY_EQ_H

#include <cmath>
#include <deque>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// GYPSY component equations (MAY 21, 2009 functions, GYPSY_run.sas)
//
// Part of the GYPSY core (src/core): plain C++ with no R or Rcpp
// dependencies, so it can also be compiled outside R (e.g. to WebAssembly).
//
// Missing values are NaN. SAS treats a missing value as smaller than any
// number, so conditions copied from the SAS go through sv(), which maps
// NaN to -Inf (e.g. "if bhage1_aw<0" is sv(bhage1) < 0).
//
// Species index: 0 = AW, 1 = SB, 2 = SW, 3 = PL (SAS column order).
// ---------------------------------------------------------------------------

namespace gypsy {

enum Spp { AW = 0, SB = 1, SW = 2, PL = 3 };

const double NaN = std::numeric_limits<double>::quiet_NaN();

constexpr double sq50    = 7.0710678118654752;  // sqrt(50)
constexpr double log51   = 3.9318256327243257;  // log(50 + 1)
constexpr double log2501 = 7.824445930877619;   // log(50^2 + 1), SW only

inline bool   miss(double x) { return std::isnan(x); }
inline double sv(double x)   { return std::isnan(x) ? -std::numeric_limits<double>::infinity() : x; }
inline double nz(double x)   { return std::isnan(x) ? 0.0 : x; }   // SAS "if x=. then x=0"

// SAS round(x, unit)
inline double sas_round(double x, double unit) {
  return std::isnan(x) ? NaN : std::round(x/unit)*unit;
}

// SAS sum(): ignores missing, missing if all missing
inline double sas_sum4(const double* x) {
  double s = 0.0; bool any = false;
  for (int i = 0; i < 4; i++) if (!std::isnan(x[i])) { s += x[i]; any = true; }
  return any ? s : NaN;
}

inline int spp_index(const std::string& s) {
  if (s == "AW") return AW;
  if (s == "SB") return SB;
  if (s == "SW") return SW;
  if (s == "PL") return PL;
  return -1;
}

// ---- Core types -------------------------------------------------------------

// Input errors are thrown; warnings are collected and returned to the caller
struct Error : std::runtime_error { using std::runtime_error::runtime_error; };
using Messages = std::vector<std::string>;

// Ordered named columns (missing = NaN); a deque keeps column references
// valid while columns are added
struct Table {
  std::vector<std::string>         names;
  std::deque<std::vector<double>>  cols;
  std::vector<double>& add(const std::string& name, std::size_t n, double fill = NaN) {
    names.push_back(name);
    cols.emplace_back(n, fill);
    return cols.back();
  }
};

// 1-based row numbers as "1, 2, 3"
inline std::string join_rows(const std::vector<int>& rows) {
  std::string out;
  for (std::size_t i = 0; i < rows.size(); i++) out += (i ? ", " : "") + std::to_string(rows[i] + 1);
  return out;
}

// ---- Top height -----------------------------------------------------------
// Height from total age SI and total age

inline double AW_Ht(double si_t, double age_t) {
  const double b1 =  9.908888;
  const double b2 = -3.92451;
  const double b3 = -0.32778;
  const double b4 =  0.134376;

  double lsit = std::log(si_t);

  // Numerator and denominator set as separate variables for readability
  double x10 = 1.0 + std::exp(b1 + b2*std::sqrt(log51)               + b3*lsit*lsit + b4*sq50);
  double x20 = 1.0 + std::exp(b1 + b2*std::sqrt(std::log(age_t+1.0)) + b3*lsit*lsit + b4*sq50);

  return si_t*x10/x20;
}

inline double SB_Ht(double si_t, double age_t) {
  const double b1 =  14.56236;
  const double b2 =  -6.04705;
  const double b3 =  -1.53715;
  const double b4 =   0.240174;

  double lsit = std::log(si_t);

  double x10 = 1.0 + std::exp(b1 + b2*std::sqrt(log51)               + b3*lsit + b4*sq50);
  double x20 = 1.0 + std::exp(b1 + b2*std::sqrt(std::log(age_t+1.0)) + b3*lsit + b4*sq50);

  return si_t*x10/x20;
}

inline double SW_Ht(double si_t, double age_t) {
  const double b1 =  12.14943;
  const double b2 =  -3.77051;
  const double b3 =  -0.28534;
  const double b4 =   0.165483;

  double lsit = std::log(si_t);

  double x10 = 1.0 + std::exp(b1 + b2*std::sqrt(log2501)                     + b3*lsit*lsit + b4*sq50);
  double x20 = 1.0 + std::exp(b1 + b2*std::sqrt(std::log(age_t*age_t+1.0))   + b3*lsit*lsit + b4*sq50);

  return si_t*x10/x20;
}

inline double PL_Ht(double si_t, double age_t) {
  const double b1 =  12.84571;
  const double b2 =  -5.73936;
  const double b3 =  -0.91312;
  const double b4 =   0.150668;

  double lsit = std::log(si_t);

  double x10 = 1.0 + std::exp(b1 + b2*std::sqrt(log51)               + b3*lsit + b4*sq50);
  double x20 = 1.0 + std::exp(b1 + b2*std::sqrt(std::log(age_t+1.0)) + b3*lsit + b4*sq50);

  return si_t*x10/x20;
}

inline double Ht(int s, double si_t, double age_t) {
  switch (s) {
    case AW: return AW_Ht(si_t, age_t);
    case SB: return SB_Ht(si_t, age_t);
    case SW: return SW_Ht(si_t, age_t);
    case PL: return PL_Ht(si_t, age_t);
  }
  return NaN;
}

// ---- Density (spatial density functions) ----------------------------------
// N(age) = SDF * kk(50) / kk(age), kk(a) = 1 + exp(b1 + b2*x + b3*log(a + 1))
//   x = SI_bh for AW, log(SI_bh) for SB, SW, PL
//   age = breast height age for AW, total age for SB, SW, PL
// b1-b3 depend on SDF, on PSI (PS > 0 branch), and for SW and PL on the
// other species' SDFs. The same coefficients are used to solve for SDF.
//
// The PS > 0 / PS <= 0 SAS branches are selected with psi > 0: psi is only
// positive when PS > 0 and the species has SI and age, which every caller
// requires anyway.

struct DenCoef { double b1, b2, b3; };

inline DenCoef AW_DenCoef(double sdf, double psi) {
  double sqrtsdf = std::sqrt(sdf);
  if (psi > 0) {
    const double c0 = 0.71123;
    const double c1 = 7.015511;
    return { -(1.0/std::sqrt(sdf/1000.0) + std::sqrt(1.0 + 50.0/(sqrtsdf*log51))
               + c0*std::log(1.0 + psi/sqrtsdf)/c1)*log51,
             c0*std::pow(sqrtsdf, 1.0/sdf)/4.0,
             (1.0 + c0)*std::pow(sdf, (c1 + std::log(sdf) + std::log(1.0 + psi)/sdf)/sdf) };
  }
  const double c0 = 0.717966;
  const double c1 = 6.67468;
  return { -(1.0/std::sqrt(sdf/1000.0) + std::sqrt(1.0 + std::sqrt(50.0/(sqrtsdf*log51))))*log51,
           c0*std::pow(sqrtsdf, 1.0/sdf)/4.0,
           (1.0 + c0)*std::pow(sdf, (c1 + std::log(sdf))/sdf) };
}

inline DenCoef SB_DenCoef(double sdf, double psi) {
  if (psi > 0) {
    const double c1 = -25.7666;
    const double c2 =  0.192351;
    const double c3 =  2.718748;
    return { c1/std::pow(std::sqrt(sdf/1000.0) + log51, c2) - std::log(1.0 + psi)/(1.0 + c3),
             c3,
             c3*std::pow(sdf, 1.0/sdf) };
  }
  const double c1 = -26.3836;
  const double c2 =  0.166483;
  const double c3 =  2.738569;
  return { c1/std::pow(std::sqrt(sdf/1000.0) + log51, c2),
           c3,
           c3*std::pow(sdf, 1.0/sdf) };
}

inline DenCoef SW_DenCoef(double sdf, double psi, double sdf_aw) {
  double z1 = (sdf_aw > 0) ? std::sqrt(1.0 + sdf_aw/1000.0) : 0.0;
  if (psi > 0) {
    const double c1 = -254.391;
    const double c2 =  1.146709;
    const double c3 =  2.202149;
    return { c1/std::pow(std::log(sdf) + log51, c2) + z1 + std::log(1.0 + sdf)/10.0,
             c3,
             c3*std::pow(sdf, 1.0/sdf)*std::pow(std::sqrt(1.0 + psi), -1.0/10.0) };
  }
  const double c1 = -231.617;
  const double c2 =  1.176995;
  const double c3 =  1.733601;
  return { c1/std::pow(std::log(sdf) + log51, c2) + z1,
           c3,
           c3*std::pow(sdf, 1.0/sdf) };
}

inline DenCoef PL_DenCoef(double sdf, double psi, double sdf_aw, double sdf_sw, double sdf_sb) {
  double sqrtsdf = std::sqrt(sdf);
  double z = ((sdf_aw > 0) ? sdf_aw/1000.0/2.0 : 0.0)
           + ((sdf_sw > 0) ? sdf_sw/1000.0/3.0 : 0.0)
           + ((sdf_sb > 0) ? sdf_sb/1000.0/4.0 : 0.0);
  if (psi > 0) {
    const double c1 = -3.09711;
    const double c2 = -512.394;
    const double c3 =  1.14981;
    const double c4 =  1.027108;
    const double c5 = -0.05479;
    const double c6 =  4.12191;
    return { c1 + z + c2/std::pow(sqrtsdf, c3) - std::log(1.0 + psi)/(1.0 + c4),
             c4/std::pow(sqrtsdf, c5),
             c4*std::pow(sdf, (1.0 + c6*sqrtsdf)/sdf) };
  }
  const double c1 = -5.25144;
  const double c2 = -483.195;
  const double c3 =  1.138167;
  const double c4 =  1.017479;
  const double c5 = -0.05471;
  const double c6 =  4.11215;
  return { c1 + z + c2/std::pow(sqrtsdf, c3),
           c4/std::pow(sqrtsdf, c5),
           c4*std::pow(sdf, (1.0 + c6*sqrtsdf)/sdf) };
}

// sdf_all = {AW, SB, SW, PL} SDFs; only the other species' entries are used
inline DenCoef DenCoef_spp(int s, double sdf, double psi, const double* sdf_all) {
  switch (s) {
    case AW: return AW_DenCoef(sdf, psi);
    case SB: return SB_DenCoef(sdf, psi);
    case SW: return SW_DenCoef(sdf, psi, sdf_all[AW]);
    default: return PL_DenCoef(sdf, psi, sdf_all[AW], sdf_all[SW], sdf_all[SB]);
  }
}

inline double Den_kk(int s, const DenCoef& c, double si_bh, double age) {
  double x = (s == AW) ? si_bh : std::log(si_bh);
  return 1.0 + std::exp(c.b1 + c.b2*x + c.b3*std::log(age + 1.0));
}

// Density (stems/ha) at age from SDF
inline double Density(int s, double age, double sdf, double si_bh, double psi, const double* sdf_all) {
  DenCoef c = DenCoef_spp(s, sdf, psi, sdf_all);
  return sdf*Den_kk(s, c, si_bh, 50.0)/Den_kk(s, c, si_bh, age);
}

// SDF (density at age 50) from an observed density, hill climbing from 1000
inline double SolveSDF(int s, double den, double age, double si_bh, double psi,
                       const double* sdf_all, int iter, double tol, bool& converged) {
  double sdf0 = 1000.0, sdf1 = 0.0;
  converged = false;
  for (int it = 0; it < iter; it++) {
    DenCoef c = DenCoef_spp(s, sdf0, psi, sdf_all);
    sdf1 = den*Den_kk(s, c, si_bh, age)/Den_kk(s, c, si_bh, 50.0);
    sdf0 = (sdf0 + sdf1)/2.0;
    if (std::abs(sdf0 - sdf1) < tol) {
      converged = true;
      break;
    }
  }
  return sdf1;
}

// ---- Basal area increment -------------------------------------------------
// Annual BA increment (m2/ha) at breast height age x
//   sc   = BH-age based species composition
//   den0 = initial density (AW at BH age 0, others at total age 0)
//   ba   = current BA of the species
//   sdf1 = {AW, SB, SW, PL} competitor SDFs for SW and PL, already set to 0
//          for species not yet at breast height (caller's job)
// At x = 0 and ba = 0 these reduce to the SAS BAinc0 equations (examp13).

inline double AW_BAinc(double x, double sc, double den0, double si_bh, double psi, double ba) {
  double ld = std::log(1.0 + den0*std::sqrt(1.0 + x));
  if (psi > 0) {
    const double a1 =  0.055906;
    const double a2 =  0.054216;
    const double a3 =  1.240035;
    const double a4 = -0.05515;
    const double a5 =  0.218326;
    const double a6 =  0.749594;
    return 1e-4*a1*x*x*std::exp(-a2*x)*std::pow(sc, a5)*ld*ld*si_bh*std::pow(psi, a6)/
           (std::pow(1.0 + ba, a3)*(1.0 + std::exp(1.0 - std::log(1.0 + sc*sc)/2.0)))
           + a4*std::log(0.01 + x/10.0);
  }
  const double a1 =  0.751313;
  const double a2 =  0.018847;
  const double a3 =  1.143762;
  const double a4 = -0.03475;
  const double a5 =  0.835189;
  return 1e-4*a1*x*x*std::exp(-a2*std::pow(x, 0.5 + a1))*std::pow(sc, a5)*ld*ld*si_bh/
         (std::pow(1.0 + ba, a3)*(1.0 + std::exp(1.0 - std::log(1.0 + sc*sc)/2.0)))
         + a4*std::log(0.01 + x/10.0);
}

inline double SB_BAinc(double x, double sc, double den0, double si_bh, double psi, double ba) {
  if (psi > 0) {
    const double a1 = 0.333889;
    const double a2 = 0.047669;
    const double a3 = 0.03635;
    const double a4 = 0.484915;
    return 1e-4*a1*std::exp(-a2*x)*std::pow(sc, a3)*std::pow(x, std::sqrt(a1) + a2)
           *(1.0 + std::sqrt(den0)*std::sqrt(1.0 + x))*std::exp(-(den0/4.0)/10000.0)
           *std::log(1.0 + si_bh)*std::pow(psi, a4)/std::pow(1.0 + ba, a2);
  }
  const double a1 = 0.966285;
  const double a2 = 0.056315;
  const double a3 = 0.17191;
  return 1e-4*a1*std::exp(-a2*x)*std::pow(sc, a3)*std::pow(x, std::sqrt(a1) + a2)
         *(1.0 + std::sqrt(den0)*std::sqrt(1.0 + x))*std::exp(-(den0/4.0)/10000.0)
         *std::log(1.0 + si_bh)/std::pow(1.0 + ba, a2);
}

inline double SW_BAinc(double x, double sc, double den0, double si_bh, double psi, double ba,
                       const double* sdf1) {
  bool p = psi > 0;
  double a1, a2, a3, a4, a5, a6, a7 = 0.0;
  if (p) {
    a1 =  0.041072; a2 = 0.067012; a3 = -0.00163; a4 = 5.701986;
    a5 =  2.222453; a6 = 0.001126; a7 =  0.368267;
  } else {
    a1 =  0.089153; a2 = 0.072171; a3 = -0.11483; a4 = 5.839408;
    a5 =  1.753002; a6 = 0.239521;
  }
  double k = a4*((sdf1[AW] > 0) ? std::log(1.0 + sdf1[AW]/10000.0) : 0.0)
           + a5*((sdf1[PL] > 0) ? std::log(1.0 + sdf1[PL]/10000.0) : 0.0)
           +    ((sdf1[SB] > 0) ? std::log(1.0 + sdf1[SB]/10000.0) : 0.0);
  double ld  = std::log(1.0 + den0*std::sqrt(1.0 + x));
  double num = 1e-4*a1*(a2 + x)*(a2 + x)*std::pow(1.0 + x, std::sqrt(a1) + a2 - a3)*std::exp(-a2*x)
               *std::sqrt(si_bh)*ld*ld*std::pow(sc, a6)*std::exp(-(den0/10.0)/10000.0);
  if (p) num *= std::pow(psi, a7);
  return num/(1.0 + std::exp(1.0 + k + std::log(1.0 + std::sqrt(den0)/10000.0)/2.0 + a3*std::log(1.0 + ba)));
}

inline double PL_BAinc(double x, double sc, double den0, double si_bh, double psi, double ba,
                       const double* sdf1) {
  bool p = psi > 0;
  double a1, a2, a3, a4, a5, a6;
  if (p) {
    a1 = 0.067391; a2 = 0.054113; a3 = 0.213767;
    a4 = 0.977498; a5 = -0.00452; a6 = 1.080421;
  } else {
    a1 = 3.923984; a2 = 0.05752;  a3 = 0.560402;
    a4 = 0.672506; a5 = -0.00358; a6 = 0.775765;
  }
  double k = ((sdf1[AW] > 0) ? std::log(1.0 + sdf1[AW]/1000.0)     : 0.0)
           + ((sdf1[SW] > 0) ? std::log(1.0 + sdf1[SW]/1000.0)/2.0 : 0.0)
           + ((sdf1[SB] > 0) ? std::log(1.0 + sdf1[SB]/1000.0)/2.0 : 0.0);
  double num = 1e-4*a1*x*std::exp(-a2*x)*(1.0 + a3 + std::pow(si_bh, a6))*std::sqrt(den0)
               *std::exp(-(den0/3.0)/10000.0)*(1.0 + std::log(1.0 + x)/2.0);
  if (p) num *= psi;
  return num/(1.0 + std::exp(std::log(1.0 + (den0/3.0)/10000.0) + k/2.0 - a3*std::sqrt(sc) + a4*std::log(1.0 + ba)))
         + a5*std::log(0.01 + x/10.0);
}

inline double BAinc(int s, double x, double sc, double den0, double si_bh, double psi, double ba,
                    const double* sdf1) {
  switch (s) {
    case AW: return AW_BAinc(x, sc, den0, si_bh, psi, ba);
    case SB: return SB_BAinc(x, sc, den0, si_bh, psi, ba);
    case SW: return SW_BAinc(x, sc, den0, si_bh, psi, ba, sdf1);
    case PL: return PL_BAinc(x, sc, den0, si_bh, psi, ba, sdf1);
  }
  return NaN;
}

// ---- Volume ----------------------------------------------------------------

// Total volume (m3/ha) from BA and top height
inline double TotVol(int s, double ba, double ht) {
  switch (s) {
    case AW: return 0.248718*std::pow(ba, 0.98568)*std::pow(ht, 0.857278)*std::exp(1.0 - 24.9961/(ht*ht + 1.0));
    case SB: return 0.48628 *std::pow(ba, 0.982962)*std::pow(ht, 0.910603);
    case SW: return 0.41104 *std::pow(ba, 0.983108)*std::pow(ht, 0.971061);
    case PL: return 0.194086*std::pow(ba, 0.988276)*std::pow(ht, 0.949346)*std::exp(1.0 - 3.39036/(ht*ht + 1.0));
  }
  return NaN;
}

struct MerchResult { double mv, mn; };

// Merchantable volume and merchantable density for one utilization standard
//   tv = total volume, qmd = sqrt(BA/(N/10000)) (SAS k_<spp>)
//   sh = stump height (m, <= 0.3), sdob = stump DOB (cm), tdib = top DIB (cm)
inline MerchResult Merch(int s, double tv, double qmd, double ht, double bhage, double den,
                         double sh, double sdob, double tdib) {
  // merchantable volume ratio b0-b6, stump DOB height correction a-c,
  // stump diameter to DBH k4-k6
  static const double MV[4][7] = {
    {0.993673, 923.5825, -3.96171, 3.366144, 0.316236, 0.968953, -1.61247},
    {0.98152,  0.678011, -1.10256, 4.148139, 0.511391, 1.484988, -3.26425},
    {0.996262, 7.021736, -1.77615, 1.91562,  0.4111,   1.024803, -0.80121},
    {0.989889, 1.055091, -0.19072, 4.915593, 0.42574,  1.006379, -4.87808}};
  static const double RT[4][3] = {
    {0.009628, 10.82514, 1.068204},
    {0.030241, 29.10004, 1.873577},
    {0.01939,  14.58411, 1.286149},
    {0.013081, 13.01001, 1.213286}};
  static const double KD[4][3] = {
    {-0.33057, 1.10664, 0.000986},
    { 0.53613, 0.99631, 0.004617},
    {-0.56778, 1.14215, 0.000429},
    {-0.24529, 1.08842, 0.000993}};

  const double* b = MV[s];
  double mv = tv*std::pow(qmd, b[0])/(b[1]*std::pow(ht, b[2])*std::pow(sdob, b[3])*std::pow(sh, b[4])
                                      *std::pow(tdib, b[5])*std::pow(qmd, b[6]) + qmd);

  double stump03;   // stump DOB at 0.3 m
  if (sh == 0.3)     stump03 = sdob;
  else if (sh < 0.3) stump03 = sdob/(RT[s][0]*std::exp(RT[s][1]*std::pow(0.3 - sh, RT[s][2])) + 1.0);
  else               stump03 = NaN;

  const double k4 = KD[s][0], k5 = KD[s][1], k6 = KD[s][2];
  double dbh = (-k5 + std::sqrt(k5*k5 - 4.0*k6*(k4 - stump03)))/(2.0*k6);

  // Cumulative % of stems below the merchantable DBH
  double cum;
  if (s == AW) {
    cum = 110.37/(1.0 + std::exp(8.9504 - 4.2732*std::log(dbh + 1.0) + 1.3280*std::log(bhage + 1.0)
                                 + 0.1255*ht - 0.6410*std::log(1.0 + den)));
    if (cum > 100) cum = 100;
  } else {
    // Weibull: a0, b0-b4, c0-c2
    static const double WB[3][9] = {
      {104.09, 3.8437,  0.5145, -0.1358,  0.2201, -0.01063, 0.5013, 0.03470, 0.2354},   // SB
      {102.99, 4.0921,  0.6369, -0.2063,  0.2555, -0.01816, 0.3233, 0.02614, 0.2457},   // SW
      {100.0,  14.4907, 0.7772, -0.2370, -0.1114,  0.03872, 1.5357, 0.02706, 0.07775}}; // PL
    const double* w = WB[s - 1];
    double scale = w[1]*std::pow(ht - 1.3, w[2])*std::pow(den, w[3])*std::pow(bhage, w[4]) + w[5]*bhage;
    double shape = w[6]*std::exp(w[7]*(ht - 1.3))*std::pow(bhage, w[8]);
    cum = w[0]*(1.0 - std::exp(-std::pow(dbh/scale, shape)));
    if (s != PL && cum > 100) cum = 100;   // SAS does not cap PL
  }

  return {mv, den*(100.0 - cum)/100.0};
}

} // namespace gypsy

#endif
