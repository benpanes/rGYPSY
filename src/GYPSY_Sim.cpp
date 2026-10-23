#include <Rcpp.h>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// ---------------------------------------------------------------------------
// Full GYPSY model logic
//
// Includes several component functions, with three exports:
//   GYPSY_Height: Solves for height using SI and age (will only accept total age).
//   GYPSY_PSI: Solves for percent stocking index
//       (Might have to consult with Robert about what these should actually be)
//   GYPSY_Sim: Full GYPSY simulation
//
// ---------------------------------------------------------------------------

namespace { // Global variables
  constexpr double sq50    = 7.0710678118654752;  // sqrt(50)
  constexpr double log51   = 3.9318256327243257;  // log(50 + 1)
  constexpr double log2501 = 7.824445930877619;   // log(50^2 + 1), SW only
}

// ---- Percent Stocking Index -----------------------------------------------

// Single return, no struct needed

// ---- AW -------------------------------------------------------------------
static inline double AW_PSI(double age_t,
                            double ps,
                            double si_b) {
  const double b1 = -14.9555;
  const double b2 = 0.45162;
  const double b3 = 1.59953;

  double lsib  = std::log(si_b);
  double lage1 = std::log(age_t+1.0);

  double psi = std::min(ps*(1.0+std::exp(b1 + b2*lage1*lage1 + b3*lsib))/
                           (1.0+std::exp(b1 + b2*log51*log51 + b3*lsib)),
                        100.0);

  return psi;
}

// ---- SB -------------------------------------------------------------------
static inline double SB_PSI(double age_t,
                            double ps,
                            double si_b) {
  const double b1 = -16.7154;
  const double b2 = 0.369843;
  const double b3 = 2.808537;

  double lsib  = std::log(si_b);
  double lage1 = std::log(age_t+1.0);

  double psi = std::min(ps*(1.0+std::exp(b1 + b2*lage1*lage1 + b3*lsib))/
                           (1.0+std::exp(b1 + b2*log51*log51 + b3*lsib)),
                        100.0);

  return psi;
}

// ---- SW -------------------------------------------------------------------
static inline double SW_PSI(double age_t,
                            double ps,
                            double psi_aw,
                            double si_b) {
  const double b1 = -21.8721;
  const double b2 = 0.395767;
  const double b3 = 4.128582;

  double lsib  = std::log(si_b);
  double lage1 = std::log(age_t+1.0);
  double z1    = psi_aw/b3/50.0; // Should this handle cases where AW PSI is NA?

  double psi = std::min(ps*(1.0+std::exp(b1 + b2*lage1*lage1 + b3*lsib + z1))/
                           (1.0+std::exp(b1 + b2*log51*log51 + b3*lsib + z1)),
                        100.0); // Returns 100 if result exceeds

  return psi;
}

// ---- PL -------------------------------------------------------------------
static inline double PL_PSI(double age_t,
                            double ps,
                            double psi_aw,
                            double psi_sb,
                            double psi_sw,
                            double si_b) {
  const double b1 = -12.2731;
  const double b2 = 0.23981;
  const double b3 = 2.17610;
  const double b4 = 0.02492;
  const double b5 = 0.08990;
  const double b6 = 0.17058;

  double lsib  = std::log(si_b);
  double lage1 = std::log(age_t+1.0);
  double z1    = b4*psi_aw;
  double z2    = b5*sqrt(psi_sb);
  double z3    = b6*sqrt(psi_sw);

  double psi = std::min(ps*(1.0+std::exp(b1 + b2*lage1*lage1 + b3*lsib + z1 + z2 + z3))/
                           (1.0+std::exp(b1 + b2*log51*log51 + b3*lsib + z1 + z2 + z3)),
                        100.0);

  return psi;
}

// ---- Stand Density Factor -------------------------------------------------
// Two function styles: with or without % stocking
// Entries need to be zero-filled (either no PS or NA species)
// Initial guess of 1000

struct SDFResult { double sdf1, sdf0; };

// ---- AW -------------------------------------------------------------------
static inline SDFResult AW_SDF(double bhage, // Only AW uses BH age for some reason
                               double den,
                               double sibh,
                               double psi_aw,
                               double sdf0) {

  double sqrtsdf0 = sqrt(sdf0);

  double b1, b2, b3;
  if(psi_aw>0){
    const double c0 = 0.71123;
    const double c1 = 7.015511;

    double b1 = -(1.0/std::sqrt(sdf0/1000.0)+std::sqrt(51.0/(sqrtsdf0*log51)) + c0*std::log(1.0+psi_aw/sqrtsdf0)/c1)*log51;
    double b2 = c0*std::pow(sqrtsdf0,1.0/sdf0)/4.0;
    double b3 = (1.0+c0)*std::pow(sdf0,(c1+std::log(sdf0) + std::log(1.0+psi_aw)/sdf0)/sdf0);
  }
  else         {
    const double c0 = 0.717966;
    const double c1 = 6.67468;

    double b1 = -(1.0/sqrt(sdf0/1000.0)+std::sqrt(51.0/(sqrtsdf0*log51)))*log51;
    double b2 = c0*std::pow(sqrtsdf0,1/sdf0)/4.0;
    double b3 = (1.0+c0)*std::pow(sdf0,(c1+std::log(sdf0)/sdf0));
  }

  double kk1   = 1.0+std::exp(b1+b2*sibh+b3*log51);
  double kk2   = 1.0+std::exp(b1+b2*sibh+b3*std::log(bhage+1));
  double sdf1  = den*kk2/kk1;
  double sdf0n = (sdf1+sdf0n)/2.0;

  return {sdf1, sdf0n};
}

// ---- SB -------------------------------------------------------------------
static inline SDFResult SB_SDF(double age,
                               double den,
                               double sibh,
                               double psi_sb,
                               double sdf0) {

  double b1, b2, b3;
  if(psi_sb>0){
    const double c1 = -25.7666;
    const double c2 = 0.192351;
    const double c3 = 2.718748;

    double b1 = c1/std::pow(sqrt(sdf0/1000.0)+log51,c2)-std::log(1.0+psi_sb)/(1.0+c3);
    double b2 = c3;
    double b3 = c3*std::pow(sdf0,1.0/sdf0);
  }
  else         {
    const double c1 = -26.3836;
    const double c2 = 0.166483;
    const double c3 = 2.738569;

    double b1 = c1/std::pow(sqrt(sdf0/1000.0)+log51,c2);
    double b2 = c3;
    double b3 = c3*std::pow(sdf0,1.0/sdf0);
  }

  double kk1   = 1.0+std::exp(b1+b2*std::log(sibh)+b3*log51);
  double kk2   = 1.0+std::exp(b1+b2*sibh+b3*std::log(age+1.0));
  double sdf1  = den*kk2/kk1;
  double sdf0n = (sdf1+sdf0n)/2.0;

  return {sdf1, sdf0n};
}

// ---- SW -------------------------------------------------------------------
static inline SDFResult SW_SDF(double age,
                               double den,
                               double sibh,
                               double psi_sw,
                               double sdf_aw,
                               double sdf0) {
  double z1;
  if(sdf_aw>0) z1 = std::sqrt(1+sdf_aw/1000);
  else         z1 = 0;

  double b1, b2, b3;
  if(psi_sw>0){
    const double c1 = -254.391;
    const double c2 = 1.146709;
    const double c3 = 2.202149;

    double b1 = c1/std::pow(std::log(sdf0)+log51,c2)+z1+std::log(1.0+sdf0)/10.0;
    double b2 = c3;
    double b3 = c3*std::pow(sdf0,1.0/sdf0)*std::pow(std::sqrt(1.0+psi_sw),-1.0/10.0);
  }
  else         {
    const double c1 = -231.617;
    const double c2 = 1.176995;
    const double c3 = 1.733601;

    double b1 = c1/std::pow(std::log(sdf0)+log51,c2) + z1;
    double b2 = c3;
    double b3 = c3*std::pow(sdf0,1.0/sdf0);
  }

  double kk1   = 1.0+std::exp(b1+b2*std::log(sibh)+b3*log51);
  double kk2   = 1.0+std::exp(b1+b2*sibh+b3*std::log(age+1.0));
  double sdf1  = den*kk2/kk1;
  double sdf0n = (sdf1+sdf0n)/2.0;

  return {sdf1, sdf0n};
}

// ---- PL -------------------------------------------------------------------
static inline SDFResult PL_SDF(double age,
                               double den,
                               double sibh,
                               double psi_pl,
                               double sdf_aw,
                               double sdf_sw,
                               double sdf_sb,
                               double sdf0) {
  double sqrtsdf0 = std::sqrt(sdf0);

  double z1;
  if(sdf_aw>0) z1 = sdf_aw/1000.0/2.0;
  else         z1 = 0;
  double z2;
  if(sdf_sw>0) z2 = sdf_sw/1000.0/3.0;
  else         z2 = 0;
  double z3;
  if(sdf_sb>0) z3 = sdf_sb/1000.0/4.0;
  else         z3 = 0;

  double b1, b2, b3;
  if(psi_pl>0){
    const double c1 = -3.09711;
    const double c2 = -512.394;
    const double c3 = 1.14981;
    const double c4 = 1.027108;
    const double c5 = -0.05479;
    const double c6 = 4.12191;

    double b1 = c1+z1+z2+z3+c2/std::pow(sqrtsdf0,c3)-std::log(1.0+psi_pl)/(1.0+c4);
    double b2 = c4/std::pow(sqrtsdf0,c5);
    double b3 = c4*std::pow(sdf0,(1.0+c6*sqrtsdf0)/sdf0);
  }
  else         {
    const double c1 = -5.25144;
    const double c2 = -483.195;
    const double c3 = 1.138167;
    const double c4 = 1.017479;
    const double c5 = -0.05471;
    const double c6 = 4.11215;

    double b1 = c1+z1+z2+z3+c2/std::pow(sqrtsdf0,c3);
    double b2 = c4/std::pow(sqrtsdf0,c5);
    double b3 = c4*std::pow(sdf0,(1+c6*sqrtsdf0)/sdf0);
  }

  double kk1   = 1.0+std::exp(b1+b2*log(sibh)+b3*log51);
  double kk2   = 1.0+std::exp(b1+b2*sibh+b3*log(age+1.0));
  double sdf1  = den*kk2/kk1;
  double sdf0n = (sdf1+sdf0n)/2.0;

  return {sdf1, sdf0n};
}

// ---- SDF hill climber -----------------------------------------------------
template <typename Func>
double solve_sdf(Func   f,
                 int    iter,
                 double tol,
                 int    row_id) {

  double sdf0 = 1000.0, sdf1 = 0.0, diff = 1000.0;
  bool converged = false;

  for (int it = 0; it < iter; it++) {
    auto res = f(sdf0);

    sdf0 = res.sdf0;
    sdf1 = res.sdf1;
    diff = std::abs(sdf0 - sdf1);

    if (diff < tol) {
      converged = true;
      break;
    }
  }

  if (!converged) {
    Rcpp::warning("Row " + std::to_string(row_id + 1) +
      " did not converge. Final diff was " +
      std::to_string(diff));
  }

  return sdf1;
}

// ---- Top Height Calculation -----------------------------------------------
// Height from total SI and total age
// ---- AW -------------------------------------------------------------------
static inline double AW_Ht(double si_t,
                           double age_t) {
  const double b1 =  9.908888;
  const double b2 = -3.92451;
  const double b3 = -0.32778;
  const double b4 =  0.134376;

  double lsit = std::log(si_t);

  // Numerator and denominator set as separate variables for readability
  double x10 = 1.0 + std::exp(b1 + b2*std::sqrt(log51)               + b3*lsit*lsit + b4*sq50);
  double x20 = 1.0 + std::exp(b1 + b2*std::sqrt(std::log(age_t+1.0)) + b3*lsit*lsit + b4*sq50);

  double ht  = si_t*x10/x20;

  return ht;
}

// ---- SB -------------------------------------------------------------------
static inline double SB_Ht(double si_t,
                           double age_t) {
  const double b1 =  14.56236;
  const double b2 =  -6.04705;
  const double b3 =  -1.53715;
  const double b4 =   0.240174;

  double lsit = std::log(si_t);

  double x10  = 1.0 + exp(b1 + b2*std::sqrt(log51)               + b3*lsit + b4*sq50);
  double x20  = 1.0 + exp(b1 + b2*std::sqrt(std::log(age_t+1.0)) + b3*lsit + b4*sq50);

  double ht   = si_t*x10/x20;

  return ht;
}

// ---- SW -------------------------------------------------------------------
static inline double SW_Ht(double si_t,
                           double age_t) {
  const double b1 =  12.14943;
  const double b2 =  -3.77051;
  const double b3 =  -0.28534;
  const double b4 =   0.165483;

  double lsit = std::log(si_t);

  double x10  = 1.0 + std::exp(b1 + b2*std::sqrt(log2501)              + b3*lsit*lsit + b4*sq50);
  double x20  = 1.0 + exp(b1 + b2*std::sqrt(std::log(age_t*age_t+1.0)) + b3*lsit*lsit + b4*sq50);

  double ht   = si_t*x10/x20;

  return ht;
}

// ---- PL -------------------------------------------------------------------
static inline double PL_Ht(double si_t,
                           double age_t) {
  const double b1 =  12.84571;
  const double b2 =  -5.73936;
  const double b3 =  -0.91312;
  const double b4 =   0.150668;

  double lsit = std::log(si_t);

  double x10  = 1.0 + exp(b1 + b2*std::sqrt(log51)               + b3*lsit + b4*sq50);
  double x20  = 1.0 + exp(b1 + b2*std::sqrt(std::log(age_t+1.0)) + b3*lsit + b4*sq50);

  double ht   = si_t*x10/x20;

  return ht;
}

// ---- Exported PSI function ------------------------------------------------
// [[Rcpp::export]]
Rcpp::List GYPSY_PSI_cpp(CharacterVector stand,

                         NumericVector   aw_ps,
                         NumericVector   aw_age_t,
                         NumericVector   aw_si_t,

                         NumericVector   sb_ps,
                         NumericVector   sb_age_t,
                         NumericVector   sb_si_t,

                         NumericVector   sw_ps,
                         NumericVector   sw_age_t,
                         NumericVector   sw_si_t,

                         NumericVector   pl_ps,
                         NumericVector   pl_age_t,
                         NumericVector   pl_si_t) {
  int n = stand.size();

  NumericVector aw_res(n), sb_res(n), sw_res(n), pl_res(n);
  for (int i = 0; i < n; i++) {
    double aw_psI    = aw_ps[i];
    double aw_age_tI = aw_age_t[i];
    double aw_si_tI  = aw_si_t[i];

    double sb_psI    = sb_ps[i];
    double sb_age_tI = sb_age_t[i];
    double sb_si_tI  = sb_si_t[i];

    double sw_psI    = sw_ps[i];
    double sw_age_tI = sw_age_t[i];
    double sw_si_tI  = sw_si_t[i];

    double pl_psI    = pl_ps[i];
    double pl_age_tI = pl_age_t[i];
    double pl_si_tI  = pl_si_t[i];

    double aw_psi = (aw_psI > 0 && aw_age_tI > 0)
      ? AW_PSI(aw_age_tI, aw_psI, aw_si_tI) : 0.0; // ? = ifelse equivalent
    double sb_psi = (sb_psI > 0 && aw_age_tI > 0)
      ? SB_PSI(sb_age_tI, sb_psI, sb_si_tI): 0.0;
    double sw_psi = (sw_psI > 0 && sw_age_tI > 0)
      ? SW_PSI(sw_age_tI, sw_psI, aw_psi, sw_si_tI): 0.0;
    double pl_psi = (pl_psI > 0 && pl_age_tI > 0)
      ? PL_PSI(pl_age_tI, pl_psI, aw_psi, sb_psi, sw_psi, pl_si_tI): 0.0;

    aw_res[i] = aw_psi;
    sb_res[i] = sb_psi;
    sw_res[i] = sw_psi;
    pl_res[i] = pl_psi;
  }
  return List::create(
    _["aw_psi"] = aw_res,
    _["sb_psi"] = sb_res,
    _["sw_psi"] = sw_res,
    _["pl_psi"] = pl_res
  );
}

// ---- Exported SDF function ------------------------------------------------
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
  int n = stand.size();

  NumericVector aw_res(n), sb_res(n), sw_res(n), pl_res(n);
  for (int i = 0; i < n; i++) {
    // --- AW ---
    double aw_sdf = (aw_n[i] > 0)
      ? solve_sdf(
          [&](double sdf0) {
            return AW_SDF(aw_age_b[i], aw_n[i], aw_si_b[i], aw_psi[i], sdf0);
          },
          iter, tol, i)
        : 0.0;

    // --- SB ---
    double sb_sdf = (sb_n[i] > 0)
      ? solve_sdf(
          [&](double sdf0) {
            return SB_SDF(sb_age_t[i], sb_n[i], sb_si_b[i], sb_psi[i], sdf0);
          },
          iter, tol, i)
        : 0.0;

    // --- SW
    double sw_sdf = (sw_n[i] > 0)
      ? solve_sdf(
          [&](double sdf0) {
            return SW_SDF(sw_age_t[i], sw_n[i], sw_si_b[i], sw_psi[i], aw_sdf, sdf0);
          },
          iter, tol, i)
        : 0.0;

    // --- PL
    double pl_sdf = (pl_n[i] > 0)
      ? solve_sdf(
          [&](double sdf0) {
            return PL_SDF(pl_age_t[i], pl_n[i], pl_si_b[i], pl_psi[i], aw_sdf, sw_sdf, sb_sdf, sdf0);
          },
          iter, tol, i)
        : 0.0;


    aw_res[i] = aw_sdf;
    sb_res[i] = sb_sdf;
    sw_res[i] = sw_sdf;
    pl_res[i] = pl_sdf;
  }
  return List::create(
    _["aw_sdf"] = aw_res,
    _["sb_sdf"] = sb_res,
    _["sw_sdf"] = sw_res,
    _["pl_sdf"] = pl_res
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
    std::string sppI    = as<std::string>(spp[i]);
    double      si_tI   = si_t[i];
    double      age_tI  = age_t[i];

    double ht;
    if        (sppI == "AW") {
      ht = AW_Ht(si_tI, age_tI);
    } else if (sppI == "SB") {
      ht = SB_Ht(si_tI, age_tI);
    } else if (sppI == "SW") {
      ht = SW_Ht(si_tI, age_tI);
    } else if (sppI == "PL") {
      ht = PL_Ht(si_tI, age_tI);
    }
    res[i] = ht;
  }
  return res;
}
