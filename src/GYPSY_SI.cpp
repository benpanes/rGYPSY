#include <Rcpp.h>
#include <cmath>
#include <string>
using namespace Rcpp;

// ---------------------------------------------------------------------------
// Species-specific SI solvers
//
// Each solver takes the current si0 guess and returns a SIResult:
//   si0   = (si0_old + si1) / 2  - averaged, becomes next iteration input
//   si1   = raw si1 before averaging - used for convergence check
//   si_bh = site index at breast height age
//   Y2BH  = years to breast height
//   Y2SH  = years to stump height
//
// Convergence check: abs(si0_averaged - si1_raw)
//
// Adjusts ages to total ages within each equation
// ---------------------------------------------------------------------------

namespace {                                       // Global variables
  constexpr double sq50    = 7.0710678118654752;  // sqrt(50)
  constexpr double log51   = 3.9318256327243257;  // log(50 + 1)
  constexpr double log2501 = 7.824445930877619;   // log(50^2 + 1), SW only
}

struct SIResult { double si0, si1, si_bh, Y2BH, Y2SH; };

// ---- AW -------------------------------------------------------------------
static inline SIResult AW_SI(double          age,
                             std::string     at,
                             double          ht,
                             double          si0,
                             double          bh,
                             double          sh) {
  const double b1 =  9.908888;
  const double b2 = -3.92451;
  const double b3 = -0.32778;
  const double b4 =  0.134376;

  double lsi0 = log(si0);
  double k1   = exp(b1 + b2*sqrt(log51) + b3*lsi0*lsi0 + b4*sq50);
  double k2   = pow(si0, b3*lsi0);
  double k3b  = (si0*(1.0+k1)/bh - 1.0) / (exp(b1)*exp(b4*sq50)*k2);
  double k3s  = (si0*(1.0+k1)/sh - 1.0) / (exp(b1)*exp(b4*sq50)*k2);
  double Y2BH = exp(pow(log(k3b)/b2, 2)) - 1.0;
  double Y2SH = exp(pow(log(k3s)/b2, 2)) - 1.0;

  double tage;
  if      (at == "T") tage = age;
  else if (at == "B") tage = age+Y2BH;
  else                tage = age+Y2SH;

  double x10  = 1.0 + exp(b1 + b2*sqrt(log(tage+1.0)) + b3*lsi0*lsi0 + b4*sq50);
  double x20  = 1.0 + exp(b1 + b2*sqrt(log51)        + b3*lsi0*lsi0 + b4*sq50);
  double si1  = ht * x10 / x20;
  double si0n = (si0 + si1) / 2.0;

  double lsi1  = log(si1);
  double si_bh = si1 * (1.0 + exp(b1 + b2*sqrt(log51)              + b3*lsi1*lsi1 + b4*sq50)) /
                       (1.0 + exp(b1 + b2*sqrt(log(50.0+Y2BH+1.0)) + b3*lsi1*lsi1 + b4*sq50));
  return {si0n, si1, si_bh, Y2BH, Y2SH};
}

// ---- SB -------------------------------------------------------------------
static inline SIResult SB_SI(double age,
                             std::string at,
                             double ht,
                             double si0,
                             double bh,
                             double sh) {
  const double b1 =  14.56236;
  const double b2 =  -6.04705;
  const double b3 =  -1.53715;
  const double b4 =   0.240174;

  double lsi0 = log(si0);
  double k1   = exp(b1 + b2*sqrt(log51) + b3*lsi0 + b4*sq50);
  double k2   = pow(si0, b3);
  double k3b   = (si0*(1.0+k1)/bh - 1.0) / (exp(b1)*exp(b4*sq50)*k2);
  double k3s   = (si0*(1.0+k1)/sh - 1.0) / (exp(b1)*exp(b4*sq50)*k2);
  double Y2BH = exp(pow(log(k3b)/b2, 2)) - 1.0;
  double Y2SH = exp(pow(log(k3s)/b2, 2)) - 1.0;

  double tage;
  if      (at == "T") tage = age;
  else if (at == "B") tage = age+Y2BH;
  else                     tage = age+Y2SH;

  double x10  = 1.0 + exp(b1 + b2*sqrt(log(tage+1.0)) + b3*lsi0 + b4*sq50);
  double x20  = 1.0 + exp(b1 + b2*sqrt(log51)        + b3*lsi0 + b4*sq50);
  double si1  = ht * x10 / x20;
  double si0n = (si0 + si1) / 2.0;

  double lsi1  = log(si1);
  double si_bh = si1 * (1.0 + exp(b1 + b2*sqrt(log51)              + b3*lsi1 + b4*sq50)) /
                       (1.0 + exp(b1 + b2*sqrt(log(50.0+Y2BH+1.0)) + b3*lsi1 + b4*sq50));
  return {si0n, si1, si_bh, Y2BH, Y2SH};
}

// ---- SW -------------------------------------------------------------------
static inline SIResult SW_SI(double age,
                             std::string at,
                             double ht,
                             double si0,
                             double bh,
                             double sh) {
  const double b1 =  12.14943;
  const double b2 =  -3.77051;
  const double b3 =  -0.28534;
  const double b4 =   0.165483;

  double lsi0 = log(si0);
  double k1   = exp(b1 + b2*sqrt(log2501)           + b3*lsi0*lsi0 + b4*sq50);
  double k2   = pow(si0, b3*lsi0);
  double k3b   = (si0*(1.0+k1)/bh - 1.0) / (exp(b1)*exp(b4*sq50)*k2);
  double k3s   = (si0*(1.0+k1)/sh - 1.0) / (exp(b1)*exp(b4*sq50)*k2);
  double Y2BH = sqrt(exp(pow(log(k3b)/b2, 2)) - 1.0);
  double Y2SH = sqrt(exp(pow(log(k3s)/b2, 2)) - 1.0);

  double tage;
  if      (at == "T") tage = age;
  else if (at == "B") tage = age+Y2BH;
  else                     tage = age+Y2SH;

  double x10  = 1.0 + exp(b1 + b2*sqrt(log(tage*tage+1.0)) + b3*lsi0*lsi0 + b4*sq50);
  double x20  = 1.0 + exp(b1 + b2*sqrt(log2501)          + b3*lsi0*lsi0 + b4*sq50);
  double si1  = ht * x10 / x20;
  double si0n = (si0 + si1) / 2.0;

  double lsi1   = log(si1);
  double y50    = 50.0 + Y2BH;
  double si_bh  = si1 * (1.0 + exp(b1 + b2*sqrt(log2501)               + b3*lsi1*lsi1 + b4*sq50)) /
                        (1.0 + exp(b1 + b2*sqrt(log(y50*y50 + 1.0))    + b3*lsi1*lsi1 + b4*sq50));
  return {si0n, si1, si_bh, Y2BH, Y2SH};
}

// ---- PL -------------------------------------------------------------------
static inline SIResult PL_SI(double age,
                             std::string at,
                             double ht,
                             double si0,
                             double bh,
                             double sh) {
  const double b1 =  12.84571;
  const double b2 =  -5.73936;
  const double b3 =  -0.91312;
  const double b4 =   0.150668;

  double lsi0 = log(si0);
  double k1   = exp(b1 + b2*sqrt(log51) + b3*lsi0 + b4*sq50);
  double k2   = pow(si0, b3);
  double k3b   = (si0*(1.0+k1)/bh - 1.0) / (exp(b1)*exp(b4*sq50)*k2);
  double k3s   = (si0*(1.0+k1)/sh - 1.0) / (exp(b1)*exp(b4*sq50)*k2);
  double Y2BH = exp(pow(log(k3b)/b2, 2)) - 1.0;
  double Y2SH = exp(pow(log(k3s)/b2, 2)) - 1.0;

  double tage;
  if      (at == "T") tage = age;
  else if (at == "B") tage = age+Y2BH;
  else                tage = age+Y2SH;

  double x10  = 1.0 + exp(b1 + b2*sqrt(log(tage+1.0)) + b3*lsi0 + b4*sq50);
  double x20  = 1.0 + exp(b1 + b2*sqrt(log51)        + b3*lsi0 + b4*sq50);
  double si1  = ht * x10 / x20;
  double si0n = (si0 + si1) / 2.0;

  double lsi1  = log(si1);
  double si_bh = si1 * (1.0 + exp(b1 + b2*sqrt(log51)              + b3*lsi1 + b4*sq50)) /
                       (1.0 + exp(b1 + b2*sqrt(log(50.0+Y2BH+1.0)) + b3*lsi1 + b4*sq50));
  return {si0n, si1, si_bh, Y2BH, Y2SH};
}

// ---------------------------------------------------------------------------
// Main exported function
//
// @param spp       Character vector: "AW", "PL", "SW", or "SB"
// @param ht        Numeric vector: top height (must be > 0.1)
// @param age       Numeric vector: age
// @param agetype   String: type of input age:
//                    "T" = Total age <- default
//                    "B" = Breast height age
//                    "S" = Stump height age
// @param ndx       Integer 1-3 - which output to return:
//                    1 = si1 (SI at total age)
//                    2 = si_bh (SI at BH age) <- default
//                    3 = Y2BH (years to breast height)
//                    4 = Y2SH (years to stump height)
// @param iter      Max convergence iterations (default 1000)
// @param tol       Convergence tolerance (default 1e-8)
// @param bh        Breast height in m (default 1.3)
// @param sh        Stump height in m (default 0.3)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector GYPSY_SI_cpp(CharacterVector spp,
                           NumericVector   ht,
                           NumericVector   age,
                           std::string     agetype = "T",
                           int             ndx  = 2,
                           int             iter = 1000,
                           double          tol  = 1e-8,
                           double          bh   = 1.3,
                           double          sh   = 0.3) {
  int n = spp.size();

  if (ht.size() != n || age.size() != n)
    stop("spp, ht, and age must all be the same length.");
  if (ndx < 1 || ndx > 4)
    stop("ndx must be 1, 2, 3, or 4.");
  if (agetype != "T" && agetype != "B" && agetype != "S")
    stop("agetype must be T, B, or S");

  std::vector<std::string> errors;
  for (int i = 0; i < n; i++) {
    std::string s = as<std::string>(spp[i]);
    if (s != "AW" && s != "PL" && s != "SW" && s != "SB")
      errors.push_back("Row " + std::to_string(i+1) + ": species must be AW, PL, SW, or SB.");
    if (agetype == "T")
      if (ht[i] < 0.1 || NumericVector::is_na(ht[i]))
        errors.push_back("Row " + std::to_string(i+1) + ": ht must be >= 0.1.");
    if (agetype == "S")
      if (ht[i] < sh || NumericVector::is_na(ht[i]))
        errors.push_back("Row " + std::to_string(i+1) + ": ht must be >= stump height.");
    if (agetype == "B")
      if (ht[i] < bh || NumericVector::is_na(ht[i]))
        errors.push_back("Row " + std::to_string(i+1) + ": ht must be >= breast height.");
    if ((s == "AW" || s == "PL") && age[i] <= 7.5)
      errors.push_back("Row " + std::to_string(i+1) + ": age must be > 7.5 for AW or PL.");
    if ((s == "SW" || s == "SB") && age[i] <= 9.5)
      errors.push_back("Row " + std::to_string(i+1) + ": age must be > 9.5 for SW or SB.");
  }
  if (!errors.empty()) {
    std::string msg;
    for (auto& e : errors) msg += e + "\n";
    stop(msg.c_str());
  }

  NumericVector result(n);

  for (int i = 0; i < n; i++) {
    std::string s    = as<std::string>(spp[i]);
    double      th   = ht[i];
    double      ag   = age[i];
    double      si0  = 10.0;    // starting guess
    double      si1  = 0.0, si_bh = 0.0, Y2BH = 0.0, Y2SH = 0.0;
    double      diff = 10.0; // Starting difference

    bool converged = false;

    for (int it = 0; it < iter; it++) {
      SIResult res;
      if      (s == "AW") res = AW_SI(ag, agetype, th, si0, bh, sh);
      else if (s == "SB") res = SB_SI(ag, agetype, th, si0, bh, sh);
      else if (s == "SW") res = SW_SI(ag, agetype, th, si0, bh, sh);
      else                res = PL_SI(ag, agetype, th, si0, bh, sh);

      si0   = res.si0;    // (si0_old + si1) / 2 - next iteration's input
      si1   = res.si1;    // raw si1
      si_bh = res.si_bh;
      Y2BH  = res.Y2BH;
      Y2SH  = res.Y2SH;

      diff = std::abs(si0 - si1);

      if (diff < tol) {
        converged = true;
        break;
      }
    }
    if (!converged){
      Rcpp::warning("Row " + std::to_string(i+1) + " did not converge. Final diff was " + std::to_string(diff)); // This might slow things considerably if several fail, but that shouldn't be the case
    }

    if      (ndx == 1) result[i] = si1;
    else if (ndx == 2) result[i] = si_bh;
    else if (ndx == 3) result[i] = Y2BH;
    else               result[i] = Y2SH;
  }

  return result;
}
