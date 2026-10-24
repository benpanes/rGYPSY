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

namespace {
  constexpr double sq50    = 7.0710678118654752;  // sqrt(50)
  constexpr double log51   = 3.9318256327243257;  // log(50 + 1)
  constexpr double log2501 = 7.824445930877619;  // log(50^2 + 1), SW only
}

struct SIResult { double si0, si1, sibh, age,  y2bh, y2sh; };

// ---- AW -------------------------------------------------------------------
static inline SIResult AW_SI(double       sibh,
                             std::string  sitype, //func = 1, from sibh; func = 2, from sit
                             int          ndx, // 1 = T, 2 = S, 3 = B
                             double       ht,
                             double       si0,
                             double       bh,
                             double       sh) {
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
  double x10 = 1.0 + exp(b1 + b2*sqrt(log(50+Y2BH+1.0)) + b3*lsi0*lsi0 + b4*sq50);
  double x20 = 1.0 + exp(b1 + b2*sqrt(log51)        + b3*lsi0*lsi0 + b4*sq50);

  double si1 = sibh*x10/x20;

  double si0n;
  double sibhn;
  if (sitype == "B") {
    si0n  = (si0 + si1) / 2.0;
    sibhn = sibh;
  }
  else           {
    si0n  = si0;
    sibhn = si0*x20/x10;
  }

  double k3x       = (si0n*(1+k1)/ht-1)/(exp(b1)*exp(b4*sq50)*k2);
  double logk3x_b2 = log(k3x)/b2;
  double tage      = exp(logk3x_b2*logk3x_b2)-1;

  double agen = tage;   // ndx 4-6 return SI or Y2BH/Y2SH, age not used
  if      (ndx == 1)
    agen = tage;
  else if (ndx == 2)
    agen = tage - Y2SH;
  else if (ndx == 3)
    agen = tage - Y2BH;

  return {si0n, si1, sibhn, agen, Y2BH, Y2SH};
}

// ---- SB -------------------------------------------------------------------
static inline SIResult SB_SI(double       sibh,
                             std::string  sitype,
                             int          ndx,
                             double       ht,
                             double       si0,
                             double       bh,
                             double       sh) {
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
  double x10  = 1.0 + exp(b1 + b2*sqrt(log(50+Y2BH+1.0)) + b3*lsi0 + b4*sq50);
  double x20  = 1.0 + exp(b1 + b2*sqrt(log51)            + b3*lsi0 + b4*sq50);

  double si1 = sibh*x10/x20;

  double si0n;
  double sibhn;
  if (sitype == "B") {
    si0n  = (si0 + si1) / 2.0;
    sibhn = sibh;
  }
  else           {
    si0n  = si0;
    sibhn = si0*x20/x10;
  }

  double k3x       = (si0n*(1+k1)/ht-1)/(exp(b1)*exp(b4*sq50)*k2);
  double logk3x_b2 = log(k3x)/b2;
  double tage      = exp(logk3x_b2*logk3x_b2)-1;

  double agen = tage;   // ndx 4-6 return SI or Y2BH/Y2SH, age not used
  if      (ndx == 1)
    agen = tage;
  else if (ndx == 2)
    agen = tage - Y2SH;
  else if (ndx == 3)
    agen = tage - Y2BH;

  return {si0n, si1, sibhn, agen, Y2BH, Y2SH};
}

// ---- SW -------------------------------------------------------------------
static inline SIResult SW_SI(double       sibh,
                             std::string  sitype,
                             int          ndx,
                             double       ht,
                             double       si0,
                             double       bh,
                             double       sh) {
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
  double x10  = 1.0 + exp(b1 + b2*sqrt(log((50+Y2BH)*(50+Y2BH)+1.0)) + b3*lsi0*lsi0 + b4*sq50);
  double x20  = 1.0 + exp(b1 + b2*sqrt(log2501)               + b3*lsi0*lsi0 + b4*sq50);

  double si1 = sibh*x10/x20;

  double si0n;
  double sibhn;
  if (sitype == "B") {
    si0n  = (si0 + si1) / 2.0;
    sibhn = sibh;
  }
  else           {
    si0n  = si0;
    sibhn = si0*x20/x10;
  }

  double k3x       = (si0n*(1+k1)/ht-1)/(exp(b1)*exp(b4*sq50)*k2);
  double logk3x_b2 = log(k3x)/b2;
  double tage      = sqrt(exp(logk3x_b2*logk3x_b2)-1);

  double agen = tage;   // ndx 4-6 return SI or Y2BH/Y2SH, age not used
  if      (ndx == 1)
    agen = tage;
  else if (ndx == 2)
    agen = tage - Y2SH;
  else if (ndx == 3)
    agen = tage - Y2BH;

  return {si0n, si1, sibhn, agen, Y2BH, Y2SH};
}

// ---- PL -------------------------------------------------------------------
static inline SIResult PL_SI(double       sibh,
                             std::string  sitype,
                             int          ndx,
                             double       ht,
                             double       si0,
                             double       bh,
                             double       sh) {
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
  double x10  = 1.0 + exp(b1 + b2*sqrt(log(50+Y2BH+1.0)) + b3*lsi0 + b4*sq50);
  double x20  = 1.0 + exp(b1 + b2*sqrt(log51)            + b3*lsi0 + b4*sq50);

  double si1 = sibh*x10/x20;

  double si0n;
  double sibhn;
  if (sitype == "B") {
    si0n  = (si0 + si1) / 2.0;
    sibhn = sibh;
  }
  else           {
    si0n  = si0;
    sibhn = si0*x20/x10;
  }

  double k3x       = (si0n*(1+k1)/ht-1)/(exp(b1)*exp(b4*sq50)*k2);
  double logk3x_b2 = log(k3x)/b2;
  double tage      = exp(logk3x_b2*logk3x_b2)-1;

  double agen = tage;   // ndx 4-6 return SI or Y2BH/Y2SH, age not used
  if      (ndx == 1)
    agen = tage;
  else if (ndx == 2)
    agen = tage - Y2SH;
  else if (ndx == 3)
    agen = tage - Y2BH;

  return {si0n, si1, sibhn, agen, Y2BH, Y2SH};
}

// ---------------------------------------------------------------------------
// Main exported function
//
// @param spp    Character vector: "AW", "PL", "SW", or "SB"
// @param ht  Numeric vector: top height (must be > 0.1)
// @param age    Numeric vector: stand age
// @param ndx    Integer 1-3 - which output to return:
//                 1 = si1 (SI at total age)
//                 2 = si_bh (SI at BH age) <- default
//                 3 = Y2BH (years to breast height)
//                 4 = Y2SH (years to stump height)
// @param iter   Max convergence iterations (default 1000)
// @param tol    Convergence tolerance (default 1e-8)
// @param bh     Breast height in m (default 1.3)
// @param sh     Stump height in m (default 0.3)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector GYPSY_Age_cpp(CharacterVector spp,
                            NumericVector   si,
                            Nullable<NumericVector> h = R_NilValue,
                            std::string     sitype = "B",
                            int             ndx  = 1, // ndx 1 = total age, 2 = stump age, 3 = breast height age, 4 = si conversion
                            int             iter = 1000,
                            double          tol  = 1e-8,
                            double          bh   = 1.3,
                            double          sh   = 0.3) {
  int n = spp.size();

  NumericVector ht;
  if (h.isNotNull()) {
    ht = NumericVector(h);
  } else {
    // fill missing height with NA
    ht = NumericVector(spp.size(), NA_REAL);
  }

  if (ht.size() != n || si.size() != n)
      stop("input vectors must all be the same length.");
  if (ndx < 1 || ndx > 6)
    stop("ndx must be between 1-6.");
  if (sitype != "B" && sitype != "T")
    stop("site index type must be B, or T");

  std::vector<std::string> errors;
  for (int i = 0; i < n; i++) {
    std::string s = as<std::string>(spp[i]);
    if (s != "AW" && s != "PL" && s != "SW" && s != "SB")
      errors.push_back("Row " + std::to_string(i+1) + ": species must be AW, PL, SW, or SB.");
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
    double      sibh = 0.0, si1 = 0.0, age = 0.0, y2bh = 0.0, y2sh = 0.0;
    if (sitype == "B"){
                  sibh = si[i];
      double      si0  = si[i];    // Use BH index for starting guess
      double      diff = 10.0;     // recalculated at end of loop

      bool converged = false;

      for (int it = 0; it < iter; it++) {
        SIResult res;
        if      (s == "AW") res = AW_SI(sibh, sitype, ndx, th, si0, bh, sh);
        else if (s == "SB") res = SB_SI(sibh, sitype, ndx, th, si0, bh, sh);
        else if (s == "SW") res = SW_SI(sibh, sitype, ndx, th, si0, bh, sh);
        else                res = PL_SI(sibh, sitype, ndx, th, si0, bh, sh);

        si0   = res.si0;    // (si0_old + si1) / 2 - next iteration's input
        si1   = res.si1;    // raw si1
        sibh  = res.sibh;
        age   = res.age;
        y2bh  = res.y2bh;
        y2sh  = res.y2sh;

        diff = std::abs(si0 - si1);

        if (diff < tol) {
          converged = true;
          break;
        }
      }
      if (!converged){
        Rcpp::warning("Row " + std::to_string(i+1) + " did not converge. Difference was " + std::to_string(diff));
      }
    }
    else            {
      double      si0  = si[i];    // Solve from known SI total

      SIResult res;
      if      (s == "AW") res = AW_SI(sibh, sitype, ndx, th, si0, bh, sh);
      else if (s == "SB") res = SB_SI(sibh, sitype, ndx, th, si0, bh, sh);
      else if (s == "SW") res = SW_SI(sibh, sitype, ndx, th, si0, bh, sh);
      else                res = PL_SI(sibh, sitype, ndx, th, si0, bh, sh);

      si0   = res.si0;
      si1   = res.si1;
      sibh  = res.sibh;
      age   = res.age;
      y2bh  = res.y2bh;
      y2sh  = res.y2sh;
    }
    if      (ndx <  4)                  result[i] = age;
    else if (ndx == 4 && sitype == "B") result[i] = si1;
    else if (ndx == 4 && sitype == "T") result[i] = sibh;
    else if (ndx == 5)                  result[i] = y2bh;
    else                                result[i] = y2sh;
  }

  return result;
}
