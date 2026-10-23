#include <Rcpp.h>
#include <cmath>
#include <string>
using namespace Rcpp;

// ---------------------------------------------------------------------------
// Species-specific SI solvers
// Each returns a length-4 vector: {si0, si1, si_bh, Y2BH}
// ---------------------------------------------------------------------------

static inline NumericVector AW_SI(double age, double topht, double si0) {
  const double b1 =  9.908888;
  const double b2 = -3.92451;
  const double b3 = -0.32778;
  const double b4 =  0.134376;

  double k1   = exp(b1 + b2*sqrt(log(51.0)) + b3*pow(log(si0), 2) + b4*sqrt(50.0));
  double k2   = pow(si0, b3*log(si0));
  double k3   = (si0*(1.0 + k1)/1.3 - 1.0) / (exp(b1)*exp(b4*sqrt(50.0))*k2);
  double Y2BH = exp(pow(log(k3)/b2, 2)) - 1.0;

  double x10  = 1.0 + exp(b1 + b2*sqrt(log(age  + 1.0)) + b3*pow(log(si0), 2) + b4*sqrt(50.0));
  double x20  = 1.0 + exp(b1 + b2*sqrt(log(51.0))       + b3*pow(log(si0), 2) + b4*sqrt(50.0));
  double si1  = topht * x10 / x20;
  si0         = (si0 + si1) / 2.0;

  double si_bh = si1 * (1.0 + exp(b1 + b2*sqrt(log(51.0))          + b3*pow(log(si1), 2) + b4*sqrt(50.0))) /
                       (1.0 + exp(b1 + b2*sqrt(log(50.0+Y2BH+1.0)) + b3*pow(log(si1), 2) + b4*sqrt(50.0)));

  return NumericVector::create(si0, si1, si_bh, Y2BH);
}

static inline NumericVector SB_SI(double age, double topht, double si0) {
  const double b1 =  14.56236;
  const double b2 =  -6.04705;
  const double b3 =  -1.53715;
  const double b4 =   0.240174;

  double k1   = exp(b1 + b2*sqrt(log(51.0)) + b3*log(si0) + b4*sqrt(50.0));
  double k2   = pow(si0, b3);
  double k3   = (si0*(1.0 + k1)/1.3 - 1.0) / (exp(b1)*exp(b4*sqrt(50.0))*k2);
  double Y2BH = exp(pow(log(k3)/b2, 2)) - 1.0;

  double x10  = 1.0 + exp(b1 + b2*sqrt(log(age  + 1.0)) + b3*log(si0) + b4*sqrt(50.0));
  double x20  = 1.0 + exp(b1 + b2*sqrt(log(51.0))       + b3*log(si0) + b4*sqrt(50.0));
  double si1  = topht * x10 / x20;
  si0         = (si0 + si1) / 2.0;

  double si_bh = si1 * (1.0 + exp(b1 + b2*sqrt(log(51.0))          + b3*log(si1) + b4*sqrt(50.0))) /
                       (1.0 + exp(b1 + b2*sqrt(log(50.0+Y2BH+1.0)) + b3*log(si1) + b4*sqrt(50.0)));

  return NumericVector::create(si0, si1, si_bh, Y2BH);
}

static inline NumericVector SW_SI(double age, double topht, double si0) {
  const double b1 =  12.14943;
  const double b2 =  -3.77051;
  const double b3 =  -0.28534;
  const double b4 =   0.165483;

  double k1   = exp(b1 + b2*sqrt(log(pow(50.0,2)+1.0)) + b3*pow(log(si0),2) + b4*sqrt(50.0));
  double k2   = pow(si0, b3*log(si0));
  double k3   = (si0*(1.0 + k1)/1.3 - 1.0) / (exp(b1)*exp(b4*sqrt(50.0))*k2);
  double Y2BH = sqrt(exp(pow(log(k3)/b2, 2)) - 1.0);

  double x10  = 1.0 + exp(b1 + b2*sqrt(log(pow(age,2)+1.0))  + b3*pow(log(si0),2) + b4*sqrt(50.0));
  double x20  = 1.0 + exp(b1 + b2*sqrt(log(pow(50.0,2)+1.0)) + b3*pow(log(si0),2) + b4*sqrt(50.0));
  double si1  = topht * x10 / x20;
  si0         = (si0 + si1) / 2.0;

  double si_bh = si1 * (1.0 + exp(b1 + b2*sqrt(log(pow(50.0,2)+1.0))          + b3*pow(log(si1),2) + b4*sqrt(50.0))) /
                       (1.0 + exp(b1 + b2*sqrt(log(pow(50.0+Y2BH,2)+1.0))     + b3*pow(log(si1),2) + b4*sqrt(50.0)));

  return NumericVector::create(si0, si1, si_bh, Y2BH);
}

static inline NumericVector PL_SI(double age, double topht, double si0) {
  const double b1 =  12.84571;
  const double b2 =  -5.73936;
  const double b3 =  -0.91312;
  const double b4 =   0.150668;

  double k1   = exp(b1 + b2*sqrt(log(51.0)) + b3*log(si0) + b4*sqrt(50.0));
  double k2   = pow(si0, b3);
  double k3   = (si0*(1.0 + k1)/1.3 - 1.0) / (exp(b1)*exp(b4*sqrt(50.0))*k2);
  double Y2BH = exp(pow(log(k3)/b2, 2)) - 1.0;

  double x10  = 1.0 + exp(b1 + b2*sqrt(log(age   + 1.0)) + b3*log(si0) + b4*sqrt(50.0));
  double x20  = 1.0 + exp(b1 + b2*sqrt(log(51.0))        + b3*log(si0) + b4*sqrt(50.0));
  double si1  = topht * x10 / x20;
  si0         = (si0 + si1) / 2.0;

  double si_bh = si1 * (1.0 + exp(b1 + b2*sqrt(log(51.0))              + b3*log(si1) + b4*sqrt(50.0))) /
                       (1.0 + exp(b1 + b2*sqrt(log(50.0+Y2BH+1.0))     + b3*log(si1) + b4*sqrt(50.0)));

  return NumericVector::create(si0, si1, si_bh, Y2BH);
}

// ---------------------------------------------------------------------------
// Main exported function — mirrors GYPSY_SI_vec() interface
//
// @param spp    Character vector: "AW", "PL", "SW", or "SB"
// @param topht  Numeric vector: top height (must be > 0.1)
// @param age    Numeric vector: stand age
// @param ndx    Integer 1-3: which output column to return
//               1 = si0 (intermediate), 2 = si1 (SI at total age), 3 = si_bh (SI at BH age)
// @param iter   Max convergence iterations (default 1000)
// @param tol    Convergence tolerance (default 1e-8)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector GYPSY_SI_cpp(CharacterVector spp,
                           NumericVector   topht,
                           NumericVector   age,
                           int             ndx  = 2,
                           int             iter = 1000,
                           double          tol  = 1e-8) {

  int n = spp.size();

  // --- Input validation ---
  if (topht.size() != n || age.size() != n)
    stop("spp, topht, and age must all be the same length.");
  if (ndx < 1 || ndx > 3)
    stop("ndx must be 1, 2, or 3.");

  std::vector<std::string> errors;
  for (int i = 0; i < n; i++) {
    std::string s = as<std::string>(spp[i]);
    if (s != "AW" && s != "PL" && s != "SW" && s != "SB")
      errors.push_back("Row " + std::to_string(i+1) + ": species must be AW, PL, SW, or SB.");
    if (topht[i] < 0.1 || NumericVector::is_na(topht[i]))
      errors.push_back("Row " + std::to_string(i+1) + ": topht must be >= 0.1.");
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

  // --- Output vector ---
  NumericVector result(n);

  // --- Solve row by row ---
  for (int i = 0; i < n; i++) {
    std::string s = as<std::string>(spp[i]);
    double th  = topht[i];
    double ag  = age[i];
    double si0 = 10.0;
    double si1 = 100.0;
    double si_bh = 0.0, Y2BH = 0.0;

    // Convergence loop
    for (int it = 0; it < iter; it++) {
      NumericVector res;
      if      (s == "AW") res = AW_SI(ag, th, si0);
      else if (s == "SB") res = SB_SI(ag, th, si0);
      else if (s == "SW") res = SW_SI(ag, th, si0);
      else                res = PL_SI(ag, th, si0);

      si0   = res[0];
      si1   = res[1];
      si_bh = res[2];
      Y2BH  = res[3];

      if (std::abs(si0 - si1) < tol) break;
    }

    // ndx: 1=si0, 2=si1, 3=si_bh  (matches R version's [[ndx+1]])
    if      (ndx == 1) result[i] = si0;
    else if (ndx == 2) result[i] = si1;
    else               result[i] = si_bh;
  }

  return result;
}
