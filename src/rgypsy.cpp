#include <Rcpp.h>
#include "core/gypsy_components.h"
#include "core/gypsy_sim.h"
#include "core/gypsy_ti.h"
using namespace Rcpp;

// ---------------------------------------------------------------------------
// R bindings for the GYPSY core (src/core)
//
// The model itself is plain C++ in src/core. These functions only convert
// between R and C++ types: NA <-> NaN, core warnings -> R warnings, core
// errors (gypsy::Error) -> R errors.
// ---------------------------------------------------------------------------

namespace {

std::vector<double>      num(const NumericVector& x)   { return as<std::vector<double>>(x); }
std::vector<std::string> chr(const CharacterVector& x) { return as<std::vector<std::string>>(x); }

NumericVector to_r(const std::vector<double>& x) {
  NumericVector out(x.size());
  for (std::size_t i = 0; i < x.size(); i++) out[i] = std::isnan(x[i]) ? NA_REAL : x[i];
  return out;
}

void emit(const gypsy::Messages& warn) {
  for (const std::string& w : warn) Rcpp::warning(w);
}

// Named list of columns; `int_cols` are returned as integer
List to_r(const gypsy::Table& t, std::initializer_list<const char*> int_cols = {}) {
  List out(t.names.size());
  for (std::size_t k = 0; k < t.names.size(); k++) {
    bool is_int = false;
    for (const char* nm : int_cols) if (t.names[k] == nm) is_int = true;
    if (is_int) out[k] = IntegerVector(t.cols[k].begin(), t.cols[k].end());
    else        out[k] = to_r(t.cols[k]);
  }
  out.attr("names") = wrap(t.names);
  return out;
}

void interrupt() { Rcpp::checkUserInterrupt(); }

} // namespace

// ---- Site index and age -------------------------------------------------------

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
  gypsy::Messages warn;
  std::vector<double> res = gypsy::SI(chr(spp), num(ht), num(age), agetype, ndx, iter, tol, bh, sh, warn);
  emit(warn);
  return to_r(res);
}

// [[Rcpp::export]]
NumericVector GYPSY_Age_cpp(CharacterVector spp,
                            NumericVector   si,
                            Nullable<NumericVector> h = R_NilValue,
                            std::string     sitype = "B",
                            int             ndx  = 1,
                            int             iter = 1000,
                            double          tol  = 1e-8,
                            double          bh   = 1.3,
                            double          sh   = 0.3) {
  std::vector<double> ht = h.isNotNull() ? num(NumericVector(h))
                                         : std::vector<double>(spp.size(), gypsy::NaN);
  gypsy::Messages warn;
  std::vector<double> res = gypsy::Age(chr(spp), num(si), ht, sitype, ndx, iter, tol, bh, sh, warn);
  emit(warn);
  return to_r(res);
}

// [[Rcpp::export]]
NumericVector GYPSY_SI_TI_cpp(CharacterVector spp, NumericVector si_t, NumericVector gw,
                              NumericVector sel_age) {
  gypsy::Messages warn;
  std::vector<double> res = gypsy::SI_TI(chr(spp), num(si_t), num(gw), num(sel_age), warn);
  emit(warn);
  return to_r(res);
}

// ---- Stand-level components -----------------------------------------------------

// [[Rcpp::export]]
Rcpp::List GYPSY_PSI_cpp(CharacterVector stand,
                         NumericVector aw_ps, NumericVector aw_age_t, NumericVector aw_si_b,
                         NumericVector sb_ps, NumericVector sb_age_t, NumericVector sb_si_b,
                         NumericVector sw_ps, NumericVector sw_age_t, NumericVector sw_si_b,
                         NumericVector pl_ps, NumericVector pl_age_t, NumericVector pl_si_b) {
  gypsy::Cols4 ps    = {num(aw_ps),    num(sb_ps),    num(sw_ps),    num(pl_ps)};
  gypsy::Cols4 age_t = {num(aw_age_t), num(sb_age_t), num(sw_age_t), num(pl_age_t)};
  gypsy::Cols4 si_b  = {num(aw_si_b),  num(sb_si_b),  num(sw_si_b),  num(pl_si_b)};
  gypsy::Cols4 psi;
  gypsy::PSI(ps, age_t, si_b, psi);
  return List::create(_["aw_psi"] = to_r(psi[0]), _["sb_psi"] = to_r(psi[1]),
                      _["sw_psi"] = to_r(psi[2]), _["pl_psi"] = to_r(psi[3]));
}

// [[Rcpp::export]]
Rcpp::List GYPSY_SDF_cpp(CharacterVector stand,
                         NumericVector aw_age_b, NumericVector aw_n, NumericVector aw_si_b, NumericVector aw_psi,
                         NumericVector sb_age_t, NumericVector sb_n, NumericVector sb_si_b, NumericVector sb_psi,
                         NumericVector sw_age_t, NumericVector sw_n, NumericVector sw_si_b, NumericVector sw_psi,
                         NumericVector pl_age_t, NumericVector pl_n, NumericVector pl_si_b, NumericVector pl_psi,
                         int iter = 1000, double tol = 1e-8) {
  gypsy::Cols4 age  = {num(aw_age_b), num(sb_age_t), num(sw_age_t), num(pl_age_t)};
  gypsy::Cols4 den  = {num(aw_n),     num(sb_n),     num(sw_n),     num(pl_n)};
  gypsy::Cols4 si_b = {num(aw_si_b),  num(sb_si_b),  num(sw_si_b),  num(pl_si_b)};
  gypsy::Cols4 psi  = {num(aw_psi),   num(sb_psi),   num(sw_psi),   num(pl_psi)};
  gypsy::Cols4 sdf, den0;
  gypsy::Messages warn;
  gypsy::SDF(age, den, si_b, psi, iter, tol, sdf, den0, warn);
  emit(warn);
  return List::create(
    _["aw_sdf"]  = to_r(sdf[0]),  _["sb_sdf"]  = to_r(sdf[1]),
    _["sw_sdf"]  = to_r(sdf[2]),  _["pl_sdf"]  = to_r(sdf[3]),
    _["aw_den0"] = to_r(den0[0]), _["sb_den0"] = to_r(den0[1]),
    _["sw_den0"] = to_r(den0[2]), _["pl_den0"] = to_r(den0[3]));
}

// ---- Species-level components ---------------------------------------------------

// [[Rcpp::export]]
NumericVector GYPSY_Ht_cpp(CharacterVector spp, NumericVector si_t, NumericVector age_t) {
  return to_r(gypsy::Ht(chr(spp), num(si_t), num(age_t)));
}

// [[Rcpp::export]]
NumericVector GYPSY_N_cpp(CharacterVector spp, NumericVector age, NumericVector sdf,
                          NumericVector si_b, NumericVector psi,
                          NumericVector sdf_aw, NumericVector sdf_sw, NumericVector sdf_sb) {
  return to_r(gypsy::N(chr(spp), num(age), num(sdf), num(si_b), num(psi),
                       num(sdf_aw), num(sdf_sw), num(sdf_sb)));
}

// [[Rcpp::export]]
NumericVector GYPSY_BAinc_cpp(CharacterVector spp, NumericVector bhage, NumericVector sc,
                              NumericVector den0, NumericVector si_b, NumericVector psi,
                              NumericVector ba, NumericVector sdf_aw, NumericVector sdf_sb,
                              NumericVector sdf_sw, NumericVector sdf_pl) {
  return to_r(gypsy::BAinc(chr(spp), num(bhage), num(sc), num(den0), num(si_b), num(psi),
                           num(ba), num(sdf_aw), num(sdf_sb), num(sdf_sw), num(sdf_pl)));
}

// [[Rcpp::export]]
Rcpp::List GYPSY_Vol_cpp(CharacterVector spp, NumericVector ba, NumericVector topht,
                         NumericVector n, NumericVector bhage, NumericVector stump_ht,
                         NumericVector stump_dob, NumericVector top_dib) {
  gypsy::VolResult v = gypsy::Vol(chr(spp), num(ba), num(topht), num(n), num(bhage),
                                  num(stump_ht), num(stump_dob), num(top_dib));
  return List::create(_["tvol"] = to_r(v.tvol), _["mvol"] = to_r(v.mvol),
                      _["mn"]   = to_r(v.mn),   _["qmd"]  = to_r(v.qmd));
}

// ---- Full simulation ---------------------------------------------------------------
// stand: 0-based stand code per row; stand_label: name of each code.
// stump_ht, stump_dob, top_dib: per species AW, SB, SW, PL.
// ti: run SimTI with gw and sel_age (ignored otherwise).
// Returns rows (one per stand and year), init (one per stand) and warnings.

// [[Rcpp::export]]
Rcpp::List GYPSY_Sim_cpp(IntegerVector   stand,
                         CharacterVector stand_label,
                         NumericVector   stand_age,
                         CharacterVector spp,
                         NumericVector   n,
                         NumericVector   si,
                         NumericVector   age,
                         NumericVector   ht,
                         NumericVector   ps,
                         NumericVector   ba,
                         std::string     sitype,
                         std::string     agetype,
                         int             simlen,
                         bool            psflag,
                         bool            baflag,
                         double          bh,
                         double          sh,
                         NumericVector   stump_ht,
                         NumericVector   stump_dob,
                         NumericVector   top_dib,
                         int             iter,
                         double          tol,
                         bool            detail,
                         bool            ti      = false,
                         NumericVector   gw      = NumericVector(),
                         NumericVector   sel_age = NumericVector()) {
  gypsy::SimInput in;
  in.stand       = as<std::vector<int>>(stand);
  in.stand_label = chr(stand_label);
  in.stand_age   = num(stand_age);
  in.spp         = chr(spp);
  in.n  = num(n);  in.si = num(si); in.age = num(age);
  in.ht = num(ht); in.ps = num(ps); in.ba  = num(ba);

  gypsy::SimOptions opt;
  opt.sitype = sitype; opt.agetype = agetype; opt.simlen = simlen;
  opt.psflag = psflag; opt.baflag = baflag; opt.bh = bh; opt.sh = sh;
  for (int s = 0; s < 4; s++) {
    opt.merch.stump_ht[s]  = stump_ht[s];
    opt.merch.stump_dob[s] = stump_dob[s];
    opt.merch.top_dib[s]   = top_dib[s];
  }
  opt.iter = iter; opt.tol = tol; opt.detail = detail;
  opt.check_interrupt = interrupt;

  gypsy::SimResult res;
  if (ti) {
    gypsy::TIInput t;
    t.gw = num(gw); t.sel_age = num(sel_age);
    res = gypsy::SimTI(in, t, opt);
  } else {
    res = gypsy::Sim(in, opt);
  }
  return List::create(_["rows"]     = to_r(res.rows, {"stand_idx", "year"}),
                      _["init"]     = to_r(res.init),
                      _["warnings"] = wrap(res.warnings));
}
