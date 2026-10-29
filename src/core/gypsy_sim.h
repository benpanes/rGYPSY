#ifndef GYPSY_SIM_H
#define GYPSY_SIM_H

#include <map>
#include "gypsy_eq.h"
#include "gypsy_si.h"
#include "gypsy_stand.h"
#include "gypsy_project.h"
#include "gypsy_ti.h"

// ---------------------------------------------------------------------------
// Full GYPSY run (GYPSY_run.sas): one row per stand and species in, one row
// per stand and year out.
//
//   1. Check inputs
//   2. SI and age per row (SAS data v1): SI from height and age, a supplied
//      SI overrides it, and total age is recomputed from SI and height when
//      both are given
//   3. Leading (maximum) total age per stand, species widened to columns
//   4. Percent stocking index, stand density factor and initial density
//   5. Yearly projection (gypsy_project.h)
//
// SimTI runs the same steps with each row's SI_t raised for tree
// improvement (gypsy_ti.h) after step 2. Total age is kept; SI_bh, years to
// breast height and breast height age follow the new SI_t.
//
// Stands are identified by a code 0 .. n_stands - 1 per row; the output is
// in code order. stand_label gives each code's name for messages. Missing
// values are NaN.
// ---------------------------------------------------------------------------

namespace gypsy {

struct SimInput {
  std::vector<int>         stand;         // stand code per row
  std::vector<std::string> stand_label;   // name of each stand code
  std::vector<double>      stand_age;
  std::vector<std::string> spp;           // "AW", "SB", "SW" or "PL"
  std::vector<double>      n, si, age, ht, ps, ba;
};

struct SimOptions {
  std::string sitype  = "B";   // SI type: B breast height age 50, T total age 50
  std::string agetype = "T";   // age type: T total, B breast height, S stump height
  int    simlen = 200;
  bool   psflag = true;        // false ignores ps
  bool   baflag = true;        // false ignores ba
  double bh = 1.3, sh = 0.3;   // breast and stump height (m) for SI and age
  MerchStd merch;
  int    iter = 1000;
  double tol  = 1e-8;
  bool   detail = false;       // also return sc_b_*, ba1_*, ba2_*
  void (*check_interrupt)() = nullptr;
};

// Tree improvement inputs for SimTI, one value per row (gw 0: no
// adjustment; NaN is treated as 0 for SW and PL, with a warning)
struct TIInput {
  std::vector<double> gw;        // genetic worth for height at the selection age (%)
  std::vector<double> sel_age;   // selection age
};

struct SimResult {
  Table    rows;       // one row per stand and year; stand_idx = code + 1
  Table    init;       // one row per stand (stand_age, totage_max, widened inputs, PSI, SDF, ...)
  Messages warnings;
};

// Sim and SimTI; ti is null for Sim
inline SimResult SimRun(const SimInput& in, const SimOptions& opt, const TIInput* ti) {
  const char* SPP[4] = {"aw", "sb", "sw", "pl"};
  std::size_t nr = in.spp.size();
  std::size_t ns = in.stand_label.size();
  SimResult out;
  Messages& warn = out.warnings;

  // ---- 1. Input checks ----------------------------------------------------------
  for (std::size_t m : {in.stand.size(), in.stand_age.size(), in.n.size(), in.si.size(),
                        in.age.size(), in.ht.size(), in.ps.size(), in.ba.size()})
    if (m != nr) throw Error("All inputs must have one value per row");
  for (int c : in.stand)
    if (c < 0 || static_cast<std::size_t>(c) >= ns) throw Error("Stand codes must be 0 to the number of stands - 1");

  if (opt.sitype != "B" && opt.sitype != "T") throw Error("sitype must be B or T");
  if (opt.agetype != "T" && opt.agetype != "B" && opt.agetype != "S") throw Error("agetype must be T, B, or S");
  char at = opt.agetype[0];

  std::vector<int> spp(nr), bad_spp, miss_input;
  for (std::size_t i = 0; i < nr; i++) {
    spp[i] = spp_index(in.spp[i]);
    if (spp[i] < 0) bad_spp.push_back(static_cast<int>(i));
  }
  if (!bad_spp.empty())
    throw Error("Species must be AW, SB, SW, or PL. Check row(s): " + join_rows(bad_spp));

  for (std::size_t i = 0; i < nr; i++)   // Check for missing input data
    if (miss(in.si[i]) + miss(in.age[i]) + miss(in.ht[i]) > 1) miss_input.push_back(static_cast<int>(i));
  if (!miss_input.empty())
    throw Error("Simulation requires 2 of SI, Age, and Top Height. Check row(s): " + join_rows(miss_input));

  std::vector<double> gw_in;   // SimTI: genetic worth, missing SW and PL values as 0
  if (ti) {
    for (std::size_t m : {ti->gw.size(), ti->sel_age.size()})
      if (m != nr) throw Error("All inputs must have one value per row");
    gw_in = check_ti(spp, ti->gw, ti->sel_age, warn);
  }

  for (int s = 0; s < 4; s++)
    if (opt.merch.stump_ht[s] > 0.3) throw Error("stump_ht must be <= 0.3");

  // Duplicate stand species groups, listed group by group in order of first
  // appearance
  {
    std::map<int, std::vector<int>> groups;
    std::vector<int> order;
    for (std::size_t i = 0; i < nr; i++) {
      int key = in.stand[i]*4 + spp[i];
      if (groups.find(key) == groups.end()) order.push_back(key);
      groups[key].push_back(static_cast<int>(i));
    }
    std::vector<int> dup;
    for (int key : order)
      if (groups[key].size() > 1) dup.insert(dup.end(), groups[key].begin(), groups[key].end());
    if (!dup.empty())
      throw Error("Duplicate stand species groups in row(s): " + join_rows(dup));
  }

  std::vector<double> ps = in.ps, ba = in.ba;
  if (!opt.psflag) ps.assign(nr, NaN);
  if (!opt.baflag) ba.assign(nr, NaN);

  // ---- 2. SI and age (SAS data v1) ----------------------------------------------
  std::vector<double> si_t(nr, NaN), si_bh(nr, NaN), y2bh(nr, NaN), totage(nr, NaN), bhage(nr, NaN);
  auto age_from = [&](std::size_t i, double si, double ht, char st) {
    AgeResult r = SolveAge(spp[i], si, ht, st, opt.iter, opt.tol, opt.bh, opt.sh);
    if (!r.converged) warn.push_back("Row " + std::to_string(i + 1) + ": age did not converge.");
    return r;
  };

  // SI from height and age
  {
    std::vector<std::size_t> rows;
    std::vector<std::string> errors;
    for (std::size_t i = 0; i < nr; i++) {
      double min_age = (spp[i] == AW || spp[i] == PL) ? 7.5 : 9.5;
      if (miss(in.age[i]) || miss(in.ht[i]) || !(in.ht[i] > 0) || !(in.age[i] >= min_age)) continue;
      rows.push_back(i);
      std::string row = "Row " + std::to_string(i + 1);
      if (at == 'T' && in.ht[i] < 0.1)      errors.push_back(row + ": ht must be >= 0.1.");
      if (at == 'S' && in.ht[i] < opt.sh)   errors.push_back(row + ": ht must be >= stump height.");
      if (at == 'B' && in.ht[i] < opt.bh)   errors.push_back(row + ": ht must be >= breast height.");
    }
    if (!errors.empty()) {
      std::string msg;
      for (const std::string& e : errors) msg += e + "\n";
      throw Error(msg);
    }
    for (std::size_t i : rows) {
      SIResult r = SolveSI(spp[i], in.ht[i], in.age[i], at, opt.iter, opt.tol, opt.bh, opt.sh);
      if (!r.converged) warn.push_back("Row " + std::to_string(i + 1) + ": SI did not converge.");
      si_t[i] = r.si_t; si_bh[i] = r.si_bh; y2bh[i] = r.y2bh;
    }
  }

  // A supplied SI overrides it
  for (std::size_t i = 0; i < nr; i++) {
    if (!(sv(in.si[i]) > 0)) continue;
    if (opt.sitype == "B") {
      si_bh[i] = in.si[i];
      si_t[i]  = age_from(i, in.si[i], NaN, 'B').si_t;
    } else {
      si_t[i]  = in.si[i];
      si_bh[i] = age_from(i, in.si[i], NaN, 'T').si_bh;
    }
    y2bh[i] = age_from(i, si_t[i], NaN, 'T').y2bh;
  }

  // Total age from the supplied age
  for (std::size_t i = 0; i < nr; i++) {
    if      (at == 'T') totage[i] = in.age[i];
    else if (at == 'B') totage[i] = in.age[i] + y2bh[i];
    else if (!miss(si_t[i])) totage[i] = in.age[i] + age_from(i, si_t[i], NaN, 'T').y2sh;
  }

  // Total age from SI and height when both are given
  for (std::size_t i = 0; i < nr; i++)
    if (sv(in.si[i]) > 0 && sv(in.ht[i]) > 0)
      totage[i] = age_from(i, si_t[i], in.ht[i], 'T').tage;

  // ---- 2b. Tree improvement (SimTI) --------------------------------------------------
  std::vector<double> si_t_base = si_t, si_bh_base = si_bh, aacor(nr, NaN);
  if (ti) {
    for (std::size_t i = 0; i < nr; i++) {
      double gw = gw_in[i];
      if (miss(si_t[i]) || miss(gw) || gw == 0) continue;
      aacor[i] = AgeAgeCor(spp[i], ti->sel_age[i]);
      si_t[i]  = SI_TI(spp[i], si_t[i], gw, ti->sel_age[i]);
      AgeResult r = age_from(i, si_t[i], NaN, 'T');
      si_bh[i] = r.si_bh; y2bh[i] = r.y2bh;
    }
  }

  std::vector<int> no_si;
  for (std::size_t i = 0; i < nr; i++) {
    if (sv(y2bh[i]) > 0) bhage[i] = totage[i] - y2bh[i];
    if (miss(si_t[i])) no_si.push_back(static_cast<int>(i));
  }
  if (!no_si.empty())
    warn.push_back("No SI for row(s) " + join_rows(no_si) +
                   " (age below 7.5 for AW/PL or 9.5 for SB/SW); these species are not grown but still count toward the leading age.");

  // ---- 3. Leading age and species columns ---------------------------------------
  std::vector<StandInit> st(ns);
  std::vector<bool> seen(ns, false);
  std::vector<double> wide_n[4], wide_ps[4], wide_ht[4];
  std::vector<double> wide_sit0[4], wide_sib0[4], wide_gw[4], wide_sel[4], wide_r[4];   // SimTI
  for (int s = 0; s < 4; s++) { wide_n[s].assign(ns, NaN); wide_ps[s].assign(ns, NaN); wide_ht[s].assign(ns, NaN); }
  if (ti)
    for (int s = 0; s < 4; s++)
      for (std::vector<double>* w : {&wide_sit0[s], &wide_sib0[s], &wide_gw[s], &wide_sel[s], &wide_r[s]})
        w->assign(ns, NaN);
  for (StandInit& x : st) {
    x.stand_age = NaN; x.totage_max = NaN;
    for (int s = 0; s < 4; s++) {
      x.si_t[s] = x.si_bh[s] = x.totage[s] = x.bhage[s] = x.y2bh[s] = x.ba[s] = NaN;
      x.psi[s] = x.sdf[s] = x.den0[s] = NaN;
    }
  }
  for (std::size_t i = 0; i < nr; i++) {
    int c = in.stand[i], s = spp[i];
    StandInit& x = st[c];
    if (!seen[c]) { x.stand_age = in.stand_age[i]; seen[c] = true; }
    if (!miss(totage[i]) && (miss(x.totage_max) || totage[i] > x.totage_max)) x.totage_max = totage[i];
    x.si_t[s] = si_t[i]; x.si_bh[s] = si_bh[i]; x.totage[s] = totage[i];
    x.bhage[s] = bhage[i]; x.y2bh[s] = y2bh[i]; x.ba[s] = ba[i];
    wide_n[s][c] = in.n[i]; wide_ps[s][c] = ps[i]; wide_ht[s][c] = in.ht[i];
    if (ti) {
      wide_sit0[s][c] = si_t_base[i]; wide_sib0[s][c] = si_bh_base[i];
      wide_gw[s][c]   = gw_in[i];    wide_sel[s][c]  = ti->sel_age[i]; wide_r[s][c] = aacor[i];
    }
  }
  for (std::size_t c = 0; c < ns; c++)
    if (!seen[c]) throw Error("Stand " + in.stand_label[c] + " has no rows");

  // As in the SAS, these stands are kept and project as all NA
  {
    std::vector<bool> listed(ns, false);
    std::string no_age;
    for (std::size_t i = 0; i < nr; i++) {
      int c = in.stand[i];
      if (miss(st[c].totage_max) && !listed[c]) {
        no_age += (no_age.empty() ? "" : ", ") + in.stand_label[c];
        listed[c] = true;
      }
    }
    if (!no_age.empty()) warn.push_back("No total age could be determined for stand(s): " + no_age);
  }

  // ---- 4. PSI, SDF and initial density --------------------------------------------
  for (std::size_t c = 0; c < ns; c++) {
    StandInit& x = st[c];
    double p[4], n[4], sdf_age[4];
    for (int s = 0; s < 4; s++) { p[s] = wide_ps[s][c]; n[s] = wide_n[s][c]; }
    StandPSI(p, x.totage, x.si_bh, x.psi);

    for (int s = 0; s < 4; s++) sdf_age[s] = (s == AW) ? x.bhage[AW] : x.totage[s];   // AW density runs on BH age
    bool conv[4];
    StandSDF(sdf_age, n, x.si_bh, x.psi, opt.iter, opt.tol, x.sdf, x.den0, conv);
    for (int s = 0; s < 4; s++)
      if (!conv[s]) warn.push_back("Stand " + in.stand_label[c] + ": SDF did not converge.");
  }

  // ---- 5. Grow forward ---------------------------------------------------------------
  ProjectResult res = Project(st, opt.merch, opt.simlen, 180, opt.detail, opt.check_interrupt);
  out.rows = std::move(res.rows);

  // Stand-level table: inputs widened by species, then PSI, SDF and the
  // projection's starting values
  Table& init = out.init;
  std::vector<double>& t_age  = init.add("stand_age",  ns);
  std::vector<double>& t_tmax = init.add("totage_max", ns);
  for (std::size_t c = 0; c < ns; c++) { t_age[c] = st[c].stand_age; t_tmax[c] = st[c].totage_max; }
  for (int s = 0; s < 4; s++) {
    std::string x = SPP[s];
    std::vector<double>& c_sit = init.add(x + "_si_t",   ns);
    std::vector<double>& c_sib = init.add(x + "_si_bh",  ns);
    std::vector<double>& c_ta  = init.add(x + "_totage", ns);
    std::vector<double>& c_bha = init.add(x + "_bhage",  ns);
    std::vector<double>& c_y2  = init.add(x + "_y2bh",   ns);
    std::vector<double>& c_n   = init.add(x + "_n",      ns);
    std::vector<double>& c_ps  = init.add(x + "_ps",     ns);
    std::vector<double>& c_ba  = init.add(x + "_ba",     ns);
    std::vector<double>& c_ht  = init.add(x + "_ht",     ns);
    for (std::size_t c = 0; c < ns; c++) {
      c_sit[c] = st[c].si_t[s];  c_sib[c] = st[c].si_bh[s]; c_ta[c] = st[c].totage[s];
      c_bha[c] = st[c].bhage[s]; c_y2[c]  = st[c].y2bh[s];  c_n[c]  = wide_n[s][c];
      c_ps[c]  = wide_ps[s][c];  c_ba[c]  = st[c].ba[s];    c_ht[c] = wide_ht[s][c];
    }
  }
  if (ti) {   // SI before the adjustment, and the adjustment's inputs
    for (int s = 0; s < 4; s++) {
      std::string x = SPP[s];
      init.add(x + "_si_t_base",  ns) = wide_sit0[s];
      init.add(x + "_si_bh_base", ns) = wide_sib0[s];
      init.add(x + "_gw",         ns) = wide_gw[s];
      init.add(x + "_sel_age",    ns) = wide_sel[s];
      init.add(x + "_aacor",      ns) = wide_r[s];
    }
  }
  for (int s = 0; s < 4; s++) {
    std::vector<double>& col = init.add(std::string(SPP[s]) + "_psi", ns);
    for (std::size_t c = 0; c < ns; c++) col[c] = st[c].psi[s];
  }
  for (int s = 0; s < 4; s++) {
    std::vector<double>& col = init.add(std::string(SPP[s]) + "_sdf", ns);
    for (std::size_t c = 0; c < ns; c++) col[c] = st[c].sdf[s];
  }
  for (int s = 0; s < 4; s++) {
    std::vector<double>& col = init.add(std::string(SPP[s]) + "_den0", ns);
    for (std::size_t c = 0; c < ns; c++) col[c] = st[c].den0[s];
  }
  for (std::size_t k = 0; k < res.init.names.size(); k++)
    init.add(res.init.names[k], ns) = std::move(res.init.cols[k]);

  return out;
}

inline SimResult Sim(const SimInput& in, const SimOptions& opt) { return SimRun(in, opt, nullptr); }

inline SimResult SimTI(const SimInput& in, const TIInput& ti, const SimOptions& opt) { return SimRun(in, opt, &ti); }

} // namespace gypsy

#endif
