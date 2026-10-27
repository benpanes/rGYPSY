#ifndef GYPSY_COMPONENTS_H
#define GYPSY_COMPONENTS_H

#include "gypsy_eq.h"
#include "gypsy_si.h"
#include "gypsy_stand.h"

// ---------------------------------------------------------------------------
// Vector versions of the GYPSY components, one value per row
//
// Inputs must already be the same length (recycling is the caller's job).
// Invalid inputs throw gypsy::Error; rows that don't converge add a message
// to `warn`. Missing values are NaN in and out.
// ---------------------------------------------------------------------------

namespace gypsy {

using Vec    = std::vector<double>;
using StrVec = std::vector<std::string>;
using Cols4  = Vec[4];   // one column per species, AW, SB, SW, PL

inline void check_same_length(std::size_t n, std::initializer_list<std::size_t> sizes, const char* msg) {
  for (std::size_t m : sizes) if (m != n) throw Error(msg);
}

inline void throw_row_errors(const StrVec& errors) {
  if (errors.empty()) return;
  std::string msg;
  for (const std::string& e : errors) msg += e + "\n";
  throw Error(msg);
}

// ---- Site index from height and age ----------------------------------------
// ndx: 1 SI_t, 2 SI_bh, 3 years to breast height, 4 years to stump height
inline Vec SI(const StrVec& spp, const Vec& ht, const Vec& age, const std::string& agetype,
              int ndx, int iter, double tol, double bh, double sh, Messages& warn) {
  std::size_t n = spp.size();
  check_same_length(n, {ht.size(), age.size()}, "spp, ht, and age must all be the same length.");
  if (ndx < 1 || ndx > 4) throw Error("ndx must be 1, 2, 3, or 4.");
  if (agetype != "T" && agetype != "B" && agetype != "S") throw Error("agetype must be T, B, or S");

  StrVec errors;
  for (std::size_t i = 0; i < n; i++) {
    std::string row = "Row " + std::to_string(i + 1);
    int s = spp_index(spp[i]);
    if (s < 0)
      errors.push_back(row + ": species must be AW, PL, SW, or SB.");
    if (agetype == "T" && (ht[i] < 0.1 || miss(ht[i])))
      errors.push_back(row + ": ht must be >= 0.1.");
    if (agetype == "S" && (ht[i] < sh || miss(ht[i])))
      errors.push_back(row + ": ht must be >= stump height.");
    if (agetype == "B" && (ht[i] < bh || miss(ht[i])))
      errors.push_back(row + ": ht must be >= breast height.");
    if ((s == AW || s == PL) && age[i] < 7.5)
      errors.push_back(row + ": age must be >= 7.5 for AW or PL.");
    if ((s == SW || s == SB) && age[i] < 9.5)
      errors.push_back(row + ": age must be >= 9.5 for SW or SB.");
  }
  throw_row_errors(errors);

  Vec out(n);
  for (std::size_t i = 0; i < n; i++) {
    SIResult r = SolveSI(spp_index(spp[i]), ht[i], age[i], agetype[0], iter, tol, bh, sh);
    if (!r.converged)
      warn.push_back("Row " + std::to_string(i + 1) + " did not converge. Final diff was " + std::to_string(r.diff));
    if      (ndx == 1) out[i] = r.si_t;
    else if (ndx == 2) out[i] = r.si_bh;
    else if (ndx == 3) out[i] = r.y2bh;
    else               out[i] = r.y2sh;
  }
  return out;
}

// ---- Age from SI and height, SI conversion ----------------------------------
// ndx: 1 total age, 2 stump height age, 3 breast height age,
//      4 the other SI type, 5 years to breast height, 6 years to stump height
// ht can be all NaN for ndx 4-6.
inline Vec Age(const StrVec& spp, const Vec& si, const Vec& ht, const std::string& sitype,
               int ndx, int iter, double tol, double bh, double sh, Messages& warn) {
  std::size_t n = spp.size();
  check_same_length(n, {ht.size(), si.size()}, "input vectors must all be the same length.");
  if (ndx < 1 || ndx > 6) throw Error("ndx must be between 1-6.");
  if (sitype != "B" && sitype != "T") throw Error("site index type must be B, or T");

  StrVec errors;
  for (std::size_t i = 0; i < n; i++)
    if (spp_index(spp[i]) < 0)
      errors.push_back("Row " + std::to_string(i + 1) + ": species must be AW, PL, SW, or SB.");
  throw_row_errors(errors);

  Vec out(n);
  for (std::size_t i = 0; i < n; i++) {
    AgeResult r = SolveAge(spp_index(spp[i]), si[i], ht[i], sitype[0], iter, tol, bh, sh);
    if (!r.converged)
      warn.push_back("Row " + std::to_string(i + 1) + " did not converge. Difference was " + std::to_string(r.diff));
    if      (ndx == 1) out[i] = r.tage;
    else if (ndx == 2) out[i] = r.tage - r.y2sh;
    else if (ndx == 3) out[i] = r.tage - r.y2bh;
    else if (ndx == 4) out[i] = (sitype == "B") ? r.si_t : r.si_bh;
    else if (ndx == 5) out[i] = r.y2bh;
    else               out[i] = r.y2sh;
  }
  return out;
}

// ---- Percent stocking index, one row per stand -------------------------------
inline void PSI(const Cols4& ps, const Cols4& age_t, const Cols4& si_b, Cols4& psi) {
  std::size_t n = ps[0].size();
  for (int s = 0; s < 4; s++) psi[s].assign(n, NaN);
  for (std::size_t i = 0; i < n; i++) {
    double p[4], a[4], si[4], out[4];
    for (int s = 0; s < 4; s++) { p[s] = ps[s][i]; a[s] = age_t[s][i]; si[s] = si_b[s][i]; }
    StandPSI(p, a, si, out);
    for (int s = 0; s < 4; s++) psi[s][i] = out[s];
  }
}

// ---- Stand density factor and initial density, one row per stand -----------
// age: AW breast height age, total age for the others
inline void SDF(const Cols4& age, const Cols4& den, const Cols4& si_b, const Cols4& psi,
                int iter, double tol, Cols4& sdf, Cols4& den0, Messages& warn) {
  std::size_t n = age[0].size();
  for (int s = 0; s < 4; s++) { sdf[s].assign(n, NaN); den0[s].assign(n, NaN); }
  for (std::size_t i = 0; i < n; i++) {
    double a[4], d[4], si[4], p[4], f[4], d0[4];
    bool conv[4];
    for (int s = 0; s < 4; s++) { a[s] = age[s][i]; d[s] = den[s][i]; si[s] = si_b[s][i]; p[s] = psi[s][i]; }
    StandSDF(a, d, si, p, iter, tol, f, d0, conv);
    for (int s = 0; s < 4; s++) {
      if (!conv[s]) warn.push_back("Row " + std::to_string(i + 1) + ": SDF did not converge.");
      sdf[s][i] = f[s]; den0[s][i] = d0[s];
    }
  }
}

// ---- Top height from SI_t and total age ---------------------------------------
inline Vec Ht(const StrVec& spp, const Vec& si_t, const Vec& age_t) {
  Vec out(spp.size(), NaN);
  for (std::size_t i = 0; i < spp.size(); i++) {
    int s = spp_index(spp[i]);
    if (s >= 0) out[i] = Ht(s, si_t[i], age_t[i]);
  }
  return out;
}

// ---- Density at age from SDF ---------------------------------------------------
// sdf_aw, sdf_sw, sdf_sb: the other species' SDFs (SW uses AW; PL uses AW,
// SW and SB)
inline Vec N(const StrVec& spp, const Vec& age, const Vec& sdf, const Vec& si_b, const Vec& psi,
             const Vec& sdf_aw, const Vec& sdf_sw, const Vec& sdf_sb) {
  Vec out(spp.size(), NaN);
  for (std::size_t i = 0; i < spp.size(); i++) {
    int s = spp_index(spp[i]);
    if (s >= 0 && sv(sdf[i]) > 0) {
      double comp[4] = {sdf_aw[i], sdf_sb[i], sdf_sw[i], NaN};
      out[i] = Density(s, age[i], sdf[i], si_b[i], psi[i], comp);
    }
  }
  return out;
}

// ---- Annual basal area increment -------------------------------------------------
inline Vec BAinc(const StrVec& spp, const Vec& bhage, const Vec& sc, const Vec& den0, const Vec& si_b,
                 const Vec& psi, const Vec& ba, const Vec& sdf_aw, const Vec& sdf_sb,
                 const Vec& sdf_sw, const Vec& sdf_pl) {
  Vec out(spp.size(), NaN);
  for (std::size_t i = 0; i < spp.size(); i++) {
    int s = spp_index(spp[i]);
    if (s >= 0) {
      double comp[4] = {nz(sdf_aw[i]), nz(sdf_sb[i]), nz(sdf_sw[i]), nz(sdf_pl[i])};
      out[i] = BAinc(s, bhage[i], sc[i], den0[i], si_b[i], psi[i], ba[i], comp);
    }
  }
  return out;
}

// ---- Total and merchantable volume (SAS examp18) -----------------------------------
// Total volume needs top height > 1.3 m; merchantable values also need n > 0
// (through QMD). Merchantable volume is rounded to 1e-5 as in the SAS.
struct VolResult { Vec tvol, mvol, mn, qmd; };

inline VolResult Vol(const StrVec& spp, const Vec& ba, const Vec& topht, const Vec& n, const Vec& bhage,
                     const Vec& stump_ht, const Vec& stump_dob, const Vec& top_dib) {
  std::size_t len = spp.size();
  VolResult out{Vec(len, NaN), Vec(len, NaN), Vec(len, NaN), Vec(len, NaN)};
  for (std::size_t i = 0; i < len; i++) {
    int s = spp_index(spp[i]);
    if (s < 0 || miss(ba[i])) continue;

    double k = (sv(n[i]) > 0) ? std::sqrt(ba[i]/(n[i]/10000.0)) : NaN;
    out.qmd[i] = k;

    if (sv(topht[i]) > 1.3) {
      double tv = TotVol(s, ba[i], topht[i]);
      out.tvol[i] = tv;
      if (sv(k) > 0) {
        MerchResult m = Merch(s, tv, k, topht[i], bhage[i], n[i], stump_ht[i], stump_dob[i], top_dib[i]);
        out.mvol[i] = sas_round(m.mv, 0.00001);
        out.mn[i]   = m.mn;
      }
    }
  }
  return out;
}

} // namespace gypsy

#endif
