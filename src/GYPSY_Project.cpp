#include <Rcpp.h>
#include <cmath>
#include <string>
#include <vector>
#include "GYPSY_eq.h"
using namespace Rcpp;

// ---------------------------------------------------------------------------
// Yearly GYPSY projection (SAS examp6 - examp19)
//
// Input `st` is one row per stand (sorted by stand) with columns
//   stand_age, totage_max, and for each species prefix aw_, sb_, sw_, pl_:
//   si_t, si_bh, totage, bhage, y2bh, ba, psi, sdf, den0
// (missing species are NA). Rows are produced for `simlen` years starting
// at the fractional part of the leading species' total age.
//
// The SAS data step keeps a variable's value from the previous year when its
// condition isn't met, and treats missing as below every number; both are
// reproduced here (see sv() in GYPSY_eq.h).
//
// Steps:
//   1. Species composition at BH age 0 of each species  (examp7 - examp11)
//   2. BA at the fractional BH age                      (examp13)
//   3. Yearly top height, density, composition, BA      (examp14)
//   4. BA localized to observed BA                      (examp14_1 - examp17)
//   5. Volume, merchantable volume and density, MAI     (examp18 - examp19)
// ---------------------------------------------------------------------------

namespace {

using namespace gypsy;

const char* PRE[4] = {"aw_", "sb_", "sw_", "pl_"};

// Upper limits on observed/predicted BA ratio for BA localization (examp17)
const double RATE_CAP[4] = {10.0, 5.0, 100.0, 10.0};

inline double na(double x) { return std::isnan(x) ? NA_REAL : x; }

} // namespace

// [[Rcpp::export]]
Rcpp::List GYPSY_Project_cpp(List          st,
                             NumericVector stump_ht,     // per species AW, SB, SW, PL
                             NumericVector stump_dob,
                             NumericVector top_dib,
                             int           simlen   = 200,
                             int           alignlen = 180,
                             bool          detail   = false) {
  NumericVector standage   = st["stand_age"];
  NumericVector totage_max = st["totage_max"];
  int ns = standage.size();

  NumericVector in_sit[4], in_sib[4], in_ta[4], in_bha[4], in_y2[4], in_ba[4], in_psi[4], in_sdf[4], in_den0[4];
  for (int s = 0; s < 4; s++) {
    std::string p = PRE[s];
    in_sit[s]  = st[p + "si_t"];
    in_sib[s]  = st[p + "si_bh"];
    in_ta[s]   = st[p + "totage"];
    in_bha[s]  = st[p + "bhage"];
    in_y2[s]   = st[p + "y2bh"];
    in_ba[s]   = st[p + "ba"];
    in_psi[s]  = st[p + "psi"];
    in_sdf[s]  = st[p + "sdf"];
    in_den0[s] = st[p + "den0"];
  }

  // ---- Output columns -------------------------------------------------------
  int nr = ns*simlen;
  IntegerVector o_idx(nr), o_sect(nr);
  NumericVector o_totage(nr), o_standage(nr);
  NumericVector o_ta[4], o_bha[4], o_ht[4], o_n[4], o_sc[4], o_ba[4], o_qmd[4],
                o_tv[4], o_mv[4], o_mn[4], o_mait[4], o_maim[4],
                o_scb[4], o_ba1[4], o_ba2[4];
  for (int s = 0; s < 4; s++) {
    o_ta[s] = NumericVector(nr); o_bha[s] = NumericVector(nr); o_ht[s]   = NumericVector(nr);
    o_n[s]  = NumericVector(nr); o_sc[s]  = NumericVector(nr); o_ba[s]   = NumericVector(nr);
    o_qmd[s]= NumericVector(nr); o_tv[s]  = NumericVector(nr); o_mv[s]   = NumericVector(nr);
    o_mn[s] = NumericVector(nr); o_mait[s]= NumericVector(nr); o_maim[s] = NumericVector(nr);
    if (detail) {
      o_scb[s] = NumericVector(nr); o_ba1[s] = NumericVector(nr); o_ba2[s] = NumericVector(nr);
    }
  }
  NumericVector o_tv_con(nr), o_mv_con(nr), o_tv_all(nr), o_mv_all(nr),
                o_mait_all(nr), o_maim_all(nr), o_mait_con(nr), o_maim_con(nr),
                o_mait_dec(nr), o_maim_dec(nr);

  // Per-stand initial values
  NumericVector i_t0(ns), i_sc0[4], i_ba0[4], i_rate[4];
  for (int s = 0; s < 4; s++) {
    i_sc0[s] = NumericVector(ns); i_ba0[s] = NumericVector(ns); i_rate[s] = NumericVector(ns);
  }

  // Yearly state for the current stand (row r, species s -> [r*4 + s])
  std::vector<double> TA1(simlen*4), BHA1(simlen*4), HT(simlen*4), DEN2(simlen*4),
                      SC(simlen*4), SCB(simlen*4), BA1(simlen*4), BA2(simlen*4), T1(simlen);

  // SAS lag1() runs across stands: the first year of a stand sees the last
  // year of the previous stand (normally zeroed because BH age <= 0)
  double prev_ba1[4] = {NaN, NaN, NaN, NaN}, prev_ba2[4] = {NaN, NaN, NaN, NaN};

  for (int i = 0; i < ns; i++) {
    Rcpp::checkUserInterrupt();

    double sit[4], sib[4], ta[4], bha[4], y2[4], baobs[4], psi[4], sdf[4], den0[4];
    for (int s = 0; s < 4; s++) {
      sit[s]  = in_sit[s][i];  sib[s]  = in_sib[s][i];  ta[s]   = in_ta[s][i];
      bha[s]  = in_bha[s][i];  y2[s]   = in_y2[s][i];   baobs[s]= in_ba[s][i];
      psi[s]  = in_psi[s][i];  sdf[s]  = in_sdf[s][i];  den0[s] = in_den0[s][i];
    }
    double tmax = totage_max[i];

    // ---- examp6: age offsets ------------------------------------------------
    double diff[4], bh0[4];
    for (int s = 0; s < 4; s++) {
      diff[s] = (sv(ta[s]) > 0) ? tmax - ta[s] : NaN;                    // years younger than leading species
      bh0[s]  = miss(bha[s]) ? NaN : bha[s] - (std::ceil(bha[s]) - 1.0); // fractional BH age, (0, 1]
    }
    double temp1 = tmax - std::ceil(tmax);
    double t0    = (temp1 != 0) ? tmax - (std::ceil(tmax) - 1.0) : 0.0;   // fractional leading total age

    // ---- examp7 - examp11: species composition at BH age 0 -------------------
    double sc0[4] = {NaN, NaN, NaN, NaN};
    {
      double ta1[4]  = {NaN, NaN, NaN, NaN};
      double bha1[4] = {NaN, NaN, NaN, NaN};
      double t1 = t0 - 1.0;
      for (int k = 0; k < alignlen; k++) {
        t1 = (k == 0) ? t0 : t1 + 1.0;   // SAS: age1 = age1 + 1
        for (int s = 0; s < 4; s++) if (sv(diff[s]) >= 0) ta1[s] = t1 - diff[s];
        for (int s = 0; s < 4; s++) if (sv(ta1[s]) >= 0)  bha1[s] = ta1[s] - y2[s];

        double bh2[4];   // integer BH ages that pass through 0 exactly
        for (int s = 0; s < 4; s++) bh2[s] = miss(bha1[s]) ? NaN : std::round(bha1[s] - bh0[s]);

        for (int f = 0; f < 4; f++) {
          if (!(sv(sdf[f]) > 0 && bh2[f] == 0.0)) continue;

          // Other species' total ages when species f is at BH age 0
          double ta2[4];
          for (int o = 0; o < 4; o++)
            ta2[o] = (miss(ta1[o]) || miss(bh0[f])) ? NaN : ta1[o] - bh0[f];

          double den[4];
          for (int o = 0; o < 4; o++) {
            if (o == AW) {   // AW density runs on BH age
              double a = (f == AW) ? bh2[AW]
                                   : ((miss(ta2[AW]) || miss(y2[AW])) ? NaN : ta2[AW] - y2[AW]);
              den[o] = (sv(a) >= 0 && sv(sdf[AW]) > 0) ? Density(AW, a, sdf[AW], sib[AW], psi[AW], sdf) : NaN;
            } else {
              den[o] = (sv(ta2[o]) >= 0 && sv(sdf[o]) > 0) ? Density(o, ta2[o], sdf[o], sib[o], psi[o], sdf) : NaN;
            }
          }
          // examp10: species not yet at breast height count as 0
          for (int o = 0; o < 4; o++) if (o != f && sv(bh2[o]) < 0) den[o] = 0.0;

          if (sv(den[f]) >= 0) sc0[f] = den[f]/sas_sum4(den);
        }
      }
    }
    for (int s = 0; s < 4; s++) sc0[s] = nz(sc0[s]);

    // ---- examp13: BA at the fractional BH age --------------------------------
    // BAinc at BH age 0 and BA 0 times the fractional BH age
    double ba0[4] = {NaN, NaN, NaN, NaN};
    double zero[4] = {0.0, 0.0, 0.0, 0.0};
    for (int s = AW; s <= SB; s++)
      if (sv(den0[s]) > 0 && sv(sib[s]) > 0)
        ba0[s] = bh0[s]*BAinc(s, 0.0, sc0[s], den0[s], sib[s], psi[s], 0.0, zero);

    // SW and PL: competitors that reached BH age later are left out
    double sdf1[4] = {nz(sdf[AW]), nz(sdf[SB]), sdf[SW], nz(sdf[PL])};
    for (int x : {AW, SB, PL}) if (sv(bha[SW]) > sv(bha[x])) sdf1[x] = 0.0;
    if (sv(den0[SW]) > 0 && sv(sib[SW]) > 0)
      ba0[SW] = bh0[SW]*BAinc(SW, 0.0, sc0[SW], den0[SW], sib[SW], psi[SW], 0.0, sdf1);

    sdf1[AW] = nz(sdf[AW]); sdf1[SB] = nz(sdf[SB]); sdf1[SW] = nz(sdf[SW]); sdf1[PL] = sdf[PL];
    for (int x : {AW, SB, SW}) if (sv(bha[PL]) > sv(bha[x])) sdf1[x] = 0.0;
    if (sv(den0[PL]) > 0 && sv(sib[PL]) > 0)
      ba0[PL] = bh0[PL]*BAinc(PL, 0.0, sc0[PL], den0[PL], sib[PL], psi[PL], 0.0, sdf1);
    // sdf1 carries into the yearly loop, as SDF_xx1 does in the SAS

    // ---- examp14: yearly projection ------------------------------------------
    // BA1 = predicted BA from BA at the fractional BH age
    // BA2 = predicted BA from observed BA, starting at the observed age
    double ba1[4], ba2[4];
    for (int s = 0; s < 4; s++) { ba1[s] = ba0[s]; ba2[s] = baobs[s]; }

    double ta1[4]  = {NaN, NaN, NaN, NaN}, bha1[4] = {NaN, NaN, NaN, NaN},
           ht[4]   = {NaN, NaN, NaN, NaN}, den2[4] = {NaN, NaN, NaN, NaN},
           sc[4]   = {NaN, NaN, NaN, NaN}, scb[4]  = {NaN, NaN, NaN, NaN};
    double age1 = t0;

    auto grow = [&](int s) {
      double inc1 = BAinc(s, bha1[s], scb[s], den0[s], sib[s], psi[s], ba1[s], sdf1);
      ba1[s] += inc1;
      if (sv(ba1[s]) < 0 && !miss(inc1)) ba1[s] = 0.0;

      if (sv(sas_round(ta1[s], 0.1)) >= sv(sas_round(ta[s], 0.1)) && sv(ba2[s]) > 0) {
        double inc2 = BAinc(s, bha1[s], scb[s], den0[s], sib[s], psi[s], ba2[s], sdf1);
        ba2[s] += inc2;
        if (sv(ba2[s]) < 0 && !miss(inc2)) ba2[s] = 0.0;
      }
    };

    for (int k = 0; k < simlen; k++) {
      double t1 = age1;
      for (int s = 0; s < 4; s++) if (sv(diff[s]) >= 0) ta1[s] = t1 - diff[s];
      for (int s = 0; s < 4; s++) if (sv(ta1[s]) >= 0)  bha1[s] = ta1[s] - y2[s];

      // Top height and density
      for (int s = 0; s < 4; s++) {
        if (!(sv(ta1[s]) >= 0)) continue;
        ht[s] = Ht(s, sit[s], ta1[s]);
        if (s == SW) sdf1[AW] = nz(sdf[AW]);
        if (s == PL) { sdf1[AW] = nz(sdf[AW]); sdf1[SB] = nz(sdf[SB]); sdf1[SW] = nz(sdf[SW]); sdf1[PL] = sdf[PL]; }
        if (s == AW) {
          if (sv(sdf[AW]) > 0 && sv(bha1[AW]) > 0) den2[AW] = Density(AW, bha1[AW], sdf[AW], sib[AW], psi[AW], sdf);
        } else if (sv(sdf[s]) > 0) {
          den2[s] = Density(s, ta1[s], sdf[s], sib[s], psi[s], sdf);
        }
      }

      // Species composition (total age and BH age based)
      for (int s = 0; s < 4; s++) den2[s] = nz(den2[s]);
      // AW density before BH age is taken as the initial density
      if (den2[AW] == 0 && sv(ta1[AW]) >= 0) den2[AW] = den0[AW];

      double tot = sas_sum4(den2);
      if (sv(tot) > 0) for (int s = 0; s < 4; s++) if (sv(den2[s]) >= 0) sc[s] = den2[s]/tot;

      double den2b[4];
      for (int s = 0; s < 4; s++) den2b[s] = (sv(bha1[s]) < 0) ? 0.0 : den2[s];
      double totb = sas_sum4(den2b);
      if (sv(totb) > 0) for (int s = 0; s < 4; s++) if (sv(den2b[s]) >= 0) scb[s] = den2b[s]/totb;

      // BA increments
      if (sv(bha1[AW]) >= 0 && sv(scb[AW]) > 0 && sv(den0[AW]) > 0) grow(AW);
      if (sv(bha1[SB]) >= 0 && sv(scb[SB]) > 0)                     grow(SB);

      // Competitor SDFs for SW and PL count only once at breast height
      for (int s : {SW, SB, AW, PL}) {
        if (sv(bha1[s]) > 0)      sdf1[s] = (s == PL) ? sdf[PL] : nz(sdf[s]);
        else if (sv(bha1[s]) < 0) sdf1[s] = NaN;
      }
      sdf1[AW] = nz(sdf1[AW]); sdf1[PL] = nz(sdf1[PL]); sdf1[SB] = nz(sdf1[SB]);
      if (sv(bha1[SW]) >= 0 && sv(scb[SW]) > 0) grow(SW);

      sdf1[AW] = nz(sdf1[AW]); sdf1[SW] = nz(sdf1[SW]); sdf1[SB] = nz(sdf1[SB]);
      if (sv(bha1[PL]) >= 0 && sv(scb[PL]) > 0) grow(PL);

      T1[k] = t1;
      age1 += 1.0;

      for (int s = 0; s < 4; s++) {
        int j = k*4 + s;
        TA1[j] = ta1[s]; BHA1[j] = bha1[s]; HT[j]  = ht[s];  DEN2[j] = den2[s];
        SC[j]  = sc[s];  SCB[j]  = scb[s];  BA1[j] = ba1[s]; BA2[j]  = ba2[s];
      }
    }

    // ---- examp14_1 - examp17: localized BA ------------------------------------
    // BA from the previous year (SAS lag1), zero before breast height
    std::vector<double> BA1L(simlen*4), BA2L(simlen*4);
    for (int k = 0; k < simlen; k++) {
      for (int s = 0; s < 4; s++) {
        int j = k*4 + s;
        BA1L[j] = (k > 0) ? BA1[j - 4] : prev_ba1[s];
        BA2L[j] = (k > 0) ? BA2[j - 4] : prev_ba2[s];
        if (sv(BHA1[j]) <= 0) { BA1L[j] = 0.0; BA2L[j] = 0.0; }
      }
    }
    for (int s = 0; s < 4; s++) {
      prev_ba1[s] = BA1[(simlen - 1)*4 + s];
      prev_ba2[s] = BA2[(simlen - 1)*4 + s];
    }

    // Observed / predicted BA at the observed age
    double rate[4] = {NaN, NaN, NaN, NaN};
    for (int s = 0; s < 4; s++) {
      if (!(sv(ta[s]) > 0)) continue;
      double key = sas_round(ta[s], 0.1);
      for (int k = 0; k < simlen; k++) {
        int j = k*4 + s;
        if (sas_round(TA1[j], 0.1) != key) continue;
        if (sv(BA1L[j]) > 0 && sv(baobs[s]) > 0) rate[s] = baobs[s]/BA1L[j];
        if (BA1L[j] == 0)                         rate[s] = 0.0;
        break;
      }
    }

    // ---- Output rows -------------------------------------------------------
    for (int k = 0; k < simlen; k++) {
      int r = i*simlen + k;
      double t1 = T1[k];
      o_idx[r] = i + 1;
      o_sect[r] = k + 1;
      o_totage[r] = t1;
      o_standage[r] = t1 + (standage[i] - tmax);

      double tv[4], mv[4], mait[4], maim[4];
      for (int s = 0; s < 4; s++) {
        int j = k*4 + s;
        double ta1s = TA1[j], bha1s = BHA1[j], hts = HT[j], dens = DEN2[j];

        // examp17: localized BA (later rules override earlier ones)
        double bl = NaN;
        bool before = sv(ta1s) <= sv(ta[s]);   // at or before the observed age
        if (sv(baobs[s]) <= 0) bl = BA1L[j];
        if (sv(ta[s]) > 0 && sv(bha1s) > 0 && before && sv(rate[s]) <= 0 && sv(baobs[s]) > 0)
          bl = baobs[s]/bha[s]*bha1s;
        if (sv(ta[s]) > 0 && sv(BA1L[j]) > 0 && before && sv(rate[s]) > 0 && sv(rate[s]) < RATE_CAP[s])
          bl = rate[s]*BA1L[j];
        if (sv(ta[s]) > 0 && sv(BA1L[j]) > 0 && before && sv(rate[s]) >= RATE_CAP[s] && sv(baobs[s]) > 0)
          bl = baobs[s]/bha[s]*bha1s;
        if (sv(ta[s]) > 0 && sv(BA2L[j]) > 0 && sv(ta1s) > sv(ta[s]))
          bl = BA2L[j];

        // examp18: volumes
        double qmd = (sv(dens) > 0 && !miss(bl)) ? std::sqrt(bl/(dens/10000.0)) : NaN;
        double tvs = NaN, mvs = NaN, mns = NaN, maits = NaN, maims = NaN;
        if (sv(hts) > 1.3 && !miss(bl)) {
          tvs = TotVol(s, bl, hts);
          if (sv(qmd) > 0) {
            MerchResult m = Merch(s, tvs, qmd, hts, bha1s, dens, stump_ht[s], stump_dob[s], top_dib[s]);
            mvs = sas_round(m.mv, 0.00001);
            mns = m.mn;
            if (sv(ta1s) > 0) { maits = tvs/ta1s; maims = mvs/ta1s; }
          }
        }
        tv[s] = nz(tvs); mv[s] = nz(mvs); mait[s] = nz(maits); maim[s] = nz(maims);

        // examp19: blank zeros and pre-breast-height BA
        double scs = SC[j];
        if (scs == 0)  scs = NaN;
        if (dens == 0) dens = NaN;
        if (sv(bha1s) < 0) bl = NaN;

        o_ta[s][r]   = na(ta1s);  o_bha[s][r] = na(bha1s); o_ht[s][r]   = na(hts);
        o_n[s][r]    = na(dens);  o_sc[s][r]  = na(scs);   o_ba[s][r]   = na(bl);
        o_qmd[s][r]  = na(qmd);   o_tv[s][r]  = tv[s];     o_mv[s][r]   = mv[s];
        o_mn[s][r]   = na(mns);   o_mait[s][r]= mait[s];   o_maim[s][r] = maim[s];
        if (detail) {
          o_scb[s][r] = na(SCB[j]); o_ba1[s][r] = na(BA1[j]); o_ba2[s][r] = na(BA2[j]);
        }
      }

      double tvc = tv[SB] + tv[SW] + tv[PL], mvc = mv[SB] + mv[SW] + mv[PL];
      o_tv_con[r] = tvc;           o_mv_con[r] = mvc;
      o_tv_all[r] = tvc + tv[AW];  o_mv_all[r] = mvc + mv[AW];
      bool pos = t1 > 0;
      o_mait_all[r] = pos ? (tvc + tv[AW])/t1 : NA_REAL;
      o_maim_all[r] = pos ? (mvc + mv[AW])/t1 : NA_REAL;
      o_mait_con[r] = pos ? tvc/t1    : NA_REAL;
      o_maim_con[r] = pos ? mvc/t1    : NA_REAL;
      o_mait_dec[r] = pos ? tv[AW]/t1 : NA_REAL;
      o_maim_dec[r] = pos ? mv[AW]/t1 : NA_REAL;
    }

    i_t0[i] = t0;
    for (int s = 0; s < 4; s++) {
      i_sc0[s][i] = sc0[s]; i_ba0[s][i] = na(ba0[s]); i_rate[s][i] = na(rate[s]);
    }
  }

  // ---- Assemble -------------------------------------------------------------
  const char* SPP[4] = {"aw", "sb", "sw", "pl"};
  List rows;
  rows["stand_idx"] = o_idx;
  rows["year"]      = o_sect;
  rows["totage"]    = o_totage;
  rows["standage"]  = o_standage;
  for (int s = 0; s < 4; s++) {
    std::string x = SPP[s];
    rows["totage_" + x] = o_ta[s];
    rows["bhage_"  + x] = o_bha[s];
    rows["topht_"  + x] = o_ht[s];
    rows["n_"      + x] = o_n[s];
    rows["sc_"     + x] = o_sc[s];
    rows["ba_"     + x] = o_ba[s];
    rows["qmd_"    + x] = o_qmd[s];
    rows["tvol_"   + x] = o_tv[s];
    rows["mvol_"   + x] = o_mv[s];
    rows["mn_"     + x] = o_mn[s];
    rows["mai_t_"  + x] = o_mait[s];
    rows["mai_m_"  + x] = o_maim[s];
    if (detail) {
      rows["sc_b_" + x] = o_scb[s];
      rows["ba1_"  + x] = o_ba1[s];
      rows["ba2_"  + x] = o_ba2[s];
    }
  }
  rows["tvol_con"]  = o_tv_con;   rows["mvol_con"]  = o_mv_con;
  rows["tvol_all"]  = o_tv_all;   rows["mvol_all"]  = o_mv_all;
  rows["mai_t_con"] = o_mait_con; rows["mai_m_con"] = o_maim_con;
  rows["mai_t_dec"] = o_mait_dec; rows["mai_m_dec"] = o_maim_dec;
  rows["mai_t_all"] = o_mait_all; rows["mai_m_all"] = o_maim_all;

  List init;
  init["totage_0"] = i_t0;
  for (int s = 0; s < 4; s++) {
    std::string x = SPP[s];
    init[x + "_sc0"]  = i_sc0[s];
    init[x + "_ba0"]  = i_ba0[s];
    init[x + "_rate"] = i_rate[s];
  }

  return List::create(_["rows"] = rows, _["init"] = init);
}
