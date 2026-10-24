# Helpers for comparing GYPSY_Sim with the SAS test run
# (fixture built by data-raw/make_test_fixture.R)

sas_fixture <- function() readRDS(test_path("fixtures", "sas_subset.rds"))

# Run GYPSY_Sim on the fixture inputs with one merchantability standard
run_sas_inputs <- function(inp, stump_dob, top_dib) {
  suppressWarnings(GYPSY_Sim(inp$Stand, inp$Standage, inp$Species, inp$N,
                             si = inp$SI, age = inp$totage, ht = inp$TopH,
                             ps = inp$PS, ba = inp$BA,
                             stump_ht = 0.3, stump_dob = stump_dob, top_dib = top_dib,
                             detail = TRUE))
}

.spp <- c("aw", "sb", "sw", "pl")

# SAS column (lower case) -> GYPSY_Sim column
stand_map <- function() {
  m <- c(totage_max = "totage_max", totage_0 = "totage_0",
         den_aw_bhage0 = "aw_den0", den_sb_totage0 = "sb_den0",
         den_sw_totage0 = "sw_den0", den_pl_totage0 = "pl_den0")
  for (s in .spp) {
    m[paste0(c("si_t_", "si_bh_", "totage_", "bhage_", "y2bh_"), s)] <-
      paste0(s, c("_si_t", "_si_bh", "_totage", "_bhage", "_y2bh"))
    m[paste0(c("psi_", "sdf_", "sc_", "ba_", "rate_"), s, c("", "", "0", "_0", ""))] <-
      paste0(s, c("_psi", "_sdf", "_sc0", "_ba0", "_rate"))
  }
  m
}

yearly_map <- function(standard = 1) {
  if (standard == 1) {
    m <- c(totage1 = "totage", standage1 = "standage",
           v_con_00 = "tvol_con", v_all_00 = "tvol_all",
           mai00_all = "mai_t_all", mai00_con = "mai_t_con", mai00_dec = "mai_t_dec")
    for (s in .spp) {
      m[paste0(c("totage1_", "bhage1_", "topht_", "den2_", "sc_", "ba_loc_", "tv_", "mai00_"),
               s, c("", "", "", "", "", "", "p", ""))] <-
        paste0(c("totage_", "bhage_", "topht_", "n_", "sc_", "ba_", "tvol_", "mai_t_"), s)
      m[paste0(c("sc_", "ba1_", "ba2_"), s, c("_b", "", ""))] <- paste0(c("sc_b_", "ba1_", "ba2_"), s)
    }
  } else m <- c()
  k <- standard
  m[paste0(c("v_con_mv", "v_all_mv", "mai", "mai", "mai"), k, c("", "", "_all", "_con", "_dec"))] <-
    c("mvol_con", "mvol_all", "mai_m_all", "mai_m_con", "mai_m_dec")
  for (s in .spp)
    m[paste0(c("mv", "mn", "mai"), k, "_", s)] <- paste0(c("mvol_", "mn_", "mai_m_"), s)
  m
}

# Join SAS and R rows; SAS value columns get a "sas." prefix so names never clash
align_sas <- function(sas, r, by_sas, by_r) {
  sas <- data.table::copy(sas)
  vals <- setdiff(names(sas), by_sas)
  data.table::setnames(sas, vals, paste0("sas.", vals))
  merge(sas, r, by.x = by_sas, by.y = by_r)
}

# Compare SAS and R columns: same missing pattern, values within 1e-6 relative.
# SAS writes merchantable volume rounded to 1e-5, so values built on it may be
# one rounding step apart.
expect_sas_match <- function(m, map, rel_tol = 1e-6) {
  for (a in names(map)) {
    x <- m[[paste0("sas.", a)]]
    y <- m[[map[[a]]]]
    expect_false(is.null(x), info = paste("SAS column", a))
    expect_false(is.null(y), info = paste("R column", map[[a]]))
    if (is.null(x) || is.null(y)) next

    expect_identical(is.na(x), is.na(y), info = paste(a, "missing pattern"))
    ok  <- !is.na(x) & !is.na(y)
    lim <- rel_tol*abs(x[ok]) + 1e-8
    if (grepl("^(mv|mai[12]|v_(con|all)_mv)", a)) lim <- lim + 1.01e-5
    bad <- sum(abs(x[ok] - y[ok]) > lim)
    expect_equal(bad, 0L, info = paste(a, "values"))
  }
}
