# Validate GYPSY_Sim() against the full SAS test run
#
# Runs every stand in test_gypsy_bp.xlsx through GYPSY_Sim() and compares the
# result with each column of for_gypsy_test_bp.csv (output of GYPSY_run.sas).
#
# Each SAS column is either
#   direct   - the same value is returned by GYPSY_Sim() or attr(, "stand_init")
#   derived  - rebuilt here from GYPSY_Sim() output (SAS intermediates such as
#              the lagged BA or predicted density at the observed age)
#   skipped  - SAS scratch variables that hold whatever the last species in the
#              row computed (CUM1, a, b, c, ...), or fixed settings
#
# A value passes when both sides are missing, or both are present and
# |SAS - R| <= rel_tol*|SAS| + 1e-8. Merchantable volume is rounded to 1e-5 in
# the SAS, so columns built on it may also be one rounding step apart.
#
# Run from the package root. Uses roughly 4-5 GB of memory.
# Writes the per-column report to scratch/data/sas_validation_report.csv.

library(data.table)
library(openxlsx)
devtools::load_all(".")   # or library(rGYPSY) to test the installed version

data_dir <- "scratch/data"
rel_tol  <- 1e-6

spp_l <- c("aw", "sb", "sw", "pl")
spp_u <- toupper(spp_l)

# ---- Run GYPSY_Sim -------------------------------------------------------------
inp <- as.data.table(read.xlsx(file.path(data_dir, "test_gypsy_bp.xlsx")))

run_sim <- function(stump_dob, top_dib) {
  suppressWarnings(GYPSY_Sim(inp$Stand, inp$Standage, inp$Species, inp$N,
                             si = inp$SI, age = inp$totage, ht = inp$TopH,
                             ps = inp$PS, ba = inp$BA,
                             stump_ht = 0.3, stump_dob = stump_dob, top_dib = top_dib,
                             detail = TRUE))
}

t0 <- Sys.time()
r1 <- run_sim(13, 7)    # SAS standard 1: StumpHT1 = 0.3, StumpDOB1 = 13, TopDIB1 = 7
r2 <- run_sim(15, 10)   # SAS standard 2: StumpHT2 = 0.3, StumpDOB2 = 15, TopDIB2 = 10
cat(sprintf("GYPSY_Sim: %d stands, %d rows per run, %.1f s for both runs\n",
            uniqueN(r1$stand), nrow(r1), as.numeric(Sys.time() - t0, units = "secs")))
init <- attr(r1, "stand_init")
setorder(r1, stand, year)
setorder(r2, stand, year)

# ---- Stand-level values with SAS names --------------------------------------------
# (SAS repeats these on every row; they are compared on the first year)
rs <- data.table(stand = init$stand)
src_s <- c()   # SAS name -> "direct" / "derived"
put_s <- function(name, value, how = "direct") {
  set(rs, j = name, value = value)
  src_s[[name]] <<- how
}

put_s("standage",   init$stand_age)
put_s("totage_max", init$totage_max)
put_s("totage_0",   init$totage_0)
put_s("temp1",      init$totage_max - ceiling(init$totage_max), "derived")

# Leading species: SAS sorts rows by stand and total age (stable, missing first)
# and takes the last row, so ties go to the species entered last
long <- rbindlist(lapply(1:4, function(k)
  data.table(stand = init$stand, Species = spp_u[k], totage = init[[paste0(spp_l[k], "_totage")]])))
long <- long[inp[, .(stand = Stand, Species, row = .I)], on = .(stand, Species)]
setorder(long, stand, totage, row, na.last = FALSE)
lead <- long[, .(species_max = Species[.N]), by = stand]
put_s("species_max", lead[rs, on = "stand", species_max], "derived")

for (k in 1:4) {
  s <- spp_l[k]; S <- spp_u[k]
  g <- function(v) init[[paste0(s, "_", v)]]

  # Inputs and SI / age (SAS data v1, examp3)
  put_s(paste0("si_bh_", s),  g("si_bh"))
  put_s(paste0("si_t_", s),   g("si_t"))
  put_s(paste0("bhage_", s),  g("bhage"))
  put_s(paste0("totage_", s), g("totage"))
  put_s(paste0("y2bh_", s),   g("y2bh"))
  put_s(paste0("den_", s),    g("n"))
  put_s(paste0("ps_", s),     g("ps"))
  put_s(paste0("ba_", s),     g("ba"))
  put_s(paste0("toph_", s),   g("ht"))

  # PSI, SDF and initial density (examp4, examp5, examp13)
  put_s(paste0("psi_", s), g("psi"))
  put_s(paste0("sdf_", s), g("sdf"))
  put_s(if (s == "aw") "den_aw_bhage0" else paste0("den_", s, "_totage0"), g("den0"))
  put_s(paste0("den_", s, "0"), g("den0"))

  # Predicted density at the observed age (AW on BH age), only where SDF was solved
  obs_age <- if (s == "aw") g("bhage") else g("totage")
  pred <- GYPSY_N(rep(S, nrow(init)), obs_age, g("sdf"), g("si_bh"), nafill(g("psi"), fill = 0),
                  sdf_aw = init$aw_sdf, sdf_sw = init$sw_sdf, sdf_sb = init$sb_sdf)
  pred[is.na(g("den0"))] <- NA
  put_s(paste0("den_", s, "_pred"), pred, "derived")

  # Age offsets (examp6)
  put_s(paste0("diff_totage_", s),
        ifelse(!is.na(g("totage")) & g("totage") > 0, init$totage_max - g("totage"), NA_real_), "derived")
  bh0 <- g("bhage") - (ceiling(g("bhage")) - 1)
  put_s(paste0("bhage_", s, "0"), bh0, "derived")

  # Species composition at BH age 0, starting BA, BA ratio (examp10 - examp15)
  put_s(paste0("sc_", s, "0"), g("sc0"))
  put_s(paste0("ba_", s, "_0"), g("ba0"))
  put_s(paste0("bainc0_", s), g("ba0")/bh0, "derived")
  put_s(paste0("rate_", s), g("rate"))
}

# ---- Yearly values with SAS names -------------------------------------------------
ry <- r1[, .(stand, sect = year)]
src_y <- c()
put_y <- function(name, value, how = "direct") {
  set(ry, j = name, value = value)
  src_y[[name]] <<- how
}

put_y("totage1",   r1$totage)
put_y("standage1", r1$standage)
put_y("age1",      r1$totage + 1, "derived")   # SAS age1 is already incremented when output

# SAS lag1(): previous row, across stands, zero before breast height
lag_ba <- function(x, bh) { y <- shift(x); y[is.na(bh) | bh <= 0] <- 0; y }

sc_cols <- sc_b_cols <- list()
for (k in 1:4) {
  s <- spp_l[k]; S <- spp_u[k]
  a <- function(v, r = r1) r[[paste0(v, "_", s)]]

  put_y(paste0("totage1_", s), a("totage"))
  put_y(paste0("bhage1_", s),  a("bhage"))
  put_y(paste0("topht_", s),   a("topht"))
  put_y(paste0("den2_", s),    a("n"))
  put_y(paste0("sc_", s),      a("sc"))
  put_y(paste0("sc_", s, "_b"), a("sc_b"))
  put_y(paste0("ba1_", s),     a("ba1"))
  put_y(paste0("ba2_", s),     a("ba2"))
  put_y(paste0("ba_loc_", s),  a("ba"))
  put_y(paste0("tv_", s, "p"), a("tvol"))
  put_y(paste0("mv1_", s),     a("mvol"))
  put_y(paste0("mv2_", s),     a("mvol", r2))
  put_y(paste0("mn1_", s),     a("mn"))
  put_y(paste0("mn2_", s),     a("mn", r2))
  put_y(paste0("mai00_", s),   a("mai_t"))
  put_y(paste0("mai1_", s),    a("mai_m"))
  put_y(paste0("mai2_", s),    a("mai_m", r2))

  # Derived intermediates (examp14, examp14_1)
  bh <- a("bhage")
  d2 <- nafill(a("n"), fill = 0)   # examp19 blanks 0 densities; SAS kept 0 here
  put_y(paste0("den2_", s, "_b"), ifelse(is.na(bh) | bh < 0, 0, d2), "derived")
  put_y(paste0("ba1_", s, "_1"), lag_ba(a("ba1"), bh), "derived")
  put_y(paste0("ba2_", s, "_1"), lag_ba(a("ba2"), bh), "derived")

  sc_cols[[s]]   <- a("sc")
  sc_b_cols[[s]] <- a("sc_b")
}
row_sum <- function(l) {
  m <- do.call(cbind, l)
  ifelse(rowSums(!is.na(m)) == 0, NA_real_, rowSums(m, na.rm = TRUE))
}
put_y("sc_all",   row_sum(sc_cols),   "derived")
put_y("sc_all_b", row_sum(sc_b_cols), "derived")

put_y("v_con_00",  r1$tvol_con)
put_y("v_all_00",  r1$tvol_all)
put_y("v_con_mv1", r1$mvol_con);  put_y("v_con_mv2", r2$mvol_con)
put_y("v_all_mv1", r1$mvol_all);  put_y("v_all_mv2", r2$mvol_all)
put_y("mai00_all", r1$mai_t_all); put_y("mai1_all", r1$mai_m_all); put_y("mai2_all", r2$mai_m_all)
put_y("mai00_con", r1$mai_t_con); put_y("mai1_con", r1$mai_m_con); put_y("mai2_con", r2$mai_m_con)
put_y("mai00_dec", r1$mai_t_dec); put_y("mai1_dec", r1$mai_m_dec); put_y("mai2_dec", r2$mai_m_dec)
rm(r1, r2); invisible(gc())

# SAS columns not compared, with the reason
skip <- c(
  stumpht1 = "setting (0.3)", stumpdob1 = "setting (13)", topdib1 = "setting (7)",
  stumpht2 = "setting (0.3)", stumpdob2 = "setting (15)", topdib2 = "setting (10)",
  a = "SAS scratch: stump DOB correction coefficient of the last species in the row",
  b = "SAS scratch: stump DOB correction coefficient of the last species in the row",
  c = "SAS scratch: stump DOB correction coefficient of the last species in the row",
  a0 = "SAS scratch: Weibull coefficient of the last species in the row",
  c0 = "SAS scratch: Weibull coefficient of the last species in the row",
  c1 = "SAS scratch: Weibull coefficient of the last species in the row",
  c2 = "SAS scratch: Weibull coefficient of the last species in the row",
  stump03_1 = "SAS scratch: stump DOB at 0.3 m, last species in the row",
  stump03_2 = "SAS scratch: stump DOB at 0.3 m, last species in the row",
  cum1 = "SAS scratch: % stems below merch DBH, last species in the row",
  cum2 = "SAS scratch: % stems below merch DBH, last species in the row"
)
for (s in spp_l) {
  skip[paste0("sdf_", s, "1")] <- "SAS scratch: competitor SDF used by the last BA increment"
  skip[paste0(c("bainc1_", "bainc2_"), s)] <- "not returned (BA increment of the last year it was computed)"
  skip[paste0(c("rate1_", "rate2_"), s)] <- "not returned (stump DOB height correction, 1 at 0.3 m)"
  skip[paste0(c("dbh1_", "dbh2_"), s)] <- "not returned (merchantable DBH from stump DOB)"
}

# ---- Read the SAS output ----------------------------------------------------------
csv <- file.path(data_dir, "for_gypsy_test_bp.csv")
hdr <- names(fread(csv, nrows = 0))
hdr_l <- tolower(hdr)
keep <- hdr[hdr_l %in% c("stand", "sect", names(src_s), names(src_y))]

t0 <- Sys.time()
sas <- fread(csv, select = keep, na.strings = c("", "."))
setnames(sas, tolower(names(sas)))
setorder(sas, stand, sect)
cat(sprintf("Read %d SAS rows x %d columns in %.1f s\n", nrow(sas), ncol(sas),
            as.numeric(Sys.time() - t0, units = "secs")))

# ---- Compare ----------------------------------------------------------------------
compare <- function(sas_dt, r_dt, keys, src) {
  m <- merge(sas_dt[, c(keys, intersect(names(src), names(sas_dt))), with = FALSE],
             r_dt, by = keys, suffixes = c(".sas", ".r"))
  rbindlist(lapply(intersect(names(src), names(sas_dt)), function(v) {
    x <- m[[paste0(v, ".sas")]]; y <- m[[paste0(v, ".r")]]
    if (is.character(x) || is.character(y)) {
      bad <- !((is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y))
      return(data.table(column = v, how = src[[v]], n = sum(!is.na(x) & !is.na(y)),
                        na_mismatch = sum(is.na(x) != is.na(y)), max_abs = NA_real_, max_rel = NA_real_,
                        n_fail = sum(bad), worst_stand = m$stand[which(bad)[1]],
                        worst_sect = if ("sect" %in% keys) m$sect[which(bad)[1]] else NA_integer_))
    }
    x <- as.numeric(x); y <- as.numeric(y)
    both <- !is.na(x) & !is.na(y)
    d <- abs(x - y); d[!both] <- NA
    lim <- rel_tol*abs(x) + 1e-8
    if (grepl("^(mv[12]|mai[12]|v_(con|all)_mv)", v)) lim <- lim + 1.01e-5
    bad <- (is.na(x) != is.na(y)) | (both & d > lim)
    rel <- d/pmax(abs(x), 1e-8)
    w <- if (any(bad)) which(bad)[which.max(replace(d[bad], is.na(d[bad]), Inf))] else NA_integer_
    data.table(column = v, how = src[[v]], n = sum(both),
               na_mismatch = sum(is.na(x) != is.na(y)),
               max_abs = if (any(both)) max(d, na.rm = TRUE) else NA_real_,
               max_rel = if (any(both)) max(rel, na.rm = TRUE) else NA_real_,
               n_fail = sum(bad),
               worst_stand = if (is.na(w)) NA else m$stand[w],
               worst_sect  = if (is.na(w) || !"sect" %in% keys) NA_integer_ else m$sect[w])
  }))
}

res_s <- compare(sas[sect == 1], rs, "stand", src_s)
res_y <- compare(sas, ry, c("stand", "sect"), src_y)
res <- rbind(cbind(level = "stand", res_s), cbind(level = "yearly", res_y))
res[, pass := n_fail == 0]

skipped  <- data.table(column = hdr_l[hdr_l %in% names(skip)], reason = skip[hdr_l[hdr_l %in% names(skip)]])
unmapped <- setdiff(hdr_l, c("stand", "sect", res$column, skipped$column))

# ---- Report -------------------------------------------------------------------
cat(sprintf("\nSAS columns: %d | compared: %d (%d direct, %d derived) | skipped: %d | unmapped: %d\n",
            length(hdr), nrow(res), sum(res$how == "direct"), sum(res$how == "derived"),
            nrow(skipped), length(unmapped)))
cat(sprintf("Passed: %d of %d compared columns\n\n", sum(res$pass), nrow(res)))

if (any(!res$pass)) {
  cat("Columns with differences:\n")
  print(res[pass == FALSE], nrows = Inf)
} else {
  cat("All compared columns match.\n")
}
cat("\nLargest relative differences:\n")
print(res[order(-max_rel)][1:10, .(level, column, how, n, max_abs, max_rel, n_fail)])

if (length(unmapped)) cat("\nUnmapped SAS columns:", paste(unmapped, collapse = ", "), "\n")
cat("\nSkipped SAS columns:\n"); print(skipped, nrows = Inf)

fwrite(res, file.path(data_dir, "sas_validation_report.csv"))
cat("\nFull report:", file.path(data_dir, "sas_validation_report.csv"), "\n")
