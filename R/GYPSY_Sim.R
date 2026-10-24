# Project stand with GYPSY
# Follows GYPSY_run.sas (MAY 21, 2009 functions)
# Variables:
#   SI, Age, TopHt, BA, N

.gypsy_spp <- c("AW", "SB", "SW", "PL")

# Recycle length-1 inputs to n, otherwise require length n
.recycle <- function(x, n, name){
  if(length(x) == 1) return(rep(x, n))
  if(length(x) != n) stop(name, " must be length 1 or length ", n)
  x
}

# Per-species value from a scalar or a vector named by species
.per_spp <- function(x, name){
  if(length(x) == 1) return(rep(as.numeric(x), 4))
  if(is.null(names(x)) || !all(.gypsy_spp %in% toupper(names(x))))
    stop(name, " must be length 1 or named with all of AW, SB, SW, PL")
  names(x) <- toupper(names(x))
  as.numeric(x[.gypsy_spp])
}

# ---- Stand-level components -------------------------------------------------

#' Percent stocking index
#'
#' Percent stocking at total age 50 for each species in a stand, from the
#' GYPSY percent stocking models. Takes one row per stand with species in
#' columns, because SW depends on AW and PL depends on AW, SB and SW.
#'
#' @param stand Stand identifiers (used for length only).
#' @param aw_ps,sb_ps,sw_ps,pl_ps Observed percent stocking (NA if unknown).
#' @param aw_age_t,sb_age_t,sw_age_t,pl_age_t Total age.
#' @param aw_si_b,sb_si_b,sw_si_b,pl_si_b Site index at breast height age 50.
#'
#' @return A list with `aw_psi`, `sb_psi`, `sw_psi` and `pl_psi`, capped at
#'   100. As in the SAS, AW, SB and SW are 0 when they can't be computed, PL is
#'   `NA`.
#' @seealso [GYPSY_SDF()], [GYPSY_Sim()]
#' @export
#' @examples
#' GYPSY_PSI(1, aw_ps = 60, aw_age_t = 30, aw_si_b = 16,
#'           sb_ps = NA, sb_age_t = NA, sb_si_b = NA,
#'           sw_ps = 40, sw_age_t = 25, sw_si_b = 15,
#'           pl_ps = NA, pl_age_t = NA, pl_si_b = NA)
GYPSY_PSI <- function(stand,
                      aw_ps, # First function that requires data widening
                      aw_age_t,
                      aw_si_b,
                      sb_ps,
                      sb_age_t,
                      sb_si_b,
                      sw_ps,
                      sw_age_t,
                      sw_si_b,
                      pl_ps,
                      pl_age_t,
                      pl_si_b){
  GYPSY_PSI_cpp(as.character(stand),
                as.numeric(aw_ps), as.numeric(aw_age_t), as.numeric(aw_si_b),
                as.numeric(sb_ps), as.numeric(sb_age_t), as.numeric(sb_si_b),
                as.numeric(sw_ps), as.numeric(sw_age_t), as.numeric(sw_si_b),
                as.numeric(pl_ps), as.numeric(pl_age_t), as.numeric(pl_si_b))
}

#' Stand density factor
#'
#' Solves the GYPSY spatial density functions for each species' stand density
#' factor (SDF, density at age 50) from its observed density, and returns the
#' initial density implied by that SDF. Takes one row per stand with species
#' in columns, because SW depends on AW and PL depends on AW, SW and SB.
#'
#' AW is computed only when its density is over 20 stems/ha, other species
#' when density is over 0, as in the SAS.
#'
#' @param stand Stand identifiers (used for length only).
#' @param aw_age_b AW breast height age (AW density runs on breast height
#'   age).
#' @param sb_age_t,sw_age_t,pl_age_t Total age.
#' @param aw_n,sb_n,sw_n,pl_n Observed density (stems/ha).
#' @param aw_si_b,sb_si_b,sw_si_b,pl_si_b Site index at breast height age 50.
#' @param aw_psi,sb_psi,sw_psi,pl_psi Percent stocking index from
#'   [GYPSY_PSI()] (0 or NA when unknown).
#' @param iter Maximum iterations.
#' @param tol Convergence tolerance.
#'
#' @return A list with `aw_sdf`, `sb_sdf`, `sw_sdf`, `pl_sdf` (AW, SB and SW
#'   are 0 when not computed, PL is `NA`) and `aw_den0`, `sb_den0`, `sw_den0`,
#'   `pl_den0`: density at age 0 (AW at breast height age 0, others at total
#'   age 0).
#' @seealso [GYPSY_N()], [GYPSY_PSI()]
#' @export
#' @examples
#' GYPSY_SDF(1, aw_age_b = 30, aw_n = 1500, aw_si_b = 16, aw_psi = 0,
#'           sb_age_t = NA, sb_n = NA, sb_si_b = NA, sb_psi = 0,
#'           sw_age_t = 30, sw_n = 400, sw_si_b = 15, sw_psi = 0,
#'           pl_age_t = NA, pl_n = NA, pl_si_b = NA, pl_psi = NA)
GYPSY_SDF <- function(stand,
                      aw_age_b, # AW uses BH age, others total age
                      aw_n,
                      aw_si_b,
                      aw_psi,
                      sb_age_t,
                      sb_n,
                      sb_si_b,
                      sb_psi,
                      sw_age_t,
                      sw_n,
                      sw_si_b,
                      sw_psi,
                      pl_age_t,
                      pl_n,
                      pl_si_b,
                      pl_psi,
                      iter = 1000,
                      tol  = 1e-8){
  GYPSY_SDF_cpp(as.character(stand),
                as.numeric(aw_age_b), as.numeric(aw_n), as.numeric(aw_si_b), as.numeric(aw_psi),
                as.numeric(sb_age_t), as.numeric(sb_n), as.numeric(sb_si_b), as.numeric(sb_psi),
                as.numeric(sw_age_t), as.numeric(sw_n), as.numeric(sw_si_b), as.numeric(sw_psi),
                as.numeric(pl_age_t), as.numeric(pl_n), as.numeric(pl_si_b), as.numeric(pl_psi),
                as.integer(iter), as.numeric(tol))
}

# ---- Species-level components -------------------------------------------------

#' Top height from site index and age
#'
#' @param spp Character vector of species: `"AW"`, `"SB"`, `"SW"` or `"PL"`.
#' @param si_t Site index at total age 50 (m).
#' @param age_t Total age.
#'
#' @return Top height (m).
#' @seealso [GYPSY_SI()] for the inverse.
#' @export
#' @examples
#' GYPSY_Ht(c("AW", "SB", "SW", "PL"), si_t = rep(15, 4), age_t = rep(40, 4))
GYPSY_Ht <- function(spp,
                     si_t,
                     age_t){
  GYPSY_Ht_cpp(toupper(as.character(spp)),
               as.numeric(si_t),
               as.numeric(age_t))
}

#' Density at age
#'
#' Density (stems/ha) at a given age from a species' stand density factor.
#'
#' @param spp Character vector of species: `"AW"`, `"SB"`, `"SW"` or `"PL"`.
#' @param age Breast height age for AW, total age for other species.
#' @param sdf Stand density factor from [GYPSY_SDF()].
#' @param si_b Site index at breast height age 50.
#' @param psi Percent stocking index (0 when unknown).
#' @param sdf_aw,sdf_sw,sdf_sb Other species' SDFs in the stand (SW uses AW;
#'   PL uses AW, SW and SB); 0 when absent.
#'
#' @return Density (stems/ha); `NA` when `sdf` is not positive. Length-1
#'   arguments are recycled to the longest input.
#' @seealso [GYPSY_SDF()]
#' @export
#' @examples
#' GYPSY_N("PL", age = c(10, 30, 60), sdf = 1800, si_b = 17)
GYPSY_N <- function(spp,
                    age,
                    sdf,
                    si_b,
                    psi    = 0,
                    sdf_aw = 0,
                    sdf_sw = 0,
                    sdf_sb = 0){
  n <- max(lengths(list(spp, age, sdf, si_b, psi, sdf_aw, sdf_sw, sdf_sb)))
  GYPSY_N_cpp(toupper(as.character(.recycle(spp, n, "spp"))),
              as.numeric(.recycle(age, n, "age")),
              as.numeric(.recycle(sdf, n, "sdf")),
              as.numeric(.recycle(si_b, n, "si_b")),
              as.numeric(.recycle(psi, n, "psi")),
              as.numeric(.recycle(sdf_aw, n, "sdf_aw")),
              as.numeric(.recycle(sdf_sw, n, "sdf_sw")),
              as.numeric(.recycle(sdf_sb, n, "sdf_sb")))
}

#' Annual basal area increment
#'
#' Basal area increment (m2/ha/yr) from the GYPSY basal area increment models.
#' [GYPSY_Sim()] adds these each year to grow BA forward.
#'
#' @param spp Character vector of species: `"AW"`, `"SB"`, `"SW"` or `"PL"`.
#' @param bhage Breast height age.
#' @param sc Species composition (0-1), based on density of species at or
#'   above breast height.
#' @param n0 Initial density (`*_den0` from [GYPSY_SDF()]).
#' @param si_b Site index at breast height age 50.
#' @param psi Percent stocking index (0 when unknown).
#' @param ba Current basal area of the species (m2/ha).
#' @param sdf_aw,sdf_sb,sdf_sw,sdf_pl Competitor SDFs used by SW and PL; 0
#'   for species that have not reached breast height.
#'
#' @return BA increment (m2/ha/yr). Length-1 arguments are recycled to the
#'   longest input.
#' @seealso [GYPSY_Sim()]
#' @export
#' @examples
#' GYPSY_BAinc("PL", bhage = 20, sc = 1, n0 = 2500, si_b = 17, ba = 15)
GYPSY_BAinc <- function(spp,
                        bhage,
                        sc,
                        n0,
                        si_b,
                        psi    = 0,
                        ba     = 0,
                        sdf_aw = 0,
                        sdf_sb = 0,
                        sdf_sw = 0,
                        sdf_pl = 0){
  n <- max(lengths(list(spp, bhage, sc, n0, si_b, psi, ba, sdf_aw, sdf_sb, sdf_sw, sdf_pl)))
  GYPSY_BAinc_cpp(toupper(as.character(.recycle(spp, n, "spp"))),
                  as.numeric(.recycle(bhage, n, "bhage")),
                  as.numeric(.recycle(sc, n, "sc")),
                  as.numeric(.recycle(n0, n, "n0")),
                  as.numeric(.recycle(si_b, n, "si_b")),
                  as.numeric(.recycle(psi, n, "psi")),
                  as.numeric(.recycle(ba, n, "ba")),
                  as.numeric(.recycle(sdf_aw, n, "sdf_aw")),
                  as.numeric(.recycle(sdf_sb, n, "sdf_sb")),
                  as.numeric(.recycle(sdf_sw, n, "sdf_sw")),
                  as.numeric(.recycle(sdf_pl, n, "sdf_pl")))
}

#' Total and merchantable volume
#'
#' Total volume, merchantable volume and merchantable density for one
#' utilization standard, from the GYPSY volume models.
#'
#' @param spp Character vector of species: `"AW"`, `"SB"`, `"SW"` or `"PL"`.
#' @param ba Basal area (m2/ha).
#' @param topht Top height (m). Volume needs top height over 1.3 m.
#' @param n Density (stems/ha).
#' @param bhage Breast height age.
#' @param stump_ht Stump height (m), at most 0.3.
#' @param stump_dob Minimum stump diameter outside bark (cm).
#' @param top_dib Minimum top diameter inside bark (cm).
#'
#' @return A list with `tvol` (total volume, m3/ha), `mvol` (merchantable
#'   volume, m3/ha, rounded to 1e-5 as in the SAS), `mn` (merchantable
#'   stems/ha) and `qmd` (quadratic mean diameter, cm). Length-1 arguments are
#'   recycled to the longest input.
#' @seealso [GYPSY_Sim()]
#' @export
#' @examples
#' GYPSY_Vol(c("SW", "PL"), ba = c(20, 30), topht = c(18, 20), n = c(800, 1500),
#'           bhage = c(55, 60))
GYPSY_Vol <- function(spp,
                      ba,
                      topht,
                      n,
                      bhage,
                      stump_ht  = 0.3, # m, <= 0.3
                      stump_dob = 15,  # cm
                      top_dib   = 10){ # cm
  len <- max(lengths(list(spp, ba, topht, n, bhage, stump_ht, stump_dob, top_dib)))
  if(any(stump_ht > 0.3)) stop("stump_ht must be <= 0.3")
  GYPSY_Vol_cpp(toupper(as.character(.recycle(spp, len, "spp"))),
                as.numeric(.recycle(ba, len, "ba")),
                as.numeric(.recycle(topht, len, "topht")),
                as.numeric(.recycle(n, len, "n")),
                as.numeric(.recycle(bhage, len, "bhage")),
                as.numeric(.recycle(stump_ht, len, "stump_ht")),
                as.numeric(.recycle(stump_dob, len, "stump_dob")),
                as.numeric(.recycle(top_dib, len, "top_dib")))
}

# ---- Full simulation ----------------------------------------------------------

#' Project stands with GYPSY
#'
#' Runs the full GYPSY growth and yield projection, following the 2009 GYPSY
#' SAS program (`GYPSY_run.sas`). Input is one row per stand and species;
#' output is one row per stand and year with species in columns.
#'
#' Each row needs two of `si`, `age` and `ht`. As in the SAS:
#' * SI is computed from `ht` and `age` (total age at least 7.5 for AW and PL,
#'   9.5 for SB and SW);
#' * a supplied `si` overrides that SI;
#' * when `si` and `ht` are both given, total age is recomputed from them and
#'   the supplied `age` is not used.
#'
#' Species without an SI are not grown, but their age still counts toward the
#' stand's leading (oldest) age. Stands with no usable age are kept and
#' project as `NA`, with a warning.
#'
#' The projection starts at the fractional part of the leading species' total
#' age and runs `simlen` years. Observed basal area, when given, is used to
#' localize the predicted basal area.
#'
#' @param stand Stand identifier.
#' @param stand_age Stand age (total), used for the `standage` output column.
#' @param spp Species: `"AW"`, `"SB"`, `"SW"` or `"PL"`.
#' @param n Density (stems/ha).
#' @param si Site index (m), of the type given by `sitype`.
#' @param age Age, of the type given by `agetype`.
#' @param ht Top height (m).
#' @param ps Percent stocking (`NA` if unknown; non-numeric values such as
#'   `"."` are treated as missing).
#' @param ba Basal area (m2/ha, `NA` if unknown).
#' @param sitype `"B"` SI at breast height age 50, or `"T"` at total age 50.
#' @param agetype `"T"` total, `"B"` breast height or `"S"` stump height age.
#' @param simlen Number of years to project.
#' @param psflag Set to 0 to ignore `ps` (turns off the percent stocking
#'   models).
#' @param baflag Set to 0 to ignore `ba` (turns off BA localization).
#' @param bh,sh Breast and stump height (m) for the SI and age equations.
#' @param stump_ht,stump_dob,top_dib Merchantability standard: stump height
#'   (m, at most 0.3), minimum stump diameter outside bark (cm) and minimum top
#'   diameter inside bark (cm). Either one value for all species or a vector
#'   named `AW`, `SB`, `SW`, `PL`.
#' @param iter Maximum iterations for the iterative solvers.
#' @param tol Convergence tolerance for the iterative solvers.
#' @param detail If `TRUE`, also return `sc_b_*` (breast height based species
#'   composition), `ba1_*` (BA grown from the start) and `ba2_*` (BA grown
#'   from the observed BA).
#'
#' @return A data.table with one row per stand and year:
#' * `stand`, `year` (1 to `simlen`), `totage` (leading species total age),
#'   `standage`;
#' * for each species (`_aw`, `_sb`, `_sw`, `_pl`): `totage_`, `bhage_`,
#'   `topht_`, `n_` (density), `sc_` (species composition), `ba_` (basal
#'   area), `qmd_`, `tvol_` (total volume), `mvol_` (merchantable volume),
#'   `mn_` (merchantable stems), `mai_t_` and `mai_m_` (total and
#'   merchantable MAI by species age);
#' * totals: `tvol_con`, `mvol_con` (conifers), `tvol_all`, `mvol_all`, and
#'   MAI by leading age `mai_t_all`, `mai_m_all`, `mai_t_con`, `mai_m_con`,
#'   `mai_t_dec`, `mai_m_dec`.
#'
#' Attribute `"stand_init"` holds one row per stand with the widened inputs,
#' SI, ages, PSI, SDF, initial densities (`*_den0`), species composition at
#' breast height age 0 (`*_sc0`), starting BA (`*_ba0`) and the observed /
#' predicted BA ratio (`*_rate`).
#' @seealso The components: [GYPSY_SI()], [GYPSY_Age()], [GYPSY_PSI()],
#'   [GYPSY_SDF()], [GYPSY_N()], [GYPSY_BAinc()], [GYPSY_Vol()].
#' @export
#' @examples
#' out <- GYPSY_Sim(stand = c(1, 1), stand_age = 20, spp = c("AW", "PL"),
#'                  n = c(375, 1825), age = c(6.2, 17.5), ht = c(3.5, 5.05),
#'                  ba = c(0.13, 3.87))
#' out[year %in% c(10, 50, 100), .(year, totage, ba_pl, tvol_pl, mvol_all)]
GYPSY_Sim <- function(stand,  # Outputs: one row per stand and year, species in columns (as SAS examp19)
                      stand_age,
                      spp,
                      n,               # Density (stems/ha)
                      si      = NA_real_,
                      age     = NA_real_,
                      ht      = NA_real_,
                      ps      = NA_real_,
                      ba      = NA_real_,
                      sitype  = "B",   # T for Total or B for Breast Height
                      agetype = "T",   # T = Total, B = BH, S = SH
                      simlen  = 200,   # Simulation length
                      psflag  = 1,     # Percent stocking override, set to 0 turn off PS model
                      baflag  = 1,     # BA localization override, set to 0 turn off BA localization
                      bh      = 1.3,   # Optional breast height specification
                      sh      = 0.3,   # Optional stump height specification (years to stump height)
                      stump_ht  = 0.3, # Merchantability standards: one value for all species,
                      stump_dob = 15,  #   or vectors named AW, SB, SW, PL
                      top_dib   = 10,
                      iter    = 1000,  # Maximum iterations for hill-climbing algorithms
                      tol     = 1e-8,  # Minimum difference between iterations before extracting results from hill-climbing algorithms
                      detail  = FALSE){ # Also return BA1 / BA2 and BH age based species composition

  # ---- Input checks -------------------------------------------------------------
  input_length <- length(stand)
  stand_age <- .recycle(stand_age, input_length, "stand_age")
  spp <- toupper(.recycle(spp, input_length, "spp"))
  n   <- .recycle(n,   input_length, "n")
  si  <- .recycle(si,  input_length, "si")
  age <- .recycle(age, input_length, "age")
  ht  <- .recycle(ht,  input_length, "ht")
  ps  <- .recycle(ps,  input_length, "ps")
  ba  <- .recycle(ba,  input_length, "ba")

  if(!sitype %in% c("B", "T")) stop("sitype must be B or T")
  if(!agetype %in% c("T", "B", "S")) stop("agetype must be T, B, or S")

  bad_spp <- which(!spp %in% .gypsy_spp)
  if(length(bad_spp) > 0)
    stop("Species must be AW, SB, SW, or PL. Check row(s): ", paste(bad_spp, collapse = ", "))

  miss_input <- which((is.na(si) + is.na(age) + is.na(ht)) > 1) # Check for missing input data
  if(length(miss_input) > 0)
    stop("Simulation requires 2 of SI, Age, and Top Height. Check row(s): ", paste(miss_input, collapse = ", "))

  merch <- list(stump_ht  = .per_spp(stump_ht,  "stump_ht"),
                stump_dob = .per_spp(stump_dob, "stump_dob"),
                top_dib   = .per_spp(top_dib,   "top_dib"))
  if(any(merch$stump_ht > 0.3)) stop("stump_ht must be <= 0.3")

  dt1 <- data.table("stand"     = stand,     # Combine inputs into single table
                    "stand_age" = as.numeric(stand_age),
                    "spp"       = spp,
                    "n"         = as.numeric(n),
                    "si"        = as.numeric(si),
                    "age"       = as.numeric(age),
                    "ht"        = as.numeric(ht),
                    "ps"        = suppressWarnings(as.numeric(ps)), # "." from SAS-style inputs becomes NA
                    "ba"        = as.numeric(ba))

  dup <- dt1[, .I[.N > 1], by = .(stand, spp)]$V1 # Find duplicate stand species groups
  if(length(dup) > 0)
    stop("Duplicate stand species groups in row(s): ", paste(dup, collapse = ", "))

  if(psflag == 0) dt1[, ps := NA_real_]
  if(baflag == 0) dt1[, ba := NA_real_]

  # ---- SI and age (SAS data v1) ---------------------------------------------------
  # Same order as the SAS: SI from height and age, then a supplied SI overrides it,
  # then total age is recomputed from SI and height when both are given.
  dt1[, `:=`(si_t = NA_real_, si_bh = NA_real_, y2bh = NA_real_)]

  min_age <- ifelse(dt1$spp %in% c("AW", "PL"), 7.5, 9.5)
  i1 <- dt1[, which(!is.na(age) & !is.na(ht) & ht > 0 & age >= min_age)]
  if(length(i1) > 0){
    dt1[i1, `:=`(si_t  = GYPSY_SI(spp, ht, age, agetype, ndx = 1, iter, tol, bh, sh),
                 si_bh = GYPSY_SI(spp, ht, age, agetype, ndx = 2, iter, tol, bh, sh),
                 y2bh  = GYPSY_SI(spp, ht, age, agetype, ndx = 3, iter, tol, bh, sh))]
  }

  i2 <- dt1[, which(!is.na(si) & si > 0)]
  if(length(i2) > 0){
    if(sitype == "B"){
      dt1[i2, `:=`(si_bh = si,
                   si_t  = GYPSY_SI_convert(spp, si, sitype = "B", ndx = 4, iter, tol, bh, sh))]
    } else {
      dt1[i2, `:=`(si_t  = si,
                   si_bh = GYPSY_SI_convert(spp, si, sitype = "T", ndx = 4, iter, tol, bh, sh))]
    }
    dt1[i2, y2bh := GYPSY_SI_convert(spp, si_t, sitype = "T", ndx = 5, iter, tol, bh, sh)]
  }

  # Total age from the supplied age
  if(agetype == "T"){
    dt1[, totage := age]
  } else if(agetype == "B"){
    dt1[, totage := age + y2bh]
  } else {
    dt1[, totage := NA_real_]
    i_s <- dt1[, which(!is.na(si_t))]
    if(length(i_s) > 0)
      dt1[i_s, totage := age + GYPSY_SI_convert(spp, si_t, sitype = "T", ndx = 6, iter, tol, bh, sh)]
  }

  i3 <- dt1[, which(!is.na(si) & si > 0 & !is.na(ht) & ht > 0)]
  if(length(i3) > 0)
    dt1[i3, totage := GYPSY_Age(spp, si_t, ht, sitype = "T", ndx = 1, iter, tol, bh, sh)]

  dt1[, bhage := NA_real_]
  dt1[!is.na(y2bh) & y2bh > 0, bhage := totage - y2bh]

  no_si <- dt1[, which(is.na(si_t))]
  if(length(no_si) > 0)
    warning("No SI for row(s) ", paste(no_si, collapse = ", "),
            " (age below 7.5 for AW/PL or 9.5 for SB/SW); these species are not grown but still count toward the leading age.")

  # Leading (maximum) total age per stand
  dt1[, totage_max := suppressWarnings(max(totage, na.rm = TRUE)), by = stand]
  dt1[!is.finite(totage_max), totage_max := NA_real_]
  no_age <- dt1[is.na(totage_max), unique(stand)]
  if(length(no_age) > 0) # As in the SAS, these stands are kept and project as all NA
    warning("No total age could be determined for stand(s): ", paste(no_age, collapse = ", "))

  # ---- Pivot species wider (SAS examp3) -----------------------------------------
  dt2 <- dt1[, .(stand_age = stand_age[1], totage_max = totage_max[1]), by = stand]
  setorder(dt2, stand) # SAS processes stands in sorted order
  src_col <- c("si_t", "si_bh", "totage", "bhage", "y2bh", "n", "ps", "ba", "ht")
  for(s in .gypsy_spp){
    add_col <- paste0(tolower(s), "_", src_col)
    dt2[, (add_col) := NA_real_]
    dt2[dt1[spp == s], on = "stand", (add_col) := mget(paste0("i.", src_col))]
  }

  # ---- Stand-level parameters (PSI, SDF, initial density) ------------------------
  dt2[, c("aw_psi", "sb_psi", "sw_psi", "pl_psi") :=
        GYPSY_PSI(stand,
                  aw_ps, aw_totage, aw_si_bh,
                  sb_ps, sb_totage, sb_si_bh,
                  sw_ps, sw_totage, sw_si_bh,
                  pl_ps, pl_totage, pl_si_bh)]

  sdf <- GYPSY_SDF(dt2$stand,
                   dt2$aw_bhage,  dt2$aw_n, dt2$aw_si_bh, dt2$aw_psi,
                   dt2$sb_totage, dt2$sb_n, dt2$sb_si_bh, dt2$sb_psi,
                   dt2$sw_totage, dt2$sw_n, dt2$sw_si_bh, dt2$sw_psi,
                   dt2$pl_totage, dt2$pl_n, dt2$pl_si_bh, dt2$pl_psi,
                   iter, tol)
  dt2[, names(sdf) := sdf]

  # ---- Grow forward (SAS examp6 - examp19) ----------------------------------------
  res <- GYPSY_Project_cpp(as.list(dt2),
                           merch$stump_ht, merch$stump_dob, merch$top_dib,
                           as.integer(simlen), 180L, detail)

  out <- as.data.table(res$rows)
  out[, stand := dt2$stand[stand_idx]]
  out[, stand_idx := NULL]
  setcolorder(out, "stand")

  init <- cbind(dt2, as.data.table(res$init))
  setattr(out, "stand_init", init) # Stand-level inputs and parameters (SI, PSI, SDF, SC0, BA0, ...)
  out[]
}
