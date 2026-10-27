# Project stand with GYPSY
# Follows GYPSY_run.sas (MAY 21, 2009 functions)
# Variables:
#   SI, Age, TopHt, BA, N

.gypsy_spp <- c("AW", "SB", "SW", "PL")

# Repeat single-value inputs to length n, otherwise require length n
.recycle <- function(x, n, name){
  if(length(x) == 1) return(rep(x, n))
  if(length(x) != n) stop(name, " must be a single value or length ", n)
  x
}

# Recycle a named list of arguments to length n. Gives one warning naming the
# arguments in `warn` (variable inputs, not settings) that the caller supplied
# as a single value and were repeated; defaults left in place are not warned about.
.recycle_args <- function(args, n = max(lengths(args)), warn = names(args),
                          env = parent.frame()){
  supplied <- !vapply(warn, function(a) eval(call("missing", as.name(a)), env), logical(1))
  rec <- warn[supplied & lengths(args[warn]) == 1]
  if(n > 1 && length(rec) > 0)
    warning("Single value(s) repeated to match ", n, " rows: ",
            paste(rec, collapse = ", "), call. = FALSE)
  Map(.recycle, args, n, names(args))
}

# Per-species value from a scalar or a vector named by species
.per_spp <- function(x, name){
  if(length(x) == 1) return(rep(as.numeric(x), 4))
  if(is.null(names(x)) || !all(.gypsy_spp %in% toupper(names(x))))
    stop(name, " must be a single value or named with all of AW, SB, SW, PL")
  names(x) <- toupper(names(x))
  as.numeric(x[.gypsy_spp])
}

# ---- Stand-level components -------------------------------------------------

#' Percent stocking index
#'
#' Percent stocking at total age 50 for each species in a stand, from the
#' GYPSY percent stocking models. Requires one row per stand with species in
#' separate columns.
#'
#' @param stand Stand identifiers (used for length only).
#' @param aw_ps,sb_ps,sw_ps,pl_ps Observed percent stocking (NA if unknown).
#' @param aw_age_t,sb_age_t,sw_age_t,pl_age_t Total age.
#' @param aw_si_b,sb_si_b,sw_si_b,pl_si_b Site index at breast height age 50.
#'
#' @return A list with `aw_psi`, `sb_psi`, `sw_psi` and `pl_psi`, capped at
#'   100. As in the SAS, AW, SB and SW are 0 when they can't be computed, PL is
#'   `NA`. Arguments given as a single value are repeated for every stand,
#'   with a warning.
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
  a <- lapply(.recycle_args(list(aw_ps = aw_ps, aw_age_t = aw_age_t, aw_si_b = aw_si_b,
                                 sb_ps = sb_ps, sb_age_t = sb_age_t, sb_si_b = sb_si_b,
                                 sw_ps = sw_ps, sw_age_t = sw_age_t, sw_si_b = sw_si_b,
                                 pl_ps = pl_ps, pl_age_t = pl_age_t, pl_si_b = pl_si_b),
                            n = length(stand)),
              as.numeric)
  GYPSY_PSI_cpp(as.character(stand),
                a$aw_ps, a$aw_age_t, a$aw_si_b,
                a$sb_ps, a$sb_age_t, a$sb_si_b,
                a$sw_ps, a$sw_age_t, a$sw_si_b,
                a$pl_ps, a$pl_age_t, a$pl_si_b)
}

#' Stand density factor
#'
#' Solves the GYPSY spatial density functions for each species' stand density
#' factor (SDF, density at age 50) from its observed density, and returns the
#' initial density implied by that SDF. Requires one row per stand with species
#' in separate columns.
#'
#' AW is computed only when its density is over 20 stems/ha, other species
#' when density is over 0.
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
#'   age 0). Arguments given as a single value are repeated for every stand,
#'   with a warning.
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
  a <- lapply(.recycle_args(list(aw_age_b = aw_age_b, aw_n = aw_n, aw_si_b = aw_si_b, aw_psi = aw_psi,
                                 sb_age_t = sb_age_t, sb_n = sb_n, sb_si_b = sb_si_b, sb_psi = sb_psi,
                                 sw_age_t = sw_age_t, sw_n = sw_n, sw_si_b = sw_si_b, sw_psi = sw_psi,
                                 pl_age_t = pl_age_t, pl_n = pl_n, pl_si_b = pl_si_b, pl_psi = pl_psi),
                            n = length(stand)),
              as.numeric)
  GYPSY_SDF_cpp(as.character(stand),
                a$aw_age_b, a$aw_n, a$aw_si_b, a$aw_psi,
                a$sb_age_t, a$sb_n, a$sb_si_b, a$sb_psi,
                a$sw_age_t, a$sw_n, a$sw_si_b, a$sw_psi,
                a$pl_age_t, a$pl_n, a$pl_si_b, a$pl_psi,
                as.integer(iter), as.numeric(tol))
}

# ---- Species-level components -------------------------------------------------

#' Top height from site index and age
#'
#' @param spp Character vector of species: `"AW"`, `"SB"`, `"SW"` or `"PL"`.
#' @param si_t Site index at total age 50 (m).
#' @param age_t Total age.
#'
#' @return Top height (m). Arguments given as a single value are repeated to
#'   match the longest input, with a warning.
#' @seealso [GYPSY_SI()] for the inverse.
#' @export
#' @examples
#' GYPSY_Ht(c("AW", "SB", "SW", "PL"), si_t = rep(15, 4), age_t = rep(40, 4))
GYPSY_Ht <- function(spp,
                     si_t,
                     age_t){
  a <- .recycle_args(list(spp = spp, si_t = si_t, age_t = age_t))
  GYPSY_Ht_cpp(toupper(as.character(a$spp)),
               as.numeric(a$si_t),
               as.numeric(a$age_t))
}

#' Density at age
#'
#' Density (stems/ha) at a given age from a species' stand density factor.
#'
#' @param spp Species group: `"AW"`, `"SB"`, `"SW"` or `"PL"`.
#' @param age Breast height age for AW, total age for other species.
#' @param sdf Stand density factor from [GYPSY_SDF()].
#' @param si_b Site index at breast height age 50.
#' @param psi Percent stocking index (0 when unknown).
#' @param sdf_aw,sdf_sw,sdf_sb Other species' SDFs in the stand (SW uses AW;
#'   PL uses AW, SW and SB); 0 when absent.
#'
#' @return Density (stems/ha); `NA` when `sdf` is not positive. Arguments
#'   given as a single value are repeated to match the longest input, with a
#'   warning for those supplied rather than left at their defaults.
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
  a <- .recycle_args(list(spp = spp, age = age, sdf = sdf, si_b = si_b, psi = psi,
                          sdf_aw = sdf_aw, sdf_sw = sdf_sw, sdf_sb = sdf_sb))
  GYPSY_N_cpp(toupper(as.character(a$spp)),
              as.numeric(a$age),
              as.numeric(a$sdf),
              as.numeric(a$si_b),
              as.numeric(a$psi),
              as.numeric(a$sdf_aw),
              as.numeric(a$sdf_sw),
              as.numeric(a$sdf_sb))
}

#' Annual basal area increment
#'
#' Basal area increment (m2/ha/yr) from the GYPSY basal area increment models.
#' [GYPSY_Sim()] adds these each year to grow BA forward.
#'
#' @param spp Species group: `"AW"`, `"SB"`, `"SW"` or `"PL"`.
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
#' @return BA increment (m2/ha/yr). Arguments given as a single value are
#'   repeated to match the longest input, with a warning for those supplied
#'   rather than left at their defaults.
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
  a <- .recycle_args(list(spp = spp, bhage = bhage, sc = sc, n0 = n0, si_b = si_b,
                          psi = psi, ba = ba, sdf_aw = sdf_aw, sdf_sb = sdf_sb,
                          sdf_sw = sdf_sw, sdf_pl = sdf_pl))
  GYPSY_BAinc_cpp(toupper(as.character(a$spp)),
                  as.numeric(a$bhage),
                  as.numeric(a$sc),
                  as.numeric(a$n0),
                  as.numeric(a$si_b),
                  as.numeric(a$psi),
                  as.numeric(a$ba),
                  as.numeric(a$sdf_aw),
                  as.numeric(a$sdf_sb),
                  as.numeric(a$sdf_sw),
                  as.numeric(a$sdf_pl))
}

#' Total and merchantable volume
#'
#' Total volume, merchantable volume and merchantable density for one
#' utilization standard.
#'
#' @param spp Species group: `"AW"`, `"SB"`, `"SW"` or `"PL"`.
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
#'   stems/ha) and `qmd` (quadratic mean diameter, cm). Arguments given as a
#'   single value are repeated to match the longest input, with a warning for
#'   the stand inputs (`spp`, `ba`, `topht`, `n`, `bhage`) but not the
#'   merchantability standard.
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
  if(any(stump_ht > 0.3)) stop("stump_ht must be <= 0.3")
  a <- .recycle_args(list(spp = spp, ba = ba, topht = topht, n = n, bhage = bhage,
                          stump_ht = stump_ht, stump_dob = stump_dob, top_dib = top_dib),
                     warn = c("spp", "ba", "topht", "n", "bhage")) # Merch standards recycle silently
  GYPSY_Vol_cpp(toupper(as.character(a$spp)),
                as.numeric(a$ba),
                as.numeric(a$topht),
                as.numeric(a$n),
                as.numeric(a$bhage),
                as.numeric(a$stump_ht),
                as.numeric(a$stump_dob),
                as.numeric(a$top_dib))
}

# ---- Full simulation ----------------------------------------------------------

#' Project stands with GYPSY
#'
#' Runs the full GYPSY growth and yield projection. Input is one row per stand
#' and species; output is one row per stand and year with species in columns.
#'
#' Each row needs two of `si`, `age` and `ht`:
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
#' Observed percent stocking (`ps`), when given, is converted to a percent
#' stocking index (PSI, percent stocking at total age 50) with the percent
#' stocking models. The PSI then adjusts the density and basal area increment
#' models; without it, those models use their no-stocking form. Set
#' `psflag = 0` to ignore `ps`.
#'
#' `stand_age`, `spp`, `n`, `si`, `age`, `ht`, `ps` and `ba` must be either a
#' single value or the same length as `stand`. A single value is repeated for
#' every row, with a warning for those supplied rather than left at their
#' defaults.
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
#' @param data Optional input table (data.frame or data.table) with one row
#'   per stand and species. Columns named `stand`, `stand_age`, `spp`, `n`,
#'   `si`, `age`, `ht`, `ps` and `ba` are used for those arguments; other
#'   columns are ignored. Arguments given directly take precedence over
#'   columns. A missing `stand`, `stand_age`, `spp` or `n` column is an
#'   error; other missing columns are treated as missing values, with a
#'   warning.
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
#' out <- GYPSY_Sim(stand = c(1, 1), stand_age = c(20, 20), spp = c("AW", "PL"),
#'                  n = c(375, 1825), age = c(6.2, 17.5), ht = c(3.5, 5.05),
#'                  ba = c(0.13, 3.87))
#' out[year %in% c(10, 50, 100), .(year, totage, ba_pl, tvol_pl, mvol_all)]
#'
#' # From an input table
#' inp <- data.table::fread(system.file("extdata", "example_input.csv", package = "rGYPSY"))
#' out <- GYPSY_Sim(data = inp, simlen = 150)
GYPSY_Sim <- function(stand,
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
                      detail  = FALSE, # Also return BA1 / BA2 and BH age based species composition
                      data    = NULL){ # Optional input table with columns named as the arguments

  # ---- Columns from data --------------------------------------------------------
  if(!missing(stand) && is.data.frame(stand))
    stop("To pass an input table, use GYPSY_Sim(data = ...)")
  if(!is.null(data)){
    if(!is.data.frame(data)) stop("data must be a data.frame or data.table")
    cols <- c("stand", "stand_age", "spp", "n", "si", "age", "ht", "ps", "ba")
    env  <- environment()
    from_data <- cols[vapply(cols, function(a) eval(call("missing", as.name(a)), env), logical(1))]
    not_found <- setdiff(from_data, names(data))
    req_not_found <- intersect(not_found, c("stand", "stand_age", "spp", "n"))
    if(length(req_not_found) > 0)
      stop("Required column(s) not found in data: ", paste(req_not_found, collapse = ", "))
    if(length(not_found) > 0)
      warning("Column(s) not found in data, treated as missing: ",
              paste(not_found, collapse = ", "), call. = FALSE)
    for(a in intersect(from_data, names(data))) assign(a, data[[a]]) # Arguments given directly take precedence
  }

  # ---- Inputs to one value per row ----------------------------------------------
  input_length <- length(stand)
  a <- .recycle_args(list(stand_age = stand_age, spp = spp, n = n, si = si,
                          age = age, ht = ht, ps = ps, ba = ba),
                     n = input_length)

  merch <- list(stump_ht  = .per_spp(stump_ht,  "stump_ht"),
                stump_dob = .per_spp(stump_dob, "stump_dob"),
                top_dib   = .per_spp(top_dib,   "top_dib"))

  # Stands are passed to C++ as codes in sorted order (radix sort orders
  # character IDs in the C locale, as data.table does)
  ids  <- sort(unique(stand), method = "radix")
  code <- match(stand, ids) - 1L

  # ---- Run the model (src/core/gypsy_sim.h) ---------------------------------------
  res <- GYPSY_Sim_cpp(code, as.character(ids),
                       stand_age = as.numeric(a$stand_age),
                       spp       = toupper(as.character(a$spp)),
                       n         = as.numeric(a$n),
                       si        = as.numeric(a$si),
                       age       = as.numeric(a$age),
                       ht        = as.numeric(a$ht),
                       ps        = suppressWarnings(as.numeric(a$ps)), # "." from SAS-style inputs becomes NA
                       ba        = as.numeric(a$ba),
                       sitype    = as.character(sitype),
                       agetype   = as.character(agetype),
                       simlen    = as.integer(simlen),
                       psflag    = psflag != 0,
                       baflag    = baflag != 0,
                       bh        = as.numeric(bh),
                       sh        = as.numeric(sh),
                       stump_ht  = merch$stump_ht,
                       stump_dob = merch$stump_dob,
                       top_dib   = merch$top_dib,
                       iter      = as.integer(iter),
                       tol       = as.numeric(tol),
                       detail    = isTRUE(detail))
  for(w in res$warnings) warning(w, call. = FALSE)

  out <- as.data.table(res$rows)
  out[, stand := ids[stand_idx]]
  out[, stand_idx := NULL]
  setcolorder(out, "stand")

  init <- as.data.table(c(list(stand = ids), res$init))
  setattr(out, "stand_init", init) # Stand-level inputs and parameters (SI, PSI, SDF, SC0, BA0, ...)
  out[]
}
