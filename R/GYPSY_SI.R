# GYPSY_SI (Calculate Site Index from GYPSY)
# Calculate site index in m from GYPSY model

#' Site index from top height and age
#'
#' Solves the GYPSY top height equations for site index (top height at age
#' 50) using fixed-point iteration given an observed top height and age.
#'
#' @param spp Species group: `"AW"`, `"SB"`, `"SW"` or `"PL"`.
#' @param ht Top height (m).
#' @param age Age in years, of the type given by `agetype`. Total age must be
#'   at least 7.5 for AW and PL and 9.5 for SB and SW.
#' @param agetype Type of `age`: `"T"` total, `"B"` breast height or `"S"`
#'   stump height age.
#' @param ndx Output: `1` SI at total age 50 (SI_t), `2` SI at breast height
#'   age 50 (SI_bh), `3` years to breast height, `4` years to stump height.
#' @param iter Maximum iterations.
#' @param tol Convergence tolerance.
#' @param bh Breast height (m).
#' @param sh Stump height (m).
#'
#' @return Numeric vector, one value per row. A warning is given for rows that
#'   don't converge within `iter` iterations. If `spp`, `ht` or `age`
#'   is given as a single value, it is repeated to match the longest input,
#'   with a warning.
#' @seealso [GYPSY_SI_convert()] to convert between SI_t and SI_bh,
#'   [GYPSY_Age()] for age from SI and height.
#' @export
#' @examples
#' spp <- c("AW", "SB", "SW", "PL")
#' GYPSY_SI(spp, ht = rep(12, 4), age = rep(40, 4), ndx = 1) # SI_t
#' GYPSY_SI(spp, ht = rep(12, 4), age = rep(40, 4))          # SI_bh
GYPSY_SI <- function(spp,
                     ht,
                     age,
                     agetype = "T",
                     ndx = 2,
                     iter = 1000,
                     tol = 1e-8,
                     bh = 1.3,
                     sh = 0.3){
  a <- .recycle_args(list(spp = spp, ht = ht, age = age))
  GYPSY_SI_cpp(
    spp     = as.character(a$spp),
    ht      = as.numeric(a$ht),
    age     = as.numeric(a$age),
    agetype = as.character(agetype),
    ndx     = as.integer(ndx),
    iter    = as.integer(iter),
    tol     = as.numeric(tol),
    bh      = as.numeric(bh),
    sh      = as.numeric(sh)
  )
}
