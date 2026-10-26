#' Age from site index and top height
#'
#' Solves the GYPSY top height equations for the age at which a stand with a
#' given site index reaches a given top height.
#'
#' @param spp Species group: `"AW"`, `"SB"`, `"SW"` or `"PL"`.
#' @param si Site index (m), of the type given by `sitype`.
#' @param ht Top height (m).
#' @param sitype Type of `si`: `"B"` SI at breast height age 50 (SI_bh), or
#'   `"T"` SI at total age 50 (SI_t). SI_bh is first converted to SI_t by
#'   iteration.
#' @param ndx Output: `1` total age, `2` stump height age, `3` breast height
#'   age.
#' @inheritParams GYPSY_SI
#'
#' @return Numeric vector of ages, one per row. If `spp`, `si` or `ht` is
#'   given as a single value, it is repeated to match the longest input, with
#'   a warning.
#' @seealso [GYPSY_SI()], [GYPSY_SI_convert()]
#' @export
#' @examples
#' spp   <- c("AW", "SB", "SW", "PL")
#' si_bh <- GYPSY_SI(spp, ht = rep(12, 4), age = rep(40, 4))
#' GYPSY_Age(spp, si_bh, ht = rep(12, 4))          # total age, 40
#' GYPSY_Age(spp, si_bh, ht = rep(12, 4), ndx = 3) # breast height age
GYPSY_Age <- function(spp,
                      si,
                      ht,
                      sitype = "B",
                      ndx = 1,
                      iter = 1000,
                      tol = 1e-8,
                      bh = 1.3,
                      sh = 0.3){
  a <- .recycle_args(list(spp = spp, si = si, ht = ht))
  GYPSY_Age_cpp(
    spp     = as.character(a$spp),
    si      = as.numeric(a$si),
    h       = as.numeric(a$ht),
    sitype  = as.character(sitype),
    ndx     = as.integer(ndx),
    iter    = as.integer(iter),
    tol     = as.numeric(tol),
    bh      = as.numeric(bh),
    sh      = as.numeric(sh)
  )
}

#' Convert between site index types
#'
#' Converts SI between total age at 50 (SI_t) or breast height age at 50
#' (SI_bh), or returns years to breast or stump height from a known site index.
#'
#' @param spp Species group: `"AW"`, `"SB"`, `"SW"` or `"PL"`.
#' @param si Site index (m), of the type given by `sitype`.
#' @param sitype Type of `si`: `"B"` SI_bh or `"T"` SI_t.
#' @param ndx Output: `4` the other SI type (SI_t from SI_bh, or SI_bh from
#'   SI_t), `5` years to breast height, `6` years to stump height.
#' @inheritParams GYPSY_SI
#'
#' @return Numeric vector, one value per row. If `spp` or `si` is given as a
#'   single value, it is repeated to match the longest input, with a warning.
#' @seealso [GYPSY_SI()], [GYPSY_Age()]
#' @export
#' @examples
#' spp <- c("AW", "SB", "SW", "PL")
#' GYPSY_SI_convert(spp, si = rep(16, 4))                        # SI_t from SI_bh
#' GYPSY_SI_convert(spp, si = rep(16, 4), sitype = "T")          # SI_bh from SI_t
#' GYPSY_SI_convert(spp, si = rep(16, 4), sitype = "T", ndx = 5) # years to BH
GYPSY_SI_convert <- function(spp,
                             si,
                             sitype = "B",
                             ndx = 4,
                             iter = 1000,
                             tol = 1e-8,
                             bh = 1.3,
                             sh = 0.3){
  a <- .recycle_args(list(spp = spp, si = si))
  GYPSY_Age_cpp(
    spp     = as.character(a$spp),
    si      = as.numeric(a$si),
    sitype  = as.character(sitype),
    ndx     = as.integer(ndx),
    iter    = as.integer(iter),
    tol     = as.numeric(tol),
    bh      = as.numeric(bh),
    sh      = as.numeric(sh)
  )
}
