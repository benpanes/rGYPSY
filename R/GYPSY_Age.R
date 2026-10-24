#' Age from site index and top height
#'
#' Solves the GYPSY top height equations for the age at which a stand with a
#' given site index reaches a given top height.
#'
#' @param spp Character vector of species: `"AW"`, `"SB"`, `"SW"` or `"PL"`.
#' @param si Site index (m), of the type given by `sitype`.
#' @param ht Top height (m).
#' @param sitype Type of `si`: `"B"` SI at breast height age 50 (SI_bh), or
#'   `"T"` SI at total age 50 (SI_t). SI_bh is first converted to SI_t by
#'   iteration.
#' @param ndx Output: `1` total age, `2` stump height age, `3` breast height
#'   age.
#' @inheritParams GYPSY_SI
#'
#' @return Numeric vector of ages, one per row.
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
                      sitype = "B", # T for Total or B for Breast Height
                      ndx = 1, # 1 for T, 2 for S, 3 for B
                      iter = 1000,
                      tol = 1e-8,
                      bh = 1.3,
                      sh = 0.3){
  GYPSY_Age_cpp(
    spp     = as.character(spp),
    si      = as.numeric(si),
    h       = as.numeric(ht),
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
#' Converts SI at breast height age 50 (SI_bh) to SI at total age 50 (SI_t)
#' or the reverse, or returns years to breast or stump height for a site
#' index.
#'
#' @param spp Character vector of species: `"AW"`, `"SB"`, `"SW"` or `"PL"`.
#' @param si Site index (m), of the type given by `sitype`.
#' @param sitype Type of `si`: `"B"` SI_bh or `"T"` SI_t.
#' @param ndx Output: `4` the other SI type (SI_t from SI_bh, or SI_bh from
#'   SI_t), `5` years to breast height, `6` years to stump height.
#' @inheritParams GYPSY_SI
#'
#' @return Numeric vector, one value per row.
#' @seealso [GYPSY_SI()], [GYPSY_Age()]
#' @export
#' @examples
#' spp <- c("AW", "SB", "SW", "PL")
#' GYPSY_SI_convert(spp, si = rep(16, 4))                        # SI_t from SI_bh
#' GYPSY_SI_convert(spp, si = rep(16, 4), sitype = "T")          # SI_bh from SI_t
#' GYPSY_SI_convert(spp, si = rep(16, 4), sitype = "T", ndx = 5) # years to BH
GYPSY_SI_convert <- function(spp,
                             si,
                             sitype = "B", # Input SI switch, use "T" for total
                             ndx = 4, # 4 for alternate SI, 5 for Y2BH, 6 for Y2SH
                             iter = 1000,
                             tol = 1e-8,
                             bh = 1.3,
                             sh = 0.3){
  GYPSY_Age_cpp(
    spp     = as.character(spp),
    si      = as.numeric(si),
    sitype  = as.character(sitype),
    ndx     = as.integer(ndx),
    iter    = as.integer(iter),
    tol     = as.numeric(tol),
    bh      = as.numeric(bh),
    sh      = as.numeric(sh)
  )
}
