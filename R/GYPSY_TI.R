# Tree improvement: site index adjusted for genetic worth
# Luo and Thomas (2021), section 2.4

#' Site index adjusted for tree improvement
#'
#' Raises a total age site index for the genetic worth of an improved
#' seedlot, using the age-age correlation between the selection age and the
#' site index reference age (total age 50), following section 2.4 of Luo and Thomas (2021)
#' \doi{10.1016/j.foreco.2020.118865}:
#'
#' \deqn{r = c + m \ln(sel\_age / 50)}{r = c + m ln(sel_age / 50)}
#' \deqn{SI_{t,improved} = ((gw / 100) r + 1) SI_t}{SI_t improved = ((gw / 100) r + 1) SI_t}
#'
#' with the adjusted Lambeth parameters c = 0.993, m = 0.255 for SW and
#' c = 1.017, m = 0.445 for PL (their Table 2). Only SW and PL have
#' parameters; AW and SB rows must have no genetic worth (`NA` or 0).
#'
#' The adjustment is always to the total age SI; a breast height age SI is
#' converted to it first (see [GYPSY_SI_convert()]). The correlation is not
#' capped, so a selection age close to or above 50 can give r > 1.
#'
#' @param spp Species group: `"AW"`, `"SB"`, `"SW"` or `"PL"`.
#' @param si_t Site index at total age 50 (m).
#' @param gw Genetic worth for height at the selection age (%). 0 leaves the
#'   SI unchanged. A missing (`NA`) SW or PL value is treated as 0, with a
#'   warning; AW and SB have no genetic worth, so `NA` is expected there.
#' @param sel_age Selection age (years) at which `gw` was estimated.
#'
#' @return Improved site index at total age 50 (m). Arguments given as a
#'   single value are repeated to match the longest input, with a warning.
#'   Use [GYPSY_SI_convert()] for the breast height age SI.
#' @references Luo, D. and Thomas, B.R. 2021. An analysis of age-age
#'   correlations in white spruce and lodgepole pine and how it applies to the
#'   growth and yield projection system (GYPSY) in Alberta. Forest Ecology
#'   and Management 482: 118865. \doi{10.1016/j.foreco.2020.118865}
#' @seealso [GYPSY_Sim_TI()] to project improved stands.
#' @export
#' @examples
#' # 10% genetic worth selected at age 15, on a site index of 20 m
#' GYPSY_SI_TI(c("SW", "PL"), si_t = c(20, 20), gw = c(10, 10), sel_age = c(15, 15))
GYPSY_SI_TI <- function(spp,
                        si_t,
                        gw,
                        sel_age){
  a <- .recycle_args(list(spp = spp, si_t = si_t, gw = gw, sel_age = sel_age))
  GYPSY_SI_TI_cpp(toupper(as.character(a$spp)),
                  as.numeric(a$si_t),
                  as.numeric(a$gw),
                  as.numeric(a$sel_age))
}

#' Project stands of improved seedlots with GYPSY
#'
#' Runs [GYPSY_Sim()] with each row's site index adjusted for tree
#' improvement. The stand is initialized from the same inputs as
#' [GYPSY_Sim()]; each SW and PL row with a genetic worth then has its total
#' age SI raised with [GYPSY_SI_TI()] before the projection, so top height,
#' density and basal area increment all follow the improved SI. The method
#' is from section 2.4 of Luo and Thomas (2021)
#' \doi{10.1016/j.foreco.2020.118865}.
#'
#' Total age is kept from the inputs. The breast height age SI, years to
#' breast height and breast height age are recomputed from the improved SI.
#'
#' @inheritParams GYPSY_Sim
#' @param gw Genetic worth for height at the selection age (%), per row. 0
#'   leaves the row unadjusted. Only SW and PL rows can have a genetic worth;
#'   a missing (`NA`) SW or PL value, including when `gw` is not given, is
#'   treated as 0, with a warning.
#' @param sel_age Selection age (years) at which `gw` was estimated, per row.
#' @param data Optional input table, as for [GYPSY_Sim()]. Columns `gw` and
#'   `sel_age` are also used.
#'
#' @return As for [GYPSY_Sim()]. The `"stand_init"` attribute's `*_si_t` and
#'   `*_si_bh` are the improved SIs, and it also has, by species, the SIs
#'   before adjustment (`*_si_t_base`, `*_si_bh_base`), `*_gw`, `*_sel_age`
#'   and the age-age correlation used (`*_aacor`).
#' @inherit GYPSY_SI_TI references
#' @seealso [GYPSY_SI_TI()]
#' @export
#' @examples
#' # A PL stand at total age 15, unimproved and with 10% genetic worth
#' # selected at age 15
#' base <- GYPSY_Sim(1, 15, "PL", n = 2000, si = 18, age = 15, simlen = 100)
#' imp  <- GYPSY_Sim_TI(1, 15, "PL", n = 2000, si = 18, age = 15, simlen = 100,
#'                      gw = 10, sel_age = 15)
#' attr(imp, "stand_init")[, .(pl_si_bh_base, pl_si_bh, pl_aacor)]
#' cbind(base[year %in% c(51, 81), .(totage, mvol_unimproved = mvol_pl)],
#'       imp[year %in% c(51, 81), .(mvol_improved = mvol_pl)])
GYPSY_Sim_TI <- function(stand,
                         stand_age,
                         spp,
                         n,
                         si      = NA_real_,
                         age     = NA_real_,
                         ht      = NA_real_,
                         ps      = NA_real_,
                         ba      = NA_real_,
                         gw      = NA_real_,
                         sel_age = NA_real_,
                         sitype  = "B",
                         agetype = "T",
                         simlen  = 200,
                         psflag  = 1,
                         baflag  = 1,
                         bh      = 1.3,
                         sh      = 0.3,
                         stump_ht  = 0.3,
                         stump_dob = 15,
                         top_dib   = 10,
                         iter    = 1000,
                         tol     = 1e-8,
                         detail  = FALSE,
                         data    = NULL){
  .run_sim(environment(), "GYPSY_Sim_TI", ti = TRUE)
}
