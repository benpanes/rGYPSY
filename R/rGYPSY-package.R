#' rGYPSY: GYPSY growth and yield functions
#'
#' A combined R and Rcpp implementation of the Alberta Growth and Yield Projection
#' System (GYPSY 2009).
#'
#' Component functions: [GYPSY_SI()], [GYPSY_SI_convert()], [GYPSY_Age()],
#' [GYPSY_Ht()], [GYPSY_PSI()], [GYPSY_SDF()], [GYPSY_N()], [GYPSY_BAinc()],
#' [GYPSY_Vol()]. Full projection: [GYPSY_Sim()].
#'
#' @keywords internal
#' @useDynLib rGYPSY, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import data.table
"_PACKAGE"

# data.table columns used in GYPSY_Sim()
utils::globalVariables(c(
  ".", "age", "ba", "bhage", "ht", "ps", "si", "si_bh", "si_t", "spp",
  "stand", "stand_age", "stand_idx", "totage", "totage_max", "y2bh",
  paste0(rep(c("aw", "sb", "sw", "pl"), each = 4), "_",
         c("ps", "totage", "si_bh", "bhage"))
))
