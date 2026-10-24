# Build the SAS validation fixture used by tests/testthat/test-sas.R
#
# Takes a fixed sample of stands from the SAS test run (GYPSY_run.sas on
# test_gypsy_bp.xlsx, output for_gypsy_test_bp.csv) and saves their inputs
# and SAS outputs to tests/testthat/fixtures/sas_subset.rds.
#
# The full data lives in scratch/data/ (git-ignored). Run from the package root.

library(data.table)
library(openxlsx)

data_dir <- "scratch/data"
n_stands <- 300

inp <- as.data.table(read.xlsx(file.path(data_dir, "test_gypsy_bp.xlsx")))

# Sample within each species mix so mixed stands are represented
set.seed(20090521)
mixes <- inp[, .(mix = paste(sort(Species), collapse = "-")), by = Stand]
n_all <- nrow(mixes)
pick <- mixes[, .(Stand = if (.N <= 3) Stand else sample(Stand, max(3, round(n_stands*.N/n_all)))), by = mix]$Stand
pick <- sort(unique(c(pick,
                      372,     # no usable total age (SAS projects all NA)
                      1721)))  # merch volume on a 1e-5 rounding boundary

# SAS output columns compared by the tests
spp <- c("aw", "sb", "sw", "pl")
stand_cols <- c("totage_max", "totage_0",
                unlist(lapply(spp, function(s) {
                  S <- toupper(s)
                  c(paste0(c("SI_t_", "SI_BH_", "totage_", "bhage_", "Y2BH_"), S),
                    paste0(c("psi_", "SDF_", "SC_", "BA_", "rate_"), s, c("", "", "0", "_0", "")))
                })),
                "den_aw_bhage0", "den_sb_totage0", "den_sw_totage0", "den_pl_totage0")
year_cols <- c("totage1", "Standage1",
               unlist(lapply(spp, function(s) {
                 S <- toupper(s)
                 c(paste0(c("totage1_", "bhage1_", "topht_", "den2_", "SC_", "BA_loc_",
                            "mv1_", "mv2_", "mn1_", "mn2_", "mai00_", "mai1_", "mai2_"), s),
                   paste0("SC_", s, "_b"), paste0("BA1_", S), paste0("BA2_", S), paste0("tv_", s, "P"))
               })),
               "v_con_00", "v_con_mv1", "v_con_mv2", "v_all_00", "v_all_mv1", "v_all_mv2",
               "mai00_all", "mai1_all", "mai2_all", "mai00_con", "mai1_con", "mai2_con",
               "mai00_dec", "mai1_dec", "mai2_dec")

hdr <- names(fread(file.path(data_dir, "for_gypsy_test_bp.csv"), nrows = 0))
ci <- function(x) {
  m <- hdr[match(tolower(x), tolower(hdr))]
  if (anyNA(m)) stop("Missing SAS columns: ", paste(x[is.na(m)], collapse = ", "))
  m
}
sas <- fread(file.path(data_dir, "for_gypsy_test_bp.csv"),
             select = c("Stand", "sect", ci(stand_cols), ci(year_cols)),
             na.strings = c("", "."))[Stand %in% pick]
setnames(sas, tolower(names(sas)))

fixture <- list(
  input  = inp[Stand %in% pick],
  stand  = sas[sect == 1, c("stand", tolower(stand_cols)), with = FALSE],
  yearly = sas[, c("stand", "sect", tolower(year_cols)), with = FALSE]
)

dir.create("tests/testthat/fixtures", recursive = TRUE, showWarnings = FALSE)
saveRDS(fixture, "tests/testthat/fixtures/sas_subset.rds", compress = "xz")
cat(length(pick), "stands,", nrow(fixture$yearly), "yearly rows,",
    round(file.size("tests/testthat/fixtures/sas_subset.rds")/1e6, 1), "MB\n")
