# Build the SAS validation fixture used by tests/testthat/test-sas.R
#
# Takes a fixed sample of PGYI interval starts (2026-06-03 export, compiled
# by rPSPab) from the GYPSY_run.sas run of data-raw/rpspab_validation/06_sas.R
# (input in_pgyi.xlsx, output out_pgyi.csv) and saves their inputs and SAS
# outputs to tests/testthat/fixtures/sas_subset.rds. The sample covers every
# species mix, plus up to 3 stands for each GYPSY DLL rule that shows in the
# set (D1, D2, D3/D9, D10), so the tests also check where rGYPSY departs from
# the SAS.
#
# The full data lives in scratch/data/rpspab_validation/sas_run/
# (git-ignored). Run from the package root after 06_sas.R.

library(data.table)
library(openxlsx)
devtools::load_all(".")
source("tests/testthat/helper-sas.R")

run_dir  <- "scratch/data/rpspab_validation/sas_run"
n_stands <- 300

inp <- as.data.table(read.xlsx(file.path(run_dir, "in_pgyi.xlsx")))
inp[, Stan := NULL]

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

sas_csv <- file.path(run_dir, "out_pgyi.csv")
hdr <- names(fread(sas_csv, nrows = 0))
ci <- function(x) {
  m <- hdr[match(tolower(x), tolower(hdr))]
  if (anyNA(m)) stop("Missing SAS columns: ", paste(x[is.na(m)], collapse = ", "))
  m
}
sas <- fread(sas_csv, select = c("Stand", "sect", ci(stand_cols), ci(year_cols)), na.strings = c("", "."))
setnames(sas, tolower(names(sas)))
sas_stand  <- sas[sect == 1, c("stand", tolower(stand_cols)), with = FALSE]
sas_yearly <- sas[, c("stand", "sect", tolower(year_cols)), with = FALSE]

# Stands where a DLL rule applies, found on the whole set
r1 <- run_sas_inputs(inp, stump_dob = 13, top_dib = 7)
rules <- list(D3_D9 = start_rule_stands(sas_stand, r1),
              D10   = beyond_sas_stands(sas_stand, r1),
              D1    = aw_below_bh_stands(r1),
              D2    = interp_stands(sas_yearly, r1))
print(lengths(rules))

# Sample within each species mix so mixed stands are represented
set.seed(20090521)
mixes <- inp[, .(mix = paste(sort(Species), collapse = "-")), by = Stand]
n_all <- nrow(mixes)
pick <- mixes[, .(Stand = if (.N <= 3) Stand else sample(Stand, max(3, round(n_stands*.N/n_all)))), by = mix]$Stand
pick <- sort(unique(c(pick, unlist(lapply(rules, function(s) if (length(s) <= 3) s else sample(s, 3))))))

fixture <- list(
  input  = inp[Stand %in% pick],
  stand  = sas_stand[stand %in% pick],
  yearly = sas_yearly[stand %in% pick]
)

dir.create("tests/testthat/fixtures", recursive = TRUE, showWarnings = FALSE)
saveRDS(fixture, "tests/testthat/fixtures/sas_subset.rds", compress = "xz")
cat(length(pick), "stands,", nrow(fixture$yearly), "yearly rows,",
    round(file.size("tests/testthat/fixtures/sas_subset.rds")/1e6, 1), "MB\n")
for (nm in names(rules)) cat(nm, ":", paste(intersect(rules[[nm]], pick), collapse = ", "), "\n")
