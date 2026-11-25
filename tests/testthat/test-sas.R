# GYPSY_Sim against the SAS test run (GYPSY_run.sas), 301 stands x 200 years

fx <- sas_fixture()
r1 <- run_sas_inputs(fx$input, stump_dob = 13, top_dib = 7)
r2 <- run_sas_inputs(fx$input, stump_dob = 15, top_dib = 10)

test_that("projection has one row per stand and year", {
  expect_equal(nrow(r1), nrow(fx$yearly))
  expect_setequal(unique(r1$stand), unique(fx$yearly$stand))
  expect_equal(range(r1$year), c(1L, 200L))
})

# Stand 8413 gives a BA for PL at a top height of 1.3 m (breast height age 0).
# GYPSY_Sim ignores it; the SAS uses it, so the stand is not compared.
skip_ba <- ba_ignored_stands(fx$input, r1)
test_that("only stand 8413 has a BA ignored below breast height", {
  expect_equal(skip_ba, 8413)
})

# rGYPSY follows the GYPSY DLL's starting BA rule (D3), not the SAS's: in these
# stands the composition at breast height age 0 counts, or the SW / PL BA
# increment's competitors include, a species the SAS leaves out. Their early
# BA and volume then differ too, so they are not compared.
skip_rule <- start_rule_stands(fx$stand, r1)
test_that("the DLL's starting BA rule changes only the expected stands", {
  expect_equal(skip_rule, c(869, 1879, 1899, 1930, 1963, 1980, 2007, 4455, 4782, 8829))
})
skip <- c(skip_ba, skip_rule)
fx$stand  <- fx$stand[!fx$stand$stand %in% skip, ]
fx$yearly <- fx$yearly[!fx$yearly$stand %in% skip, ]

test_that("stand-level SI, ages, PSI, SDF, initial density, SC0, BA0 match SAS", {
  m <- align_sas(fx$stand, attr(r1, "stand_init"), "stand", "stand")
  expect_equal(nrow(m), nrow(fx$stand))
  expect_sas_match(m, stand_map())
})

test_that("yearly ages, height, density, composition, BA and volume match SAS (standard 0.3/13/7)", {
  m <- align_sas(fx$yearly, r1, c("stand", "sect"), c("stand", "year"))
  expect_equal(nrow(m), nrow(fx$yearly))
  expect_sas_match(m, yearly_map(1))
})

test_that("merchantable volume, density and MAI match SAS (standard 0.3/15/10)", {
  m <- align_sas(fx$yearly, r2, c("stand", "sect"), c("stand", "year"))
  expect_sas_match(m, yearly_map(2))
})
