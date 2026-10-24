# GYPSY_Sim against the SAS test run (GYPSY_run.sas), 301 stands x 200 years

fx <- sas_fixture()
r1 <- run_sas_inputs(fx$input, stump_dob = 13, top_dib = 7)
r2 <- run_sas_inputs(fx$input, stump_dob = 15, top_dib = 10)

test_that("projection has one row per stand and year", {
  expect_equal(nrow(r1), nrow(fx$yearly))
  expect_setequal(unique(r1$stand), unique(fx$yearly$stand))
  expect_equal(range(r1$year), c(1L, 200L))
})

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
