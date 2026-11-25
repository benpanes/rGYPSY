# GYPSY_Sim against GYPSY_run.sas on a sample of PGYI interval starts
# (2026-06-03 export, compiled by rPSPab), 312 stands x 200 years

fx <- sas_fixture()
r1 <- run_sas_inputs(fx$input, stump_dob = 13, top_dib = 7)
r2 <- run_sas_inputs(fx$input, stump_dob = 15, top_dib = 10)

test_that("projection has one row per stand and year", {
  expect_equal(nrow(r1), nrow(fx$yearly))
  expect_setequal(unique(r1$stand), unique(fx$yearly$stand))
  expect_equal(range(r1$year), c(1L, 200L))
})

# No input BA is given for a species below breast height (D4), which
# GYPSY_Sim would ignore and the SAS would use
skip_ba <- ba_ignored_stands(fx$input, r1)
test_that("no BA is ignored below breast height", {
  expect_length(skip_ba, 0)
})

# Where rGYPSY follows the GYPSY DLL instead of the SAS, the stands differ, so
# they are left out of the comparison. Each check pins the stands of the
# sample where the rule applies:
# * D3, D9: the starting BA (composition and competitors at breast height age
#   0), or a species reaching breast height after the 180 years the SAS
#   searches; their early BA and volume differ too
skip_rule <- start_rule_stands(fx$stand, r1)
test_that("the DLL's starting BA rule changes only the expected stands", {
  expect_equal(skip_rule, c(67, 801, 2431, 2871, 3086, 3152, 3313, 3399, 3805, 4266, 4298))
})
# * D10: leading species older than the SAS's 200 projected years: the SAS has
#   no observed / predicted BA ratio, rGYPSY (as the DLL) has
skip_beyond <- beyond_sas_stands(fx$stand, r1)
test_that("only stands older than the SAS projection get a ratio the SAS lacks", {
  expect_equal(skip_beyond, c(18, 181, 2042, 2466, 2728, 2792))
})
# * D1: aspen just below breast height, grown from breast height age -1
skip_aw <- aw_below_bh_stands(r1)
test_that("only the expected stands have aspen grown below breast height", {
  expect_equal(skip_aw, c(67, 3152, 3399))
})
# * D2: BA interpolated over the ratio cap while the lagged predicted BA is 0
skip_interp <- interp_stands(fx$yearly, r1)
test_that("only the expected stands get BA the DLL interpolates", {
  expect_equal(skip_interp, c(67, 2570, 2782, 3152, 3399))
})
skip <- c(skip_ba, skip_rule, skip_beyond, skip_aw, skip_interp)
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
