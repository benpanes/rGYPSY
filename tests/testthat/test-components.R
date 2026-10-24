spp <- c("AW", "SB", "SW", "PL")
ht  <- rep(12, 4)
age <- rep(40, 4)

test_that("SI conversions round trip", {
  si_t  <- GYPSY_SI(spp, ht, age, ndx = 1)
  si_bh <- GYPSY_SI(spp, ht, age, ndx = 2)
  expect_equal(GYPSY_SI_convert(spp, si_bh, sitype = "B", ndx = 4), si_t, tolerance = 1e-8)
  expect_equal(GYPSY_SI_convert(spp, si_t,  sitype = "T", ndx = 4), si_bh, tolerance = 1e-8)
})

test_that("age from SI and height recovers the input age", {
  si_t  <- GYPSY_SI(spp, ht, age, ndx = 1)
  si_bh <- GYPSY_SI(spp, ht, age, ndx = 2)
  expect_equal(GYPSY_Age(spp, si_t,  ht, sitype = "T", ndx = 1), age, tolerance = 1e-8)
  expect_equal(GYPSY_Age(spp, si_bh, ht, sitype = "B", ndx = 1), age, tolerance = 1e-8)
  y2bh <- GYPSY_SI(spp, ht, age, ndx = 3)
  expect_equal(GYPSY_Age(spp, si_t, ht, sitype = "T", ndx = 3), age - y2bh, tolerance = 1e-8)
})

test_that("top height is 1.3 m at Y2BH and 0.3 m at Y2SH", {
  si_t <- GYPSY_SI(spp, ht, age, ndx = 1)
  expect_equal(GYPSY_Ht(spp, si_t, GYPSY_SI(spp, ht, age, ndx = 3)), rep(1.3, 4), tolerance = 1e-8)
  expect_equal(GYPSY_Ht(spp, si_t, GYPSY_SI(spp, ht, age, ndx = 4)), rep(0.3, 4), tolerance = 1e-8)
  expect_equal(GYPSY_Ht(spp, si_t, age), ht, tolerance = 1e-8)
})

test_that("density at the SDF age returns the observed density", {
  sdf <- GYPSY_SDF(1, aw_age_b = 30, aw_n = 1500, aw_si_b = 16, aw_psi = 0,
                   sb_age_t = NA, sb_n = NA, sb_si_b = NA, sb_psi = 0,
                   sw_age_t = NA, sw_n = NA, sw_si_b = NA, sw_psi = 0,
                   pl_age_t = NA, pl_n = NA, pl_si_b = NA, pl_psi = NA)
  expect_equal(GYPSY_N("AW", 30, sdf$aw_sdf, 16), 1500, tolerance = 1e-8)
  expect_equal(GYPSY_N("AW", 0, sdf$aw_sdf, 16), sdf$aw_den0, tolerance = 1e-10)
})

test_that("GYPSY_Sim checks its inputs", {
  expect_error(GYPSY_Sim(1, 20, "XX", 1000, age = 20, ht = 8), "Species must be")
  expect_error(GYPSY_Sim(1, 20, "PL", 1000, age = 20), "2 of SI, Age, and Top Height")
  expect_error(GYPSY_Sim(c(1, 1), 20, "PL", 1000, age = 20, ht = 8), "Duplicate")
  expect_error(GYPSY_Sim(1, 20, "PL", 1000, age = 20, ht = 8, stump_ht = 0.5), "stump_ht")
})

test_that("GYPSY_Sim takes per-species merchantability standards", {
  base <- GYPSY_Sim(1, 60, "PL", 1500, age = 60, ht = 18, simlen = 100)
  alt  <- GYPSY_Sim(1, 60, "PL", 1500, age = 60, ht = 18, simlen = 100,
                    stump_dob = c(AW = 15, SB = 15, SW = 15, PL = 13), top_dib = 10)
  expect_equal(nrow(base), 100)
  expect_equal(base$tvol_pl, alt$tvol_pl)
  expect_true(all(alt$mvol_pl >= base$mvol_pl))
  expect_true(any(alt$mvol_pl > base$mvol_pl))
})
