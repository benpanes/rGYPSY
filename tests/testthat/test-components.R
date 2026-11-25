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
  expect_error(GYPSY_Sim(c(1, 1), c(20, 20), c("PL", "PL"), c(1000, 1000),
                         age = c(20, 20), ht = c(8, 8)), "Duplicate")
  expect_error(GYPSY_Sim(1, 20, "PL", 1000, age = 20, ht = 8, stump_ht = 0.5), "stump_ht")
})

test_that("single values are repeated, warning only for supplied data inputs", {
  expect_warning(h <- GYPSY_Ht("SW", si_t = 15, age_t = c(20, 40)), "spp, si_t")
  expect_equal(h, GYPSY_Ht(c("SW", "SW"), c(15, 15), c(20, 40)))
  expect_warning(GYPSY_SI(spp, 12, age), "ht")
  expect_warning(GYPSY_Age("PL", si = c(15, 18), ht = c(10, 10)), "spp")
  expect_warning(GYPSY_SI_convert("PL", si = c(15, 18)), "spp")
  # Defaults (psi, competitor SDFs) and merch standards repeat silently
  expect_no_warning(GYPSY_N(c("PL", "PL"), c(10, 30), c(1800, 1800), c(17, 17)))
  expect_no_warning(GYPSY_Vol(c("SW", "PL"), c(20, 30), c(18, 20), c(800, 1500),
                              c(55, 60), stump_dob = 13))
  expect_warning(GYPSY_Vol(c("SW", "PL"), 20, c(18, 20), c(800, 1500), c(55, 60)),
                 "2 rows: ba$")
  expect_error(suppressWarnings(GYPSY_Ht(spp, si_t = c(15, 16), age_t = 40)),
               "si_t must be a single value or length 4")
})

test_that("GYPSY_Sim reads inputs from a data table", {
  inp  <- data.table::fread(system.file("extdata", "example_input.csv", package = "rGYPSY"))
  base <- GYPSY_Sim(inp$stand, inp$stand_age, inp$spp, inp$n, si = inp$si,
                    age = inp$age, ht = inp$ht, ps = inp$ps, ba = inp$ba, simlen = 50)
  # Extra columns are ignored
  expect_equal(GYPSY_Sim(data = cbind(inp, plot = "x"), simlen = 50), base)
  # Arguments given directly override columns
  expect_equal(GYPSY_Sim(data = inp, ba = rep(NA, nrow(inp)), simlen = 50),
               GYPSY_Sim(data = inp, baflag = 0, simlen = 50))
  # Missing optional columns warn, missing required columns stop
  expect_warning(GYPSY_Sim(data = inp[, !c("ps", "ba")], simlen = 50),
                 "not found in data, treated as missing: ps, ba")
  expect_error(GYPSY_Sim(data = inp[, !"n"]), "Required column\\(s\\) not found in data: n")
  expect_error(GYPSY_Sim(inp), "use GYPSY_Sim\\(data = ...\\)")
})

test_that("GYPSY_Sim ignores BA for species below breast height", {
  # PL at total age 5 with SI 18 is below breast height; SW at 30 is above it
  args <- list(stand = c(1, 1), stand_age = c(30, 30), spp = c("SW", "PL"),
               n = c(900, 1500), si = c(16, 18), age = c(30, 5), simlen = 60)
  init <- attr(do.call(GYPSY_Sim, c(args, list(baflag = 0))), "stand_init")
  expect_gt(init$sw_bhage, 0)
  expect_lte(init$pl_bhage, 0)

  expect_warning(out <- do.call(GYPSY_Sim, c(args, list(ba = c(8, 0.5)))),
                 "BA ignored for row\\(s\\) 2: .*breast height")
  expect_true(is.na(attr(out, "stand_init")$pl_ba))
  expect_equal(attr(out, "stand_init")$sw_ba, 8)
  expect_equal(out, suppressWarnings(do.call(GYPSY_Sim, c(args, list(ba = c(8, NA))))),
               ignore_attr = TRUE)
  # A zero BA is consistent with being below breast height: no warning
  expect_no_warning(do.call(GYPSY_Sim, c(args, list(ba = c(8, 0)))))
})

test_that("GYPSY_Sim warns about densities of 30 stems/ha or fewer, and still projects them", {
  args <- list(stand = c(1, 1), stand_age = c(60, 60), spp = c("SW", "PL"),
               si = c(16, 18), age = c(60, 55), simlen = 40)
  expect_warning(out <- do.call(GYPSY_Sim, c(args, list(n = c(30, 1200)))),
                 "30 stems/ha or fewer in row\\(s\\) 1:")
  expect_true(any(out$tvol_sw > 0))
  expect_no_warning(do.call(GYPSY_Sim, c(args, list(n = c(31, 1200)))))
})

test_that("qmd_ is the quadratic mean diameter", {
  out <- GYPSY_Sim(1, 60, "PL", 1500, age = 60, ht = 18, simlen = 50)
  ok  <- !is.na(out$qmd_pl) & !is.na(out$ba_pl)
  expect_true(any(ok))
  expect_equal(pi*out$qmd_pl[ok]^2/40000*out$n_pl[ok], out$ba_pl[ok])
})

test_that("aspen just below breast height is grown, as the GYPSY DLL does", {
  # SI_bh 18 gives 3.22096 years to breast height. Reference values from the
  # GYPSY 2009 DLL (aspen alone, 2000 stems/ha, stand age 20)
  y2bh <- 3.22096052110905
  out  <- GYPSY_Sim(1, 20, "AW", 2000, si = 18, age = y2bh - 0.5, simlen = 60)
  init <- attr(out, "stand_init")
  expect_equal(init$aw_bhage, -0.5, tolerance = 1e-6)
  expect_equal(init$aw_sdf, 612.14631912940092, tolerance = 1e-8)
  expect_equal(init$aw_den0, 1998.54367505258, tolerance = 1e-8)
  expect_equal(out[abs(standage - 60) < 1e-6, mvol_aw], 53.588098992342836, tolerance = 1e-8)
  # Below BH age -1 the density model is undefined (log(age + 1)): not grown
  out <- GYPSY_Sim(1, 20, "AW", 2000, si = 18, age = y2bh - 1.001, simlen = 60)
  expect_true(is.na(attr(out, "stand_init")$aw_sdf))
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
