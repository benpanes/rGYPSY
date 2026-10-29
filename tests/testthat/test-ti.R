# Tree improvement (Luo and Thomas 2021, section 2.4)

# Rough version the core function follows
si_adj_ref <- function(s, i, t, gw, si){
  aa <- list(SW = c(m = 0.255, c = 0.993), PL = c(m = 0.445, c = 1.017))[[s]]
  ((gw/100)*(aa[["c"]] + aa[["m"]]*log(i/t)) + 1)*si
}

test_that("GYPSY_SI_TI follows the adjusted Lambeth equations", {
  for(s in c("SW", "PL"))
    for(sel in c(6, 15, 25))
      expect_equal(GYPSY_SI_TI(s, 20, 10, sel), si_adj_ref(s, sel, 50, 10, 20))
})

test_that("GYPSY_SI_TI leaves row.s. without a genetic worth unchanged", {
  expect_equal(GYPSY_SI_TI(c("AW", "SB", "SW", "PL"), rep(16, 4), c(NA, 0, 0, 0), rep(NA, 4)), rep(16, 4))
  expect_no_warning(GYPSY_SI_TI(c("AW", "SB"), c(16, 16), c(NA, NA), c(NA, NA)))
})

test_that("a missing SW or PL genetic worth is 0, with a warning", {
  expect_warning(si <- GYPSY_SI_TI(c("AW", "SW", "PL"), rep(16, 3), c(NA, NA, 10), c(NA, NA, 15)),
                 "No genetic worth for SW / PL row(s) 2; treated as 0", fixed = TRUE)
  expect_equal(si[1:2], c(16, 16))
})

test_that("genetic worth is checked", {
  expect_error(GYPSY_SI_TI("AW", 16, 10, 15), "only be given for SW and PL")
  expect_error(GYPSY_SI_TI("SW", 16, 10, NA), "selection age > 0")
  expect_error(GYPSY_Sim_TI(1, 20, "SB", 1000, age = 20, ht = 8, gw = 5, sel_age = 15),
               "only be given for SW and PL")
})

# Mixed stand: AW unimproved, SW and PL improved
stand <- list(stand = c(1, 1, 1), stand_age = rep(30, 3), spp = c("AW", "SW", "PL"),
              n = c(800, 600, 900), age = c(30, 22, 28), ht = c(14, 6.5, 11),
              ps = c(40, 30, 50), ba = c(9, 2.5, 8), simlen = 120)

test_that("GYPSY_Sim_TI without genetic worth is GYPSY_Sim", {
  base <- do.call(GYPSY_Sim, stand)
  expect_warning(imp <- do.call(GYPSY_Sim_TI, stand), "row(s) 2, 3; treated as 0", fixed = TRUE)
  expect_equal(imp, base, ignore_attr = TRUE)
  expect_equal(attr(imp, "stand_init")[, c(aw_gw, sw_gw, pl_gw)], c(NA, 0, 0))
  expect_no_warning(do.call(GYPSY_Sim, stand))   # only TI runs check genetic worth
})

test_that("GYPSY_Sim_TI projects from the improved SI at the observed age", {
  gw <- c(NA, 8, 12); sel <- c(NA, 15, 20)
  imp  <- do.call(GYPSY_Sim_TI, c(stand, list(gw = gw, sel_age = sel)))
  init <- attr(imp, "stand_init")

  base <- attr(do.call(GYPSY_Sim, stand), "stand_init")
  expect_equal(init$sw_si_t_base, base$sw_si_t)
  expect_equal(init$pl_si_bh_base, base$pl_si_bh)
  expect_equal(init$aw_si_t, base$aw_si_t)
  expect_equal(c(init$sw_si_t, init$pl_si_t),
               GYPSY_SI_TI(c("SW", "PL"), c(base$sw_si_t, base$pl_si_t), gw[2:3], sel[2:3]))
  expect_equal(init$pl_aacor, 1.017 + 0.445*log(20/50))
  expect_equal(c(init$sw_totage, init$pl_totage), c(22, 28))

  # Same as GYPSY_Sim given the improved SIs and the observed total ages
  si_b <- c(base$aw_si_bh, init$sw_si_bh, init$pl_si_bh)
  ref  <- do.call(GYPSY_Sim, c(stand[setdiff(names(stand), "ht")], list(si = si_b)))
  expect_equal(imp, ref, ignore_attr = TRUE, tolerance = 1e-8)

  # Improvement raises SW and PL growth
  b <- do.call(GYPSY_Sim, stand)
  grown <- function(x) !is.na(x) & x > 0.1
  expect_true(all((imp$topht_pl > b$topht_pl)[grown(b$topht_pl)]))
  expect_true(all((imp$topht_sw > b$topht_sw)[grown(b$topht_sw)]))
  expect_true(all(imp$mvol_sw[imp$year > 40] > b$mvol_sw[b$year > 40]))
})

test_that("GYPSY_Sim_TI reads gw and sel_age from data", {
  d <- data.frame(stand = 1, stand_age = 15, spp = "PL", n = 2000, si = 18, age = 15,
                  ht = NA, ps = NA, ba = NA, gw = 10, sel_age = 15)
  expect_equal(GYPSY_Sim_TI(data = d, simlen = 50),
               GYPSY_Sim_TI(1, 15, "PL", 2000, si = 18, age = 15, gw = 10, sel_age = 15, simlen = 50))
  expect_warning(expect_warning(GYPSY_Sim_TI(data = d[setdiff(names(d), "gw")], simlen = 5),
                                "not found in data, treated as missing: gw"),
                 "No genetic worth")
})

test_that("a breast height SI is adjusted through the total age SI", {
  si_t <- GYPSY_SI_convert("PL", 18, sitype = "B", ndx = 4)
  b <- GYPSY_Sim_TI(1, 15, "PL", 2000, si = 18,   age = 15, gw = 10, sel_age = 15, simlen = 80)
  t <- GYPSY_Sim_TI(1, 15, "PL", 2000, si = si_t, age = 15, gw = 10, sel_age = 15, simlen = 80,
                    sitype = "T")
  expect_equal(b, t, tolerance = 1e-8)
  expect_equal(attr(b, "stand_init")$pl_si_t, GYPSY_SI_TI("PL", si_t, 10, 15), tolerance = 1e-8)
})
