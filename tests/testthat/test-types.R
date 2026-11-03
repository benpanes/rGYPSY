# SI, age and density types per row (sitype, agetype, ntype)

mix <- list(stand = c(1, 1), stand_age = c(48, 48), spp = c("AW", "SW"),
            n = c(1100, 450), ht = c(19.8, 9.6), ps = c(70, 25), ba = c(21.5, 4.1), simlen = 150)
quiet <- function(f, ...) suppressWarnings(f(...))

test_that("one type per row, all the same, matches the run-wide type", {
  a <- do.call(GYPSY_Sim, c(mix, list(age = c(40, 23), agetype = "B")))
  b <- do.call(GYPSY_Sim, c(mix, list(age = c(40, 23), agetype = c("B", "B"))))
  expect_equal(b, a)
})

test_that("species can use different age types", {
  # SW given at breast height age 23.4; the same stand with SW's total age.
  # (A whole BH age would put a year exactly at BH age 0 in one run and a
  # hair past it in the other, where merchantable stems jump.)
  y2bh <- GYPSY_SI("SW", 9.6, 23.4, agetype = "B", ndx = 3)
  mixed <- do.call(GYPSY_Sim, c(mix, list(age = c(48, 23.4), agetype = c("T", "B"))))
  total <- do.call(GYPSY_Sim, c(mix, list(age = c(48, 23.4 + y2bh), agetype = "T")))
  expect_equal(mixed, total, tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(attr(mixed, "stand_init")$sw_totage, 23.4 + y2bh)
})

test_that("species can use different SI types", {
  si_t <- GYPSY_SI_convert("SW", 16, sitype = "B", ndx = 4)
  base <- mix[setdiff(names(mix), "ht")]
  mixed <- do.call(GYPSY_Sim, c(base, list(si = c(18, si_t), age = c(48, 31), sitype = c("B", "T"))))
  bh    <- do.call(GYPSY_Sim, c(base, list(si = c(18, 16),   age = c(48, 31), sitype = "B")))
  expect_equal(mixed, bh, tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("an initial density (N0) gives that den0 and the matching observed density", {
  base <- mix[setdiff(names(mix), "n")]
  n0 <- do.call(GYPSY_Sim, c(base, list(n = c(1100, 600), age = c(48, 31), ntype = c("N", "N0"))))
  i0 <- attr(n0, "stand_init")
  expect_equal(i0$sw_den0, 600, tolerance = 1e-6)
  # Observed density from that SDF at SW's age, run as an ordinary density
  n_obs <- GYPSY_N("SW", 31, i0$sw_sdf, i0$sw_si_bh, i0$sw_psi, sdf_aw = i0$aw_sdf)
  obs <- do.call(GYPSY_Sim, c(base, list(n = c(1100, n_obs), age = c(48, 31))))
  expect_equal(n0, obs, tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(attr(obs, "stand_init")$sw_den0, 600, tolerance = 1e-6)
})

test_that("types are checked", {
  expect_error(do.call(GYPSY_Sim, c(mix, list(age = c(48, 31), agetype = c("T", "X")))), "agetype must be T, B, or S")
  expect_error(do.call(GYPSY_Sim, c(mix, list(age = c(48, 31), ntype = "M"))), "ntype must be N or N0")
  expect_error(do.call(GYPSY_Sim, c(mix, list(age = c(48, 31), sitype = c("B", "T", "B")))), "sitype must be a single value or length 2")
})
