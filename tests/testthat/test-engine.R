# Hand-calculated toy model: three annual cycles.
#   PFS at boundaries: 1, 0.8, 0.5, 0.3; OS: 1, 0.9, 0.7, 0.5
#   Trapezoidal occupancy: PF 0.90, 0.65, 0.40; PD 0.05, 0.15, 0.20
#   Utilities PF 0.8, PD 0.6 -> QALYs per cycle 0.75, 0.61, 0.44 (sum 1.80)
#   State costs PF 100, PD 50 per cycle, plus a one-off 1000 in cycle 1
#     -> costs per cycle 1092.5, 72.5, 50.0 (sum 1215)
#   Discount factors at midpoints, 3.5%: 1.035^-0.5, 1.035^-1.5, 1.035^-2.5
#     = 0.98294637, 0.94970664, 0.91759096
#   Discounted QALYs 1.72027085; discounted costs 1188.60219345;
#   LYs 2.35 (discounted 2.24411894)

toy <- function() {
  run_psm(s_pfs = c(1, 0.8, 0.5, 0.3), s_os = c(1, 0.9, 0.7, 0.5), cycle_years = 1,
          u_pf = 0.8, u_pd = 0.6, c_pf = 100, c_pd = 50,
          cost_streams = list(one_off = c(1000, 0, 0)), rate_costs = 0.035, rate_qalys = 0.035)
}

test_that("three-cycle toy PSM matches hand calculation", {
  r <- toy()
  expect_equal(r$trace$occ_pf, c(0.90, 0.65, 0.40))
  expect_equal(r$trace$occ_pd, c(0.05, 0.15, 0.20))
  expect_equal(r$trace$qaly, c(0.75, 0.61, 0.44))
  expect_equal(r$trace$cost_total, c(1092.5, 72.5, 50.0))
  expect_equal(r$totals[["qaly_undisc"]], 1.80)
  expect_equal(r$totals[["qaly_disc"]], 1.72027085, tolerance = 1e-8)
  expect_equal(r$totals[["cost_undisc"]], 1215)
  expect_equal(r$totals[["cost_disc"]], 1188.60219345, tolerance = 1e-8)
  expect_equal(r$totals[["ly_undisc"]], 2.35)
  expect_equal(r$totals[["ly_disc"]], 2.24411894, tolerance = 1e-8)
})

test_that("discounting and half-cycle correction match hand values", {
  expect_equal(discount_midpoint(3, 1, 0.035), c(0.98294637, 0.94970664, 0.91759096),
               tolerance = 1e-8)
  expect_equal(discount_midpoint(2, 0.5, 0), c(1, 1))
  expect_equal(half_cycle_trapezoid(c(1, 0.8, 0.5)), c(0.9, 0.65))
})

test_that("state membership sums to 1 in every cycle of the base case", {
  skip_if_not(file.exists(testthat::test_path("..", "..", "_targets", "meta", "meta")),
              "pipeline has not been run")
  store <- testthat::test_path("..", "..", "_targets")
  base <- targets::tar_read(base_case, store = store)
  for (a in c("pembro", "chemo")) {
    tr <- base$run[[a]]$trace
    expect_equal(tr$pf_start + tr$pd_start + tr$dead_start, rep(1, nrow(tr)), tolerance = 1e-12)
    expect_true(all(tr$pf_start >= -1e-12 & tr$pd_start >= -1e-12 & tr$dead_start >= -1e-12))
  }
})

test_that("psm_states rejects PFS above OS and sums to one", {
  st <- psm_states(c(1, 0.6, 0.2), c(1, 0.7, 0.4))
  expect_equal(rowSums(st), c(1, 1, 1))
  expect_error(psm_states(c(1, 0.8), c(1, 0.7)))
})

test_that("dosing stream doses on schedule and respects the stopping rule", {
  t_weeks <- 0:12
  s <- rep(1, 13)
  d <- dosing_stream(s, t_weeks, start_week = 0, interval_weeks = 3, max_doses = 3)
  expect_equal(which(d > 0), c(1, 4, 7))
  d2 <- dosing_stream(s, t_weeks, start_week = 6, interval_weeks = 3)
  expect_equal(which(d2 > 0), c(7, 10))
})

test_that("unit costs follow from the inputs", {
  data <- load_model_data(testthat::test_path("..", "..", "data-raw"))
  p <- model_params(data)
  uc <- unit_costs(p)
  expect_equal(uc$pembro_drug_dose, 2 * 2630 * 0.9921)
  expect_equal(uc$pembro_dose, 2 * 2630 * 0.9921 + p$admin_sb12z_op)
  p$price_discount_pembro <- 0.5
  expect_equal(unit_costs(p)$pembro_drug_dose, 2630 * 0.9921)
  expect_equal(p$inflation_201415, prod(1 + c(0.40, 2.09, 1.24, 1.60, 2.14, 2.49, 2.58, 7.32,
                                               2.47, 4.02) / 100))
})

test_that("inputs.csv prices match the prepared eMIT file", {
  dir <- testthat::test_path("..", "..", "data-raw")
  inp <- read_inputs(file.path(dir, "inputs.csv"))
  emit <- utils::read.csv(file.path(dir, "emit_2025_26.csv"))
  for (d in emit$drug) {
    expect_equal(input_value(inp, paste0("price_per_mg_", d)), emit$price_per_mg[emit$drug == d],
                 tolerance = 1e-4, info = d)
  }
})

test_that("time-to-death utility equals the band utility when death timing is certain", {
  t <- seq(0, 40, by = 0.25)
  u <- list(u_ttd_0 = 0.4, u_ttd_30 = 0.6, u_ttd_180 = 0.7, u_ttd_360 = 0.8)
  # No deaths at all: everyone is 360+ days from death.
  expect_equal(ttd_utility(rep(1, length(t)), t, c(1, 5), u), c(0.8, 0.8))
  # Constant hazard: probabilities follow from exp(-h * d).
  h <- 0.05
  s <- exp(-h * t)
  d2m <- function(d) d / (365.25 / 12)
  pr <- c(1 - exp(-h * d2m(30)), exp(-h * d2m(30)) - exp(-h * d2m(180)),
          exp(-h * d2m(180)) - exp(-h * d2m(360)), exp(-h * d2m(360)))
  expect_equal(ttd_utility(s, t, 2, u), sum(pr * unlist(u)), tolerance = 1e-6)
})

test_that("waning sets the pembrolizumab hazard to HR x comparator hazard", {
  t <- seq(0, 120, by = 1)
  s_int <- exp(-0.01 * t)
  s_comp <- exp(-0.03 * t)
  w <- apply_waning(s_int, s_comp, t, t_wane = 36, hr = 1)
  expect_equal(w[t <= 36], s_int[t <= 36])
  dh <- diff(-log(w))
  expect_equal(dh[t[-length(t)] >= 36], rep(0.03, sum(t[-length(t)] >= 36)))
})

test_that("PSA is reproducible with a fixed seed and gives finite results", {
  skip_if_not(file.exists(testthat::test_path("..", "..", "_targets", "meta", "meta")),
              "pipeline has not been run")
  store <- testthat::test_path("..", "..", "_targets")
  p <- targets::tar_read(params, store = store)
  data <- targets::tar_read(model_data, store = store)
  fits <- targets::tar_read(fits_by_curve, store = store)
  settings <- targets::tar_read(base_settings, store = store)
  ctx <- targets::tar_read(ctx, store = store)
  d1 <- draw_psa(p, data, fits, settings$curves, n = 5, seed = 99)
  d2 <- draw_psa(p, data, fits, settings$curves, n = 5, seed = 99)
  expect_identical(d1, d2)
  r1 <- run_psa(d1, settings, ctx)
  r2 <- run_psa(d2, settings, ctx)
  expect_identical(r1, r2)
  expect_true(all(is.finite(as.matrix(r1))))
  psa <- targets::tar_read(psa, store = store)
  expect_true(all(is.finite(psa$inc_cost)) && all(is.finite(psa$inc_qaly)))
})
