toy_lifetable <- function() {
  ages <- 0:100
  rbind(data.frame(age = ages, sex = "male", qx = pmin(0.001 * 1.1^ages, 0.5)),
        data.frame(age = ages, sex = "female", qx = pmin(0.0006 * 1.1^ages, 0.5)))
}

test_that("general-population survival matches a hand calculation", {
  lt <- toy_lifetable()
  # From exact age 65 for 18 months: 12 months at the age-65 hazard, 6 at age 66.
  h65m <- -log(1 - lt$qx[lt$sex == "male" & lt$age == 65])
  h66m <- -log(1 - lt$qx[lt$sex == "male" & lt$age == 66])
  h65f <- -log(1 - lt$qx[lt$sex == "female" & lt$age == 65])
  h66f <- -log(1 - lt$qx[lt$sex == "female" & lt$age == 66])
  expected <- 0.6 * exp(-(h65m + 0.5 * h66m)) + 0.4 * exp(-(h65f + 0.5 * h66f))
  expect_equal(gpm_survival(lt, 65, 0.6, 18), expected, tolerance = 1e-12)
  expect_equal(gpm_survival(lt, 65, 0.6, 0), 1)
  # A non-integer starting age splits the first year correctly.
  s_half <- gpm_survival(lt, 65.5, 1, 6)
  expect_equal(s_half, exp(-0.5 * h65m), tolerance = 1e-12)
})

test_that("constraints: OS hazard >= general-population hazard, PFS <= OS, monotone", {
  grid <- extrapolation_grid(20)
  s_gpm <- gpm_survival(toy_lifetable(), 65, 0.646, grid)
  # A model whose hazard falls towards zero, so the floor must bind.
  s_os <- exp(-0.5 * (1 - exp(-grid / 12)))
  s_pfs <- pmin(1, exp(-0.02 * grid) * 1.05)
  s_pfs[1] <- 1
  res <- apply_survival_constraints(s_os, s_pfs, s_gpm)

  dh_os <- diff(-log(res$os))
  dh_gpm <- diff(-log(s_gpm))
  expect_true(all(dh_os >= dh_gpm - 1e-12))
  expect_true(any(res$os_gpm_binding))
  expect_true(all(res$pfs <= res$os))
  expect_true(all(diff(res$os) <= 0))
  expect_true(all(diff(res$pfs) <= 0))
  expect_equal(res$os[1], 1)
})

test_that("constraints leave curves unchanged when they do not bind", {
  grid <- seq(0, 60, by = 1)
  s_gpm <- exp(-0.0001 * grid)
  s_os <- exp(-0.05 * grid)
  s_pfs <- exp(-0.1 * grid)
  res <- apply_survival_constraints(s_os, s_pfs, s_gpm)
  expect_equal(res$os, s_os, tolerance = 1e-12)
  expect_equal(res$pfs, s_pfs)
  expect_false(any(res$os_gpm_binding))
})

test_that("piecewise model equals KM up to the cut-off and is continuous after it", {
  set.seed(42)
  ipd <- data.frame(time = rexp(200, 0.05), status = rbinom(200, 1, 0.8))
  cutoff <- weeks_to_months(33)
  pw <- fit_piecewise(ipd, cutoff)
  km <- survival::survfit(survival::Surv(time, status) ~ 1, data = ipd)
  t_pre <- seq(0, cutoff, length.out = 25)
  expect_equal(predict_piecewise(pw, t_pre), summary(km, times = t_pre, extend = TRUE)$surv)
  expect_equal(predict_piecewise(pw, cutoff + 1e-8), pw$s_cutoff, tolerance = 1e-6)
  t_post <- cutoff + c(1, 10, 50)
  rate <- pw$tail$res["rate", "est"]
  expect_equal(predict_piecewise(pw, t_post), pw$s_cutoff * exp(-rate * (t_post - cutoff)),
               tolerance = 1e-10)
})

test_that("all candidate models fit and extrapolations are valid survival curves", {
  set.seed(7)
  ipd <- data.frame(time = rweibull(300, 1.2, 20), status = rbinom(300, 1, 0.85))
  fits <- fit_survival_models(ipd)
  expect_length(fits, 9)
  expect_true(all(is.na(attr(fits, "errors"))))
  grid <- extrapolation_grid(20)
  for (m in names(fits)) {
    s <- predict_fit(fits[[m]], grid)
    expect_true(all(diff(s) <= 1e-12), info = m)
    expect_true(all(s >= 0 & s <= 1), info = m)
  }
})

test_that("weeks convert to months consistently with the weekly grid", {
  expect_equal(weeks_to_months(365.25 / 7), 12)
  grid <- extrapolation_grid(20)
  expect_true(all(c(60, 120, 240) %in% grid))
  expect_equal(max(grid), 240)
})
