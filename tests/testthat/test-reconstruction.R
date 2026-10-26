# Round-trip test: simulate IPD with known properties, turn it into a KM
# curve and numbers-at-risk table as a digitiser would see them, reconstruct,
# and check the result passes the same acceptance rules used for KEYNOTE-024.

simulate_km_inputs <- function(n, rate, seed, admin_cens = 60, nar_times = seq(0, 60, 6)) {
  set.seed(seed)
  t_event <- stats::rexp(n, rate)
  t_cens <- stats::runif(n, 40, admin_cens)
  time <- pmin(t_event, t_cens)
  status <- as.integer(t_event <= t_cens)
  fit <- survival::survfit(survival::Surv(time, status) ~ 1)
  ev <- fit$n.event > 0
  # Corner points either side of each drop, as the README instructs.
  s_before <- c(1, fit$surv[ev])[seq_len(sum(ev))]
  curve <- data.frame(
    time_months = c(0, rep(fit$time[ev], each = 2), max(time)),
    surv = c(1, as.vector(rbind(s_before, fit$surv[ev])), utils::tail(fit$surv, 1))
  )
  nar <- data.frame(time_months = nar_times,
                    n = vapply(nar_times, function(t) sum(time >= t), 1L))
  list(curve = curve, nar = nar, ipd = data.frame(time = time, status = status),
       events = sum(status))
}

test_that("Guyot reconstruction recovers simulated data within acceptance thresholds", {
  skip_if_not_installed("IPDfromKM")
  a <- simulate_km_inputs(154, rate = log(2) / 26, seed = 1)
  b <- simulate_km_inputs(151, rate = log(2) / 13, seed = 2)
  nar <- data.frame(time_months = a$nar$time_months, pembro = a$nar$n, chemo = b$nar$n)

  ipd <- rbind(
    cbind(reconstruct_ipd(a$curve, nar, "pembro", 154, a$events)[, c("time", "status")], arm = "pembro"),
    cbind(reconstruct_ipd(b$curve, nar, "chemo", 151, b$events)[, c("time", "status")], arm = "chemo")
  )
  truth <- rbind(cbind(a$ipd, arm = "pembro"), cbind(b$ipd, arm = "chemo"))
  truth$arm <- factor(truth$arm, levels = c("chemo", "pembro"))
  tfit <- survival::survfit(survival::Surv(time, status) ~ arm, data = truth)
  tmed <- summary(tfit)$table[, "median"]
  thr <- exp(stats::coef(survival::coxph(survival::Surv(time, status) ~ arm, data = truth)))

  targets <- list(
    median = c(pembro = tmed[["arm=pembro"]], chemo = tmed[["arm=chemo"]]),
    hr = unname(thr), hr_lower = NA_real_, hr_upper = NA_real_,
    landmarks = data.frame(arm = character(), time_months = numeric(), surv = numeric())
  )
  v <- validate_reconstruction(ipd, nar, targets, acceptance_thresholds(), "sim")
  expect_true(v$pass)
  expect_equal(sum(ipd$status[ipd$arm == "pembro"]), a$events)
  expect_equal(sum(ipd$arm == "pembro"), 154)
})

test_that("a failed check stops the pipeline with an informative error", {
  checks <- data.frame(endpoint = "os", check = "Median (pembro), months",
                       published = 26.3, reconstructed = 24.0, abs_diff = 2.3,
                       threshold = 1, report_only = FALSE, pass = FALSE)
  expect_error(assert_reconstruction(list(os = list(checks = checks))), "Median \\(pembro\\)")
  checks$pass <- TRUE
  expect_true(assert_reconstruction(list(os = list(checks = checks))))
})

test_that("report-only landmarks never fail the acceptance check", {
  checks <- data.frame(endpoint = "os", check = "1-year survival (pembro)",
                       published = NA, reconstructed = 0.7, abs_diff = NA,
                       threshold = NA, report_only = TRUE, pass = NA)
  expect_true(assert_reconstruction(list(os = list(checks = checks))))
})

test_that("published targets match docs/source_checks.md", {
  p <- published_targets()
  expect_equal(p$os$median, c(pembro = 26.3, chemo = 13.4))
  expect_equal(p$pfs$median, c(pembro = 7.7, chemo = 5.5))
  expect_equal(c(p$os$hr, p$pfs$hr), c(0.62, 0.50))
  expect_equal(nrow(p$pfs$landmarks), 3)
})

test_that("KEYNOTE-024 reconstruction passes acceptance checks (real digitised data)", {
  skip_if_not_installed("IPDfromKM")
  dig <- testthat::test_path("..", "..", "data-raw", "digitised")
  files <- file.path(dig, c("os_pembro.csv", "os_chemo.csv", "pfs_pembro.csv", "pfs_chemo.csv",
                            "nar_os.csv", "nar_pfs.csv", "axis_limits.csv"))
  skip_if_not(all(file.exists(files)), "digitised data not yet available")

  axis <- utils::read.csv(file.path(dig, "axis_limits.csv"))
  inputs <- utils::read.csv(testthat::test_path("..", "..", "data-raw", "reconstruction_inputs.csv"))
  vals <- lapply(c("os", "pfs"), function(ep) {
    nar <- read_nar(file.path(dig, paste0("nar_", ep, ".csv")))
    ipd <- do.call(rbind, lapply(c("pembro", "chemo"), function(a) {
      raw <- read_digitised_curve(file.path(dig, paste0(ep, "_", a, ".csv")))
      cv <- clean_digitised_curve(raw, axis_max(axis, ep), paste0(ep, "_", a))$curve
      reconstruct_from_inputs(cv, nar, a, inputs, ep)
    }))
    validate_reconstruction(ipd, nar, published_targets()[[ep]], acceptance_thresholds(), ep)
  })
  expect_true(assert_reconstruction(vals))
})
