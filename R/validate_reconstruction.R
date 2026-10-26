#' Published KEYNOTE-024 five-year targets for reconstruction validation
#'
#' Values from Reck et al. (2021), quoted in `docs/source_checks.md`
#' (R2-R8). Landmarks are included only where the paper reports them
#' explicitly (D-003).
#'
#' @return A list with elements `os` and `pfs`. Each holds `median` and
#'   `events` (named by arm), `hr`, `hr_lower`, `hr_upper` and `landmarks`
#'   (a data frame with `arm`, `time_months` and `surv`).
#' @export
published_targets <- function() {
  list(
    os = list(
      median = c(pembro = 26.3, chemo = 13.4),
      hr = 0.62, hr_lower = 0.48, hr_upper = 0.81,
      events = c(pembro = 103, chemo = 123),
      landmarks = data.frame(arm = c("pembro", "chemo"), time_months = 60,
                             surv = c(0.319, 0.163))
    ),
    pfs = list(
      median = c(pembro = 7.7, chemo = 5.5),
      hr = 0.50, hr_lower = 0.39, hr_upper = 0.65,
      events = c(pembro = 126, chemo = 141),
      landmarks = data.frame(arm = c("pembro", "pembro", "chemo"),
                             time_months = c(36, 60, 36),
                             surv = c(0.228, 0.128, 0.041))
    )
  )
}

#' Acceptance thresholds for the reconstruction
#'
#' Pre-specified in `DECISIONS.md` (D-001, D-003).
#'
#' @return A list with `median_months`, `hr` and `landmark_abs`.
#' @export
acceptance_thresholds <- function() {
  list(median_months = 1, hr = 0.03, landmark_abs = 0.03)
}

#' Validate reconstructed IPD against published results
#'
#' Compares medians, the unstratified Cox HR (Efron ties; chemotherapy as
#' reference, D-005) and the reported landmarks with published values.
#' Also reports unchecked 1-, 2- and 5-year survival and the numbers at
#' risk. Total events are reported without a pass/fail check.
#'
#' @param ipd Data frame with `time`, `status` and `arm` (both arms stacked).
#' @param nar Digitised numbers-at-risk table.
#' @param targets One element of [published_targets()].
#' @param thresholds Output of [acceptance_thresholds()].
#' @param endpoint `"os"` or `"pfs"`, used for labelling.
#' @return A list with `checks` (data frame, one row per check, with
#'   logical `pass`, `NA` for report-only rows), `nar` (published vs
#'   reconstructed numbers at risk) and `pass` (all checks passed).
#' @export
validate_reconstruction <- function(ipd, nar, targets, thresholds, endpoint) {
  ipd$arm <- factor(ipd$arm, levels = c("chemo", "pembro"))
  fit <- survival::survfit(survival::Surv(time, status) ~ arm, data = ipd)
  strata_arm <- sub("^arm=", "", names(fit$strata))

  med <- stats::setNames(summary(fit)$table[, "median"], strata_arm)
  rows <- list()
  for (a in c("pembro", "chemo")) {
    rows[[length(rows) + 1]] <- data.frame(
      endpoint = endpoint, check = sprintf("Median (%s), months", a),
      published = targets$median[[a]], reconstructed = med[[a]],
      abs_diff = abs(med[[a]] - targets$median[[a]]),
      threshold = thresholds$median_months, report_only = FALSE
    )
  }

  cox <- survival::coxph(survival::Surv(time, status) ~ arm, data = ipd, ties = "efron")
  hr <- unname(exp(stats::coef(cox)))
  ci <- unname(exp(stats::confint(cox)))
  rows[[length(rows) + 1]] <- data.frame(
    endpoint = endpoint, check = "HR pembro vs chemo (unstratified Cox)",
    published = targets$hr, reconstructed = hr, abs_diff = abs(hr - targets$hr),
    threshold = thresholds$hr, report_only = FALSE
  )
  rows[[length(rows) + 1]] <- data.frame(
    endpoint = endpoint, check = "HR lower 95% CI",
    published = targets$hr_lower, reconstructed = ci[1],
    abs_diff = abs(ci[1] - targets$hr_lower), threshold = NA_real_, report_only = TRUE
  )
  rows[[length(rows) + 1]] <- data.frame(
    endpoint = endpoint, check = "HR upper 95% CI",
    published = targets$hr_upper, reconstructed = ci[2],
    abs_diff = abs(ci[2] - targets$hr_upper), threshold = NA_real_, report_only = TRUE
  )

  if (!is.null(targets$events)) {
    for (a in c("pembro", "chemo")) {
      ev <- sum(ipd$status[ipd$arm == a])
      rows[[length(rows) + 1]] <- data.frame(
        endpoint = endpoint, check = sprintf("Total events (%s)", a),
        published = targets$events[[a]], reconstructed = ev,
        abs_diff = abs(ev - targets$events[[a]]), threshold = NA_real_, report_only = TRUE
      )
    }
  }

  lm_times <- sort(unique(c(12, 24, 60, targets$landmarks$time_months)))
  for (a in c("pembro", "chemo")) {
    for (t in lm_times) {
      s <- km_at(fit, strata_arm, a, t)
      pub <- targets$landmarks$surv[targets$landmarks$arm == a &
                                      targets$landmarks$time_months == t]
      reported <- length(pub) == 1
      rows[[length(rows) + 1]] <- data.frame(
        endpoint = endpoint, check = sprintf("%g-year survival (%s)", t / 12, a),
        published = if (reported) pub else NA_real_, reconstructed = s,
        abs_diff = if (reported) abs(s - pub) else NA_real_,
        threshold = if (reported) thresholds$landmark_abs else NA_real_,
        report_only = !reported
      )
    }
  }

  checks <- do.call(rbind, rows)
  checks$pass <- ifelse(checks$report_only, NA,
                        !is.na(checks$reconstructed) & checks$abs_diff <= checks$threshold)

  nar_cmp <- do.call(rbind, lapply(c("pembro", "chemo"), function(a) {
    d <- ipd[ipd$arm == a, ]
    data.frame(endpoint = endpoint, arm = a, time_months = nar$time_months,
               published = nar[[a]],
               reconstructed = vapply(nar$time_months, function(t) sum(d$time >= t), 1L))
  }))

  list(checks = checks, nar = nar_cmp, km_fit = fit,
       pass = all(checks$pass, na.rm = TRUE))
}

#' @noRd
km_at <- function(fit, strata_arm, arm, t) {
  idx <- which(strata_arm == arm)
  s <- summary(fit[idx], times = t, extend = TRUE)
  if (length(s$n.risk) == 0 || s$n.risk == 0) return(NA_real_)
  s$surv
}

#' Stop the pipeline if any acceptance check fails
#'
#' @param validations List of outputs from [validate_reconstruction()].
#' @return `TRUE`, invisibly, if all checks pass; otherwise an error listing
#'   the failed checks.
#' @export
assert_reconstruction <- function(validations) {
  checks <- do.call(rbind, lapply(validations, `[[`, "checks"))
  failed <- checks[!is.na(checks$pass) & !checks$pass, ]
  if (nrow(failed) > 0) {
    stop("Reconstruction failed acceptance checks; see output/validation/reconstruction.md:\n",
         paste(sprintf("  %s: %s (published %.3f, reconstructed %.3f, threshold %.3f)",
                       failed$endpoint, failed$check, failed$published,
                       failed$reconstructed, failed$threshold), collapse = "\n"),
         call. = FALSE)
  }
  invisible(TRUE)
}
