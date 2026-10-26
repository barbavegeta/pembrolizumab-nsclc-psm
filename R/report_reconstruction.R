#' Plot reconstructed Kaplan-Meier curves over the digitised points
#'
#' @param ipd Stacked reconstructed IPD (`time`, `status`, `arm`).
#' @param digitised Stacked cleaned digitised curves (`time_months`, `surv`,
#'   `arm`).
#' @param endpoint_label Axis label, e.g. `"Overall survival"`.
#' @return A ggplot object.
#' @export
plot_reconstruction <- function(ipd, digitised, endpoint_label) {
  fit <- survival::survfit(survival::Surv(time, status) ~ arm, data = ipd)
  km <- data.frame(
    time_months = c(0, 0, fit$time),
    surv = c(1, 1, fit$surv),
    arm = c(sub("^arm=", "", names(fit$strata)),
            rep(sub("^arm=", "", names(fit$strata)), fit$strata))
  )
  km <- km[order(km$arm, km$time_months, -km$surv), ]
  arm_labels <-c(pembro = "Pembrolizumab", chemo = "Chemotherapy")
  ggplot2::ggplot() +
    ggplot2::geom_step(data = km, ggplot2::aes(time_months, surv, colour = arm),
                       linewidth = 0.6) +
    ggplot2::geom_point(data = digitised, ggplot2::aes(time_months, surv, colour = arm),
                        size = 0.7, alpha = 0.6, shape = 16) +
    ggplot2::scale_colour_manual(values = c(pembro = "#1b6ca8", chemo = "#c8553d"),
                                 labels = arm_labels, name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(x * 100, "%")) +
    ggplot2::labs(x = "Time since randomisation (months)", y = endpoint_label,
                  caption = paste("Lines: KM of reconstructed IPD. Points: digitised",
                                  "coordinates from Reck et al. 2021, Fig 2.")) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "top")
}

#' Write the reconstruction validation report
#'
#' Writes a Markdown file with the acceptance thresholds, the result of
#' every check, the numbers-at-risk comparison, the cleaning log and links
#' to the overlay plots. The file is written whether or not the checks
#' pass. [assert_reconstruction()] then stops the pipeline on failure.
#'
#' @param validations Named list (`os`, `pfs`) of [validate_reconstruction()]
#'   outputs.
#' @param cleaning_log Combined cleaning log from [clean_digitised_curve()].
#' @param methods Named character vector giving the implementation used for
#'   each curve.
#' @param plot_files Named character vector of overlay plot paths, relative
#'   to the report.
#' @param path Output path.
#' @return `path`, for use as a `targets` file target.
#' @export
write_reconstruction_report <- function(validations, cleaning_log, methods,
                                        plot_files, path) {
  thr <- acceptance_thresholds()
  fmt <- function(x, d = 3) {
    d <- rep_len(d, length(x))
    vapply(seq_along(x), function(i) {
      if (is.na(x[i])) "–" else formatC(x[i], format = "f", digits = d[i])
    }, character(1))
  }
  status <- function(p) ifelse(is.na(p), "report only", ifelse(p, "PASS", "**FAIL**"))
  all_pass <- all(vapply(validations, `[[`, logical(1), "pass"))

  lines <- c(
    "# IPD reconstruction: validation",
    "",
    "Reconstructed from digitised Figure 2A (OS) and 2B (PFS) of Reck et al.,",
    "*J Clin Oncol* 2021;39(21):2339–2349 (data cut-off 1 June 2020), with the",
    "Guyot et al. (2012) algorithm.",
    "",
    sprintf("**Overall result: %s**", if (all_pass) "all acceptance checks passed" else
      "one or more acceptance checks FAILED, so the pipeline stops"),
    "",
    "## Acceptance thresholds (pre-specified; DECISIONS.md D-001, D-003)",
    "",
    sprintf("- Median survival: within %g month of the published value, in each arm.", thr$median_months),
    sprintf("- Hazard ratio: within %g of the published HR.", thr$hr),
    sprintf("- Landmark survival: within %g percentage points, only where Reck 2021 reports the value.",
            thr$landmark_abs * 100),
    "",
    "The published HRs come from a Cox model stratified by ECOG PS, histology and",
    "region. The reconstructed data contain no stratification factors, so the",
    "reconstructed HR comes from an unstratified Cox model (Efron ties). Both",
    "estimates are intention-to-treat and neither adjusts for crossover (D-005).",
    "",
    "## Implementation",
    "",
    paste0("- ", names(methods), ": `", methods, "`"),
    ""
  )

  for (ep in names(validations)) {
    v <- validations[[ep]]
    ck <- v$checks
    is_pct <- grepl("survival", ck$check)
    is_n <- grepl("events", ck$check)
    tab <- c(
      "| Check | Published | Reconstructed | Abs. difference | Threshold | Result |",
      "|---|---|---|---|---|---|",
      sprintf("| %s | %s | %s | %s | %s | %s |", ck$check,
              ifelse(is_pct, fmt(ck$published * 100, 1), fmt(ck$published, ifelse(is_n, 0, 2))),
              ifelse(is_pct, fmt(ck$reconstructed * 100, 1), fmt(ck$reconstructed, ifelse(is_n, 0, 2))),
              ifelse(is_pct, fmt(ck$abs_diff * 100, 1), fmt(ck$abs_diff, ifelse(is_n, 0, 3))),
              ifelse(is_pct, fmt(ck$threshold * 100, 1), fmt(ck$threshold, 2)),
              status(ck$pass))
    )
    nr <- v$nar
    nar_tab <- c(
      "| Arm | Time (months) | Published at risk | Reconstructed at risk |",
      "|---|---|---|---|",
      sprintf("| %s | %g | %d | %d |", nr$arm, nr$time_months, nr$published, nr$reconstructed)
    )
    lines <- c(lines,
      sprintf("## %s", toupper(ep)), "",
      "Survival values are percentages.", "",
      tab, "",
      sprintf("![%s overlay](%s)", toupper(ep), plot_files[[ep]]), "",
      sprintf("### %s numbers at risk", toupper(ep)), "",
      nar_tab, "")
  }

  lines <- c(lines, "## Cleaning log (DECISIONS.md D-002, D-006)", "")
  if (nrow(cleaning_log) == 0) {
    lines <- c(lines, "No cleaning changes were needed.")
  } else {
    lines <- c(lines,
      "| Curve | Rule | Points | Detail |", "|---|---|---|---|",
      sprintf("| %s | %s | %d | %s |", cleaning_log$curve, cleaning_log$rule,
              cleaning_log$n_points, cleaning_log$detail))
  }

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path)
  path
}
