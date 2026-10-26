#' Write the survival model-selection evidence report
#'
#' @param selection Stacked output of [model_selection_table()].
#' @param piecewise Stacked output of [piecewise_table()].
#' @param ph Named list (`os`, `pfs`) of [ph_test()] outputs.
#' @param population List with `start_age`, `prop_male` and `lifetable_edition`.
#' @param figures Named character vector of figure paths relative to the report.
#' @param path Output path.
#' @return `path`.
#' @export
write_survival_report <- function(selection, piecewise, ph, population, figures, path) {
  pct <- function(x) ifelse(is.na(x), "–", sprintf("%.1f", 100 * x))
  num <- function(x, d = 1) ifelse(is.na(x), "–", formatC(x, format = "f", digits = d))
  arm_lab <- c(pembro = "Pembrolizumab", chemo = "Chemotherapy")
  ep_lab <- c(os = "Overall survival", pfs = "Progression-free survival")

  lines <- c(
    "# Survival extrapolation: model-selection evidence",
    "",
    "Models were fitted separately to each arm of the reconstructed KEYNOTE-024 IPD",
    "(Reck et al. 2021, data cut-off 1 June 2020), following NICE DSU TSD 14 and TSD 21.",
    "Candidates: exponential, Weibull, Gompertz, log-logistic, lognormal, generalised",
    "gamma, and Royston–Parmar splines on the hazard scale with 1–3 internal knots.",
    "",
    sprintf(paste("General-population mortality: ONS National Life Tables, England, %s,",
                  "for a cohort aged %g at baseline with %.1f%% male."),
            population$lifetable_edition, population$start_age, 100 * population$prop_male),
    "OS columns marked *constrained* apply this mortality floor. PFS is shown",
    "unconstrained, because its cap at OS depends on the OS model chosen.",
    "",
    "Survival values are percentages. RMST is restricted mean survival to 20 years, in months.",
    "",
    "## Proportional hazards (Grambsch–Therneau test on the treatment term)",
    "",
    "| Endpoint | Chi-squared | p-value |",
    "|---|---|---|",
    vapply(names(ph), function(e) sprintf("| %s | %.2f | %.3f |", ep_lab[[e]], ph[[e]]$chisq,
                                          ph[[e]]$p), character(1)),
    "",
    sprintf("![Log-cumulative hazard, OS](%s)", figures[["lch_os"]]), "",
    sprintf("![Log-cumulative hazard, PFS](%s)", figures[["lch_pfs"]]), ""
  )

  for (ep in c("os", "pfs")) {
    for (a in c("pembro", "chemo")) {
      d <- selection[selection$endpoint == ep & selection$arm == a, ]
      d <- d[order(d$aic), ]
      key <- paste(ep, a, sep = "_")
      lines <- c(lines, sprintf("## %s: %s", ep_lab[[ep]], arm_lab[[a]]), "",
                 sprintf("Observed KM at 5 years: %s%%.", pct(d$km_5y[1])), "")
      if (ep == "os") {
        lines <- c(lines,
          "| Model | AIC | ΔAIC | BIC | ΔBIC | 5y | 10y | 20y | RMST | 5y constrained | 10y constrained | 20y constrained | RMST constrained | Floor binds from (years) |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
          sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
                  d$model, num(d$aic), num(d$delta_aic), num(d$bic), num(d$delta_bic),
                  pct(d$s_5y), pct(d$s_10y), pct(d$s_20y), num(d$rmst_horizon),
                  pct(d$s_5y_constrained), pct(d$s_10y_constrained), pct(d$s_20y_constrained),
                  num(d$rmst_horizon_constrained), num(d$gpm_binding_from_years)))
      } else {
        lines <- c(lines,
          "| Model | AIC | ΔAIC | BIC | ΔBIC | 5y | 10y | 20y | RMST |",
          "|---|---|---|---|---|---|---|---|---|",
          sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s |",
                  d$model, num(d$aic), num(d$delta_aic), num(d$bic), num(d$delta_bic),
                  pct(d$s_5y), pct(d$s_10y), pct(d$s_20y), num(d$rmst_horizon)))
      }
      failed <- d[!d$converged, ]
      if (nrow(failed) > 0) {
        lines <- c(lines, "", paste0("Did not converge: ", failed$model, " (", failed$note, ")"))
      }
      lines <- c(lines, "",
                 sprintf("![Extrapolations](%s)", figures[[paste0("extrap_", key)]]), "",
                 sprintf("![Hazards](%s)", figures[[paste0("haz_", key)]]), "")
    }
  }

  lines <- c(lines, "## TA531-style piecewise scenario (KM to cut-off, exponential tail)", "",
    "| Endpoint | Arm | Cut-off (weeks) | At risk at cut-off | Events after | KM at cut-off | 5y | 10y | 20y | RMST | 5y constrained | 10y constrained | 20y constrained | RMST constrained |",
    "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    sprintf("| %s | %s | %d | %d | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
            toupper(piecewise$endpoint), arm_lab[piecewise$arm], as.integer(piecewise$cutoff_weeks),
            as.integer(piecewise$n_at_cutoff), as.integer(piecewise$events_after),
            pct(piecewise$s_at_cutoff), pct(piecewise$s_5y), pct(piecewise$s_10y),
            pct(piecewise$s_20y), num(piecewise$rmst_horizon),
            pct(piecewise$s_5y_constrained), pct(piecewise$s_10y_constrained),
            pct(piecewise$s_20y_constrained), num(piecewise$rmst_horizon_constrained)),
    "")

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path)
  path
}
