#' Survival from a flexsurv fit at given parameters
#'
#' Evaluates the fitted distribution's survival function directly, so that
#' PSA parameter draws can be used in place of the estimates.
#'
#' @param fit A `flexsurvreg` object (including `flexsurvspline` fits).
#' @param t Times in months.
#' @param pars Named parameters on the natural scale, as in `fit$res`
#'   (default: the maximum-likelihood estimates).
#' @return Survival at `t`.
#' @export
surv_from_fit <- function(fit, t, pars = NULL) {
  if (is.null(pars)) pars <- stats::setNames(fit$res[, "est"], rownames(fit$res))
  s <- do.call(fit$dfns$p, c(list(q = t, lower.tail = FALSE), as.list(pars), fit$aux))
  s[t <= 0] <- 1
  s
}

#' Specification of the survival curves for one model run
#'
#' @param os_model,pfs_model Model labels from [candidate_models()], or
#'   `"Piecewise"` (OS only).
#' @param piecewise_weeks Cut-off for the piecewise OS scenario (weeks).
#' @param waning_years Start of treatment-effect waning (years), or `NA`.
#' @param waning_hr Hazard ratio applied to the chemotherapy OS hazard for
#'   pembrolizumab after `waning_years` (1 = full waning).
#' @return A list.
#' @export
curve_spec <- function(os_model = "Lognormal", pfs_model = "Lognormal",
                       piecewise_weeks = NA, waning_years = NA, waning_hr = 1) {
  list(os_model = os_model, pfs_model = pfs_model, piecewise_weeks = piecewise_weeks,
       waning_years = waning_years, waning_hr = waning_hr)
}

#' Apply treatment-effect waning to the pembrolizumab OS curve
#'
#' From `t_wane`, the pembrolizumab OS hazard in each interval equals
#' `hr` times the chemotherapy OS hazard in the same interval.
#'
#' @param s_int,s_comp OS of the intervention and comparator on grid `t`.
#' @param t Grid (months).
#' @param t_wane Waning start (months).
#' @param hr Hazard ratio after waning.
#' @return Intervention OS with waning applied.
#' @export
apply_waning <- function(s_int, s_comp, t, t_wane, hr = 1) {
  dh_int <- diff(-log(pmax(s_int, .Machine$double.xmin)))
  dh_comp <- diff(-log(pmax(s_comp, .Machine$double.xmin)))
  after <- utils::head(t, -1) >= t_wane
  dh_int[after] <- hr * dh_comp[after]
  exp(-c(0, cumsum(dh_int)))
}

#' Build constrained OS and PFS curves for both arms
#'
#' Order of operations: unconstrained curves from the chosen models;
#' treatment-effect waning on pembrolizumab OS (if specified); the
#' general-population mortality floor on each arm's OS; PFS capped at OS
#' ([apply_survival_constraints()]).
#'
#' @param spec Output of [curve_spec()].
#' @param fits Named list (`os_pembro`, `os_chemo`, `pfs_pembro`,
#'   `pfs_chemo`) of [fit_survival_models()] outputs.
#' @param ipd Named list (`os`, `pfs`) of stacked IPD, used for the
#'   piecewise scenario.
#' @param t Grid in months (starting at 0).
#' @param s_gpm General-population survival on `t`.
#' @param pars Optional named list of parameter vectors keyed like `fits`,
#'   for PSA draws of the chosen models.
#' @return A list with `pembro` and `chemo`, each a list with `os`, `pfs`
#'   (constrained) and `os_raw`, `pfs_raw` (unconstrained).
#' @export
build_curves <- function(spec, fits, ipd, t, s_gpm, pars = NULL) {
  raw <- list()
  for (a in c("pembro", "chemo")) {
    for (ep in c("os", "pfs")) {
      key <- paste(ep, a, sep = "_")
      model <- if (ep == "os") spec$os_model else spec$pfs_model
      if (identical(model, "Piecewise")) {
        d <- ipd[[ep]]
        pw <- fit_piecewise(d[d$arm == a, ], weeks_to_months(spec$piecewise_weeks))
        raw[[key]] <- predict_piecewise(pw, t)
      } else {
        fit <- fits[[key]][[model]]
        if (is.null(fit)) stop("Model '", model, "' is not available for ", key)
        raw[[key]] <- surv_from_fit(fit, t, pars[[key]])
      }
    }
  }
  os_p <- raw$os_pembro
  if (!is.na(spec$waning_years)) {
    os_p <- apply_waning(os_p, raw$os_chemo, t, spec$waning_years * 12, spec$waning_hr)
  }
  out <- list()
  for (a in c("pembro", "chemo")) {
    os_raw <- if (a == "pembro") os_p else raw$os_chemo
    cons <- apply_survival_constraints(os_raw, raw[[paste0("pfs_", a)]], s_gpm)
    out[[a]] <- list(os = cons$os, pfs = cons$pfs, os_raw = raw[[paste0("os_", a)]],
                     pfs_raw = raw[[paste0("pfs_", a)]])
  }
  out
}
