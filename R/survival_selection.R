#' Time grid for extrapolation
#'
#' Weekly grid from 0 to `horizon_years`, with the landmark times added.
#'
#' @param horizon_years Horizon in years.
#' @param landmarks_years Landmark times in years.
#' @return Sorted grid in months.
#' @export
extrapolation_grid <- function(horizon_years = 20, landmarks_years = c(5, 10, 20)) {
  h <- horizon_years * 12
  sort(unique(c(seq(0, h, by = weeks_to_months(1)), landmarks_years * 12, h)))
}

#' Model-selection table for one endpoint and arm
#'
#' Reports AIC, BIC and their differences from the best model, plus
#' extrapolated survival at the landmarks and restricted mean survival to
#' the horizon. OS is shown both unconstrained and with the
#' general-population mortality floor. PFS is shown unconstrained, because
#' its cap at OS depends on which OS model is chosen.
#'
#' @param fits Output of [fit_survival_models()].
#' @param ipd IPD for the arm (`time`, `status`).
#' @param endpoint `"os"` or `"pfs"`.
#' @param arm Arm label.
#' @param s_gpm General-population survival on `grid`.
#' @param grid Output of [extrapolation_grid()].
#' @param landmarks_years Landmark times in years.
#' @return A data frame with one row per model.
#' @export
model_selection_table <- function(fits, ipd, endpoint, arm, s_gpm, grid,
                                  landmarks_years = c(5, 10, 20)) {
  km <- survival::survfit(survival::Surv(time, status) ~ 1, data = ipd)
  km5 <- summary(km, times = 60, extend = TRUE)
  km5 <- if (km5$n.risk > 0) km5$surv else NA_real_
  errors <- attr(fits, "errors")
  rows <- lapply(names(fits), function(m) {
    f <- fits[[m]]
    base <- data.frame(endpoint = endpoint, arm = arm, model = m)
    if (is.null(f)) {
      return(cbind(base, aic = NA, bic = NA, converged = FALSE, note = errors[[m]]))
    }
    s <- predict_fit(f, grid)
    lm_idx <- match(landmarks_years * 12, grid)
    out <- cbind(base, aic = stats::AIC(f), bic = stats::BIC(f), converged = TRUE, note = NA)
    for (k in seq_along(landmarks_years)) {
      out[[sprintf("s_%gy", landmarks_years[k])]] <- s[lm_idx[k]]
    }
    out$rmst_horizon <- rmst_trapezoid(s, grid)
    if (endpoint == "os") {
      sc <- apply_survival_constraints(s, s_gpm = s_gpm)$os
      for (k in seq_along(landmarks_years)) {
        out[[sprintf("s_%gy_constrained", landmarks_years[k])]] <- sc[lm_idx[k]]
      }
      out$rmst_horizon_constrained <- rmst_trapezoid(sc, grid)
      binding <- apply_survival_constraints(s, s_gpm = s_gpm)$os_gpm_binding
      out$gpm_binding_from_years <- if (any(binding)) grid[which(binding)[1]] / 12 else NA_real_
    }
    out
  })
  cols <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(r) { r[setdiff(cols, names(r))] <- NA; r[cols] })
  tab <- do.call(rbind, rows)
  tab$delta_aic <- tab$aic - min(tab$aic, na.rm = TRUE)
  tab$delta_bic <- tab$bic - min(tab$bic, na.rm = TRUE)
  tab$km_5y <- km5
  tab
}

#' Piecewise (TA531-style) scenario table for one endpoint and arm
#'
#' @param ipd IPD for the arm.
#' @param endpoint,arm Labels.
#' @param cutoffs_weeks Cut-offs in weeks.
#' @param s_gpm,grid As in [model_selection_table()].
#' @param dist Tail distribution.
#' @param landmarks_years Landmark times in years.
#' @return A data frame with one row per cut-off.
#' @export
piecewise_table <- function(ipd, endpoint, arm, cutoffs_weeks = c(23, 33, 43), s_gpm, grid,
                            dist = "exp", landmarks_years = c(5, 10, 20)) {
  do.call(rbind, lapply(cutoffs_weeks, function(w) {
    pw <- fit_piecewise(ipd, weeks_to_months(w), dist = dist)
    s <- predict_piecewise(pw, grid)
    lm_idx <- match(landmarks_years * 12, grid)
    out <- data.frame(endpoint = endpoint, arm = arm, cutoff_weeks = w,
                      tail_dist = dist, n_at_cutoff = pw$n_tail, events_after = pw$events_tail,
                      s_at_cutoff = pw$s_cutoff)
    for (k in seq_along(landmarks_years)) out[[sprintf("s_%gy", landmarks_years[k])]] <- s[lm_idx[k]]
    out$rmst_horizon <- rmst_trapezoid(s, grid)
    if (endpoint == "os") {
      sc <- apply_survival_constraints(s, s_gpm = s_gpm)$os
      for (k in seq_along(landmarks_years)) {
        out[[sprintf("s_%gy_constrained", landmarks_years[k])]] <- sc[lm_idx[k]]
      }
      out$rmst_horizon_constrained <- rmst_trapezoid(sc, grid)
    }
    out
  }))
}

#' Proportional-hazards diagnostics for one endpoint
#'
#' @param ipd Stacked IPD with `arm`.
#' @return A data frame with the Grambsch–Therneau (Schoenfeld residual)
#'   test statistic and p-value for the treatment term.
#' @export
ph_test <- function(ipd) {
  ipd$arm <- factor(ipd$arm, levels = c("chemo", "pembro"))
  z <- survival::cox.zph(survival::coxph(survival::Surv(time, status) ~ arm, data = ipd))
  data.frame(chisq = z$table["arm", "chisq"], df = z$table["arm", "df"],
             p = z$table["arm", "p"])
}
