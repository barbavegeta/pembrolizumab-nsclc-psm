#' Candidate parametric survival models
#'
#' The model set follows NICE DSU TSD 14 (six standard distributions) and
#' TSD 21 (Royston–Parmar splines on the hazard scale with 1–3 internal
#' knots; D-014).
#'
#' @return A data frame with `model` (label), `type` (`"flexsurvreg"` or
#'   `"spline"`), `dist` (flexsurv distribution name, or `NA`) and `knots`
#'   (number of internal knots, or `NA`).
#' @export
candidate_models <- function() {
  data.frame(
    model = c("Exponential", "Weibull", "Gompertz", "Log-logistic", "Lognormal",
              "Generalised gamma", "RP spline (1 knot)", "RP spline (2 knots)",
              "RP spline (3 knots)"),
    type = c(rep("flexsurvreg", 6), rep("spline", 3)),
    dist = c("exp", "weibull", "gompertz", "llogis", "lnorm", "gengamma", NA, NA, NA),
    knots = c(rep(NA, 6), 1:3)
  )
}

#' Fit all candidate models to one arm
#'
#' @param ipd Data frame with `time` (months) and `status`, for one arm.
#' @param models Output of [candidate_models()].
#' @return A named list of fitted `flexsurvreg` objects, named by model.
#'   A model that fails to converge is `NULL`. The list's attribute
#'   `"errors"` is a named character vector of error messages, with `NA`
#'   for models that fitted.
#' @export
fit_survival_models <- function(ipd, models = candidate_models()) {
  fits <- lapply(seq_len(nrow(models)), function(i) {
    m <- models[i, ]
    tryCatch({
      if (m$type == "flexsurvreg") {
        flexsurv::flexsurvreg(survival::Surv(time, status) ~ 1, data = ipd, dist = m$dist)
      } else {
        flexsurv::flexsurvspline(survival::Surv(time, status) ~ 1, data = ipd,
                                 k = m$knots, scale = "hazard")
      }
    }, error = function(e) structure(list(), error = conditionMessage(e)))
  })
  errors <- vapply(fits, function(f) {
    e <- attr(f, "error")
    if (is.null(e)) NA_character_ else e
  }, character(1))
  fits <- lapply(fits, function(f) if (inherits(f, "flexsurvreg")) f else NULL)
  structure(stats::setNames(fits, models$model),
            errors = stats::setNames(errors, models$model))
}

#' Survival, hazard and cumulative hazard from a fitted model
#'
#' @param fit A `flexsurvreg` object.
#' @param t Numeric vector of times (months).
#' @param type `"survival"`, `"hazard"` or `"cumhaz"`.
#' @return Numeric vector the same length as `t`.
#' @export
predict_fit <- function(fit, t, type = c("survival", "hazard", "cumhaz")) {
  type <- match.arg(type)
  out <- summary(fit, t = t, type = type, ci = FALSE, tidy = TRUE)
  out$est[match(t, out$time)]
}

#' Restricted mean survival time by numerical integration
#'
#' @param surv Survival values on the grid `t`.
#' @param t Grid of times starting at 0.
#' @return Area under the survival curve up to `max(t)` (trapezium rule).
#' @export
rmst_trapezoid <- function(surv, t) {
  sum(diff(t) * (utils::head(surv, -1) + utils::tail(surv, -1)) / 2)
}
