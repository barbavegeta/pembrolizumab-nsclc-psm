#' Convert weeks to months
#'
#' @param weeks Numeric weeks.
#' @return Months, using 365.25 / 12 days per month.
#' @export
weeks_to_months <- function(weeks) weeks * 7 / (365.25 / 12)

#' Fit a piecewise KM + parametric tail model (TA531 approach)
#'
#' Uses the Kaplan–Meier estimate up to `cutoff`. After that, a parametric
#' model is fitted to the patients still at risk at the cut-off, with time
#' measured from the cut-off. Survival after the cut-off is
#' `S_KM(cutoff) * S_tail(t - cutoff)`. TA531 used an exponential tail with
#' cut-offs of 23, 33 and 43 weeks (`docs/source_checks.md`, N4).
#'
#' @param ipd Data frame with `time` (months) and `status`, for one arm.
#' @param cutoff Cut-off in months.
#' @param dist flexsurv distribution for the tail (default `"exp"`).
#' @return An object of class `piecewise_fit`: a list with `km` (the
#'   `survfit` object), `cutoff`, `s_cutoff` (KM survival at the cut-off),
#'   `tail` (the `flexsurvreg` fit) and `n_tail` and `events_tail`.
#' @export
fit_piecewise <- function(ipd, cutoff, dist = "exp") {
  km <- survival::survfit(survival::Surv(time, status) ~ 1, data = ipd)
  s_cut <- summary(km, times = cutoff, extend = TRUE)$surv
  after <- ipd[ipd$time > cutoff, ]
  after$time <- after$time - cutoff
  tail <- flexsurv::flexsurvreg(survival::Surv(time, status) ~ 1, data = after, dist = dist)
  structure(list(km = km, cutoff = cutoff, s_cutoff = s_cut, tail = tail,
                 n_tail = nrow(after), events_tail = sum(after$status)),
            class = "piecewise_fit")
}

#' Survival from a piecewise fit
#'
#' @param pw A `piecewise_fit`.
#' @param t Times in months.
#' @return Survival at `t`.
#' @export
predict_piecewise <- function(pw, t) {
  s <- numeric(length(t))
  pre <- t <= pw$cutoff
  if (any(pre)) s[pre] <- summary(pw$km, times = t[pre], extend = TRUE)$surv
  if (any(!pre)) s[!pre] <- pw$s_cutoff * predict_fit(pw$tail, t[!pre] - pw$cutoff)
  s
}
