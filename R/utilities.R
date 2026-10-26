#' General-population EQ-5D for a cohort of given sex mix
#'
#' Linear interpolation between single years of age in the norms table,
#' weighting females and males by the baseline sex mix. Ages above the
#' table use the last row.
#'
#' @param norms Data frame with `age`, `female`, `male`.
#' @param age Age(s) in years.
#' @param prop_male Proportion male.
#' @return Expected EQ-5D at `age`.
#' @export
norm_utility <- function(norms, age, prop_male) {
  u <- prop_male * norms$male + (1 - prop_male) * norms$female
  stats::approx(norms$age, u, xout = age, rule = 2)$y
}

#' Survival at arbitrary times by log-linear interpolation on a grid
#'
#' @param s Survival on grid `t`.
#' @param t Grid (months).
#' @param x Times to evaluate (months), within the grid.
#' @return Survival at `x`.
#' @noRd
interp_surv <- function(s, t, x) {
  exp(stats::approx(t, log(pmax(s, 1e-300)), xout = x, rule = 2)$y)
}

#' Time-to-death utility per cycle
#'
#' For people alive at time `x`, the probability that death falls in each
#' time-to-death band follows from the OS curve: for a band from `a` to `b`
#' days, `(S(x + a) - S(x + b)) / S(x)`. Utility is the probability-weighted
#' mean of the band utilities (TA531 approach; D-019).
#'
#' @param s_os OS on the extended grid `t` (must reach at least 360 days
#'   past the last evaluation time).
#' @param t Extended grid (months).
#' @param x Evaluation times (months), usually cycle midpoints.
#' @param u Named utilities `u_ttd_0` (<30 days), `u_ttd_30` (30–180),
#'   `u_ttd_180` (180–360) and `u_ttd_360` (360+).
#' @return Utility at each `x`.
#' @export
ttd_utility <- function(s_os, t, x, u) {
  d2m <- function(days) days / (365.25 / 12)
  stopifnot(max(x) + d2m(360) <= max(t) + 1e-9)
  s0 <- interp_surv(s_os, t, x)
  s30 <- interp_surv(s_os, t, x + d2m(30))
  s180 <- interp_surv(s_os, t, x + d2m(180))
  s360 <- interp_surv(s_os, t, x + d2m(360))
  p <- cbind((s0 - s30), (s30 - s180), (s180 - s360), s360) / s0
  p[!is.finite(p)] <- 0
  drop(p %*% c(u[["u_ttd_0"]], u[["u_ttd_30"]], u[["u_ttd_180"]], u[["u_ttd_360"]]))
}
