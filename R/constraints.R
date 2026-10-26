#' General-population survival for a cohort of given age and sex mix
#'
#' Converts life-table `qx` to a constant annual hazard within each year of
#' age, `-log(1 - qx)`, and computes survival from the starting age
#' separately for males and females. Cohort survival is the sex-weighted
#' mixture of the two curves. This lets the sex mix of survivors shift over
#' time, rather than holding the hazard weights fixed. Ages beyond the last
#' row of the life table keep the hazard of that last row.
#'
#' @param lifetable Data frame with `age`, `sex` (`"male"`/`"female"`) and
#'   `qx`.
#' @param start_age Cohort age at time 0 (years).
#' @param prop_male Proportion male at time 0.
#' @param t Times since baseline, in months.
#' @return General-population survival at `t`.
#' @export
gpm_survival <- function(lifetable, start_age, prop_male, t) {
  cumhaz_sex <- function(sex) {
    lt <- lifetable[lifetable$sex == sex, ]
    lt <- lt[order(lt$age), ]
    h <- -log(1 - lt$qx)
    max_age <- max(lt$age)
    age_at <- start_age + t / 12
    vapply(age_at, function(a) {
      # Integrate the piecewise-constant hazard from start_age to a.
      edges <- seq(floor(start_age), max(floor(a), floor(start_age)) + 1)
      lo <- pmax(edges[-length(edges)], start_age)
      hi <- pmin(edges[-1], a)
      dur <- pmax(hi - lo, 0)
      ages <- pmin(edges[-length(edges)], max_age)
      sum(dur * h[match(ages, lt$age)])
    }, numeric(1))
  }
  prop_male * exp(-cumhaz_sex("male")) + (1 - prop_male) * exp(-cumhaz_sex("female"))
}

#' Apply survival constraints on a time grid
#'
#' Applies, interval by interval on the grid:
#'
#' 1. OS hazard is never below general-population mortality. In each
#'    interval the OS cumulative-hazard increment is the larger of the
#'    model increment and the general-population increment.
#' 2. PFS never exceeds OS: `S_pfs = min(S_pfs, S_os)` after step 1.
#'
#' The function is generic: it takes survival values on any grid starting
#' at 0, so the partitioned survival engine can reuse it.
#'
#' @param s_os,s_pfs Unconstrained OS and PFS survival on the grid `t`.
#'   `s_pfs` may be `NULL`.
#' @param s_gpm General-population survival on the same grid.
#' @return A list with constrained `os` and `pfs`, plus logical vectors
#'   `os_gpm_binding` (the general-population hazard was used in the
#'   interval ending at that grid point) and `pfs_capped` (PFS was capped at
#'   OS).
#' @export
apply_survival_constraints <- function(s_os, s_pfs = NULL, s_gpm) {
  stopifnot(length(s_os) == length(s_gpm), s_os[1] == 1, s_gpm[1] == 1)
  dh_model <- diff(-log(pmax(s_os, .Machine$double.xmin)))
  dh_gpm <- diff(-log(s_gpm))
  binding <- dh_gpm > dh_model
  os <- exp(-c(0, cumsum(pmax(dh_model, dh_gpm))))
  out <- list(os = os, os_gpm_binding = c(FALSE, binding))
  if (!is.null(s_pfs)) {
    stopifnot(length(s_pfs) == length(s_os))
    out$pfs <- pmin(s_pfs, os)
    out$pfs_capped <- s_pfs > os
  }
  out
}
