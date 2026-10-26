#' Reconstruct individual patient data from a digitised Kaplan-Meier curve
#'
#' Implements the Guyot et al. (2012) algorithm using
#' [IPDfromKM::getIPD()]. If that call fails, the function falls back to
#' [survHE::digitise()], an independent implementation of the same
#' algorithm (D-008).
#'
#' @param curve Cleaned curve: a data frame with `time_months` and `surv`
#'   (from [clean_digitised_curve()]).
#' @param nar Numbers-at-risk data frame with `time_months` and one column
#'   per arm.
#' @param arm Name of the arm column in `nar` (`"pembro"` or `"chemo"`).
#' @param n_total Number randomised to the arm.
#' @param tot_events Total number of events in the arm at the same data cut,
#'   or `NULL` if not reported (D-004).
#' @return A data frame with columns `time` (months), `status` (1 = event,
#'   0 = censored) and `arm`, with attribute `method` recording the
#'   implementation used.
#' @references Guyot P, Ades AE, Ouwens MJ, Welton NJ. Enhanced secondary
#'   analysis of survival data: reconstructing the data from published
#'   Kaplan-Meier survival curves. BMC Med Res Methodol 2012;12:9.
#' @export
reconstruct_ipd <- function(curve, nar, arm, n_total, tot_events = NULL) {
  if (nar[[arm]][1] != n_total) {
    stop(sprintf("Numbers at risk at time 0 (%d) do not match n_total (%d) for arm '%s'",
                 nar[[arm]][1], n_total, arm))
  }
  ipd <- tryCatch(
    reconstruct_ipd_ipdfromkm(curve, nar, arm, n_total, tot_events),
    error = function(e) {
      warning(sprintf("IPDfromKM failed for arm '%s' (%s); falling back to survHE::digitise()",
                      arm, conditionMessage(e)))
      reconstruct_ipd_survhe(curve, nar, arm, n_total, tot_events)
    }
  )
  ipd$arm <- arm
  ipd
}

#' Build the input object for `IPDfromKM::getIPD()` from a cleaned curve
#'
#' `IPDfromKM::preprocess()` applies its own cleaning: an IQR-based outlier
#' removal and a second monotonicity pass, neither of them logged. That
#' would break D-002 (drop only points beyond the axis, and log every drop)
#' and D-006 (apply one monotonicity method). This function therefore does
#' only the indexing part of `preprocess()`: it assigns each point to a
#' numbers-at-risk interval and builds the interval table. The cleaning
#' itself is left to [clean_digitised_curve()] (D-006a).
#'
#' @param curve Cleaned curve (`time_months`, `surv`), starting at (0, 1).
#' @param trisk,nrisk Numbers-at-risk times and counts for the arm.
#' @return An object of class `preKM`, laid out as `IPDfromKM::preprocess()`
#'   returns it.
#' @export
build_prekm <- function(curve, trisk, nrisk) {
  stopifnot(curve$time_months[1] == 0, curve$surv[1] == 1,
            !is.unsorted(curve$time_months), !is.unsorted(rev(curve$surv)))
  # Drop trailing zero counts after the first zero, as preprocess() does.
  nint <- length(nrisk)
  while (nint > 1 && nrisk[nint] == 0 && nrisk[nint - 1] == 0) nint <- nint - 1
  trisk <- trisk[seq_len(nint)]
  nrisk <- nrisk[seq_len(nint)]

  dat <- data.frame(time = curve$time_months, surv = curve$surv)
  dat$id <- seq_len(nrow(dat))
  dat$interval <- findInterval(dat$time, trisk)
  ints <- sort(unique(dat$interval))
  riskmat <- data.frame(
    interval = ints,
    lower = vapply(ints, function(i) min(dat$id[dat$interval == i]), 1L),
    upper = vapply(ints, function(i) max(dat$id[dat$interval == i]), 1L),
    t.risk = trisk[ints],
    n.risk = nrisk[ints]
  )
  endpts <- if (max(riskmat$t.risk) < max(trisk)) min(nrisk) else NULL
  structure(list(preprocessdat = dat, intervalIndex = riskmat, endpts = endpts,
                 inputdat = dat[, c("time", "surv")]),
            class = "preKM")
}

#' @noRd
reconstruct_ipd_ipdfromkm <- function(curve, nar, arm, n_total, tot_events) {
  prep <- build_prekm(curve, nar$time_months, nar[[arm]])
  est <- IPDfromKM::getIPD(prep = prep, armID = 1, tot.events = tot_events)
  ipd <- data.frame(time = est$IPD$time, status = est$IPD$status)
  attr(ipd, "method") <- "IPDfromKM::getIPD"
  ipd
}

#' @noRd
reconstruct_ipd_survhe <- function(curve, nar, arm, n_total, tot_events) {
  if (!requireNamespace("survHE", quietly = TRUE)) {
    stop("survHE is not installed, so the fallback reconstruction is unavailable")
  }
  dir <- tempfile("survhe_")
  dir.create(dir)
  km_file <- file.path(dir, "km.txt")
  nar_file <- file.path(dir, "nar.txt")
  ipd_file <- file.path(dir, "ipd.txt")
  ki_file <- file.path(dir, "ki.txt")

  km <- data.frame(k = seq_len(nrow(curve)), time = curve$time_months, surv = curve$surv)
  utils::write.table(km, km_file, sep = "\t", row.names = FALSE)

  # survHE expects, for each NaR interval, the first and last index of the
  # digitised points falling in that interval.
  t_nar <- nar$time_months
  lower <- vapply(seq_along(t_nar), function(i) min(which(curve$time_months >= t_nar[i])), 1L)
  upper <- c(lower[-1] - 1L, nrow(curve))
  keep <- lower <= upper
  nar_tab <- data.frame(interval = seq_len(sum(keep)), time = t_nar[keep],
                        lower = lower[keep], upper = upper[keep], nrisk = nar[[arm]][keep])
  utils::write.table(nar_tab, nar_file, sep = "\t", row.names = FALSE)

  survHE::digitise(km_file, nar_file, km_output = ki_file, ipd_output = ipd_file,
                   tot.events = if (is.null(tot_events)) "NA" else tot_events)
  out <- utils::read.table(ipd_file, header = TRUE)
  ipd <- data.frame(time = out[[1]], status = out[[2]])
  attr(ipd, "method") <- "survHE::digitise"
  ipd
}
