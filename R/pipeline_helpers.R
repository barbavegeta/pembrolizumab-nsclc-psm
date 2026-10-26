#' Plotted x-axis limit for an endpoint
#'
#' @param axis_limits Data frame read from `axis_limits.csv`.
#' @param endpoint `"os"` or `"pfs"`.
#' @return The last labelled x-axis value in months.
#' @export
axis_max <- function(axis_limits, endpoint) {
  x <- axis_limits$x_max_months[axis_limits$endpoint == endpoint]
  if (length(x) != 1) stop(sprintf("axis_limits.csv needs exactly one row for '%s'", endpoint))
  x
}

#' Reconstruct one arm using the documented study inputs
#'
#' Looks up the number randomised and the total events in
#' `data-raw/reconstruction_inputs.csv` and calls [reconstruct_ipd()].
#'
#' @param curve Cleaned curve.
#' @param nar Numbers-at-risk table.
#' @param arm `"pembro"` or `"chemo"`.
#' @param recon_inputs Data frame read from `reconstruction_inputs.csv`.
#' @param endpoint `"os"` or `"pfs"`.
#' @return Reconstructed IPD with columns `time`, `status`, `arm`, `method`.
#' @export
reconstruct_from_inputs <- function(curve, nar, arm, recon_inputs, endpoint) {
  row <- recon_inputs[recon_inputs$endpoint == endpoint & recon_inputs$arm == arm, ]
  if (nrow(row) != 1) stop(sprintf("No unique reconstruction input row for %s/%s", endpoint, arm))
  ev <- if (is.na(row$total_events)) NULL else row$total_events
  ipd <- reconstruct_ipd(curve, nar, arm, n_total = row$n_randomised, tot_events = ev)
  ipd$method <- attr(ipd, "method")
  ipd
}

#' Summarise the implementation used for each endpoint and arm
#'
#' @param ipd_list Named list of stacked IPD data frames.
#' @return Named character vector, e.g. `c("os / pembro" = "IPDfromKM::getIPD")`.
#' @export
describe_methods <- function(ipd_list) {
  out <- unlist(lapply(names(ipd_list), function(ep) {
    d <- unique(ipd_list[[ep]][, c("arm", "method")])
    stats::setNames(d$method, paste(ep, "/", d$arm))
  }))
  out
}

#' Stack two cleaned curves with an arm label
#'
#' @param pembro,chemo Outputs of [clean_digitised_curve()].
#' @return Data frame with `time_months`, `surv`, `arm`.
#' @export
stack_curves <- function(pembro, chemo) {
  rbind(cbind(pembro$curve, arm = "pembro"), cbind(chemo$curve, arm = "chemo"))
}

#' Save a ggplot to a file and return the path
#'
#' @param plot A ggplot object.
#' @param path Output path (PNG).
#' @param width,height Size in inches.
#' @return `path`.
#' @export
save_plot <- function(plot, path, width = 7, height = 4.5) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = 150, bg = "white")
  path
}

#' Write reconstructed IPD once the reconstruction has been accepted
#'
#' @param ipd Reconstructed IPD.
#' @param path Output CSV path.
#' @param accepted Output of [assert_reconstruction()]. Taking it as an
#'   argument makes this target depend on acceptance.
#' @return `path`.
#' @export
write_ipd <- function(ipd, path, accepted) {
  stopifnot(isTRUE(accepted))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(ipd[, c("arm", "time", "status")], path, row.names = FALSE)
  path
}

#' Plot title for a curve identifier
#'
#' @param id Curve identifier, e.g. `"os_pembro"`.
#' @return A title such as `"Overall survival: Pembrolizumab"`.
#' @export
curve_title <- function(id) {
  ep <- c(os = "Overall survival", pfs = "Progression-free survival")[[sub("_.*", "", id)]]
  arm <- c(pembro = "Pembrolizumab", chemo = "Chemotherapy")[[sub(".*_", "", id)]]
  paste0(ep, ": ", arm)
}

#' Index of survival figures for the report
#'
#' @param extrap,haz Paths of the per-curve extrapolation and hazard plots.
#' @param lch Paths of the log-cumulative hazard plots.
#' @return Character vector of file names relative to `output/survival/`,
#'   named by file stem (e.g. `"extrap_os_pembro"`, `"lch_os"`).
#' @export
survival_figure_index <- function(extrap, haz, lch) {
  files <- basename(c(extrap, haz, lch))
  stats::setNames(files, sub("\\.png$", "", files))
}

#' Write the per-cycle traces of a model run
#'
#' @param run Output of [run_model()].
#' @param dir Output directory.
#' @return Paths of the CSV files written.
#' @export
write_traces <- function(run, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  vapply(c("pembro", "chemo"), function(a) {
    f <- file.path(dir, sprintf("trace_%s.csv", a))
    utils::write.csv(run[[a]]$trace, f, row.names = FALSE)
    f
  }, character(1), USE.NAMES = FALSE)
}

#' Write a CSV and return its path
#'
#' @param x Data frame.
#' @param path Output path.
#' @param row.names Passed to [utils::write.csv()].
#' @return `path`.
#' @export
save_csv <- function(x, path, row.names = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = row.names)
  path
}

#' Model diagnostics quoted in the report
#'
#' @param curves Base-case curves ([build_curves()]).
#' @param p Parameters.
#' @param time Output of [model_time()].
#' @param norms EQ-5D norms.
#' @param fits Fits by curve.
#' @param owsa Output of [run_owsa()].
#' @return A list: `mean_chemo_cycles`, `mean_pembro_doses`, `norm_at_start`,
#'   `gengamma_q` (data frame of Q estimates and 95% CIs by curve) and
#'   `owsa_top` (descriptions of the two most influential inputs) and
#'   `acceptance` ([acceptance_thresholds()]).
#' @export
model_diagnostics <- function(curves, p, time, norms, fits, owsa) {
  idx <- seq_len(time$n_cycles + 1)
  q <- do.call(rbind, lapply(names(fits), function(k) {
    g <- fits[[k]][["Generalised gamma"]]
    data.frame(curve = k, est = g$res["Q", "est"], lower = g$res["Q", "L95%"], upper = g$res["Q", "U95%"])
  }))
  o <- owsa[order(-abs(owsa$icer_high - owsa$icer_low)), ]
  list(
    mean_chemo_cycles = sum(dosing_stream(curves$chemo$pfs[idx], time$t_weeks, 0,
                                          p$chemo_interval_weeks, p$chemo_max_cycles)),
    mean_pembro_doses = sum(dosing_stream(curves$pembro$pfs[idx], time$t_weeks, 0,
                                          p$pembro_interval_weeks, p$pembro_max_doses)),
    norm_at_start = norm_utility(norms, p$start_age, p$prop_male),
    gengamma_q = q,
    owsa_top = utils::head(o$description, 2),
    acceptance = acceptance_thresholds()
  )
}
