#' Read a digitised Kaplan-Meier curve
#'
#' Reads a WebPlotDigitizer export with columns `time_months` and `surv`
#' (survival as a proportion), as specified in `data-raw/digitised/README.md`.
#'
#' @param path Path to the CSV file.
#' @return A data frame with numeric columns `time_months` and `surv`.
#' @export
read_digitised_curve <- function(path) {
  dat <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c("time_months", "surv")
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0) {
    stop(sprintf("%s is missing column(s): %s", path, paste(missing, collapse = ", ")))
  }
  dat <- dat[, required]
  if (!is.numeric(dat$time_months) || !is.numeric(dat$surv)) {
    stop(sprintf("%s: time_months and surv must be numeric", path))
  }
  if (anyNA(dat)) stop(sprintf("%s contains missing values", path))
  if (max(dat$surv) > 1.5) {
    stop(sprintf("%s: surv looks like a percentage; export as a proportion (0-1)", path))
  }
  dat
}

#' Read a digitised numbers-at-risk table
#'
#' @param path Path to the CSV file with columns `time_months`, `pembro`, `chemo`.
#' @return A data frame with those three columns, sorted by time.
#' @export
read_nar <- function(path) {
  dat <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c("time_months", "pembro", "chemo")
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0) {
    stop(sprintf("%s is missing column(s): %s", path, paste(missing, collapse = ", ")))
  }
  dat <- dat[order(dat$time_months), required]
  if (dat$time_months[1] != 0) stop(sprintf("%s: first row must be time 0", path))
  for (arm in c("pembro", "chemo")) {
    if (any(diff(dat[[arm]]) > 0)) {
      stop(sprintf("%s: numbers at risk increase over time in arm '%s'", path, arm))
    }
  }
  dat
}

#' Clean a digitised Kaplan-Meier curve
#'
#' Applies the cleaning rules in `DECISIONS.md` (D-002, D-006), in this
#' order, and logs every change:
#'
#' 1. drop points with time beyond the plotted x-axis (`x_max`);
#' 2. set negative times to 0;
#' 3. clip survival to `[0, 1]`;
#' 4. sort by time (ties keep the higher survival first, so steps stay
#'    vertical);
#' 5. enforce non-increasing survival with a running minimum;
#' 6. remove exact duplicate points (no information is lost);
#' 7. add the point (0, 1) at the start if it is not already there.
#'
#' Points beyond the last numbers-at-risk time are kept (D-002).
#'
#' @param dat Data frame from [read_digitised_curve()].
#' @param x_max Last labelled x-axis value of the panel, in months.
#' @param label Curve label used in the log, e.g. `"os_pembro"`.
#' @return A list with `curve` (the cleaned data frame) and `log` (a data
#'   frame with one row per change: `curve`, `rule`, `n_points`, `detail`).
#' @export
clean_digitised_curve <- function(dat, x_max, label) {
  log <- data.frame(curve = character(), rule = character(),
                    n_points = integer(), detail = character())
  add_log <- function(rule, n, detail) {
    if (n > 0) {
      log[nrow(log) + 1, ] <<- list(label, rule, as.integer(n), detail)
    }
  }

  beyond <- dat$time_months > x_max
  add_log("beyond_axis", sum(beyond),
          paste(sprintf("t=%.3f", dat$time_months[beyond]), collapse = "; "))
  dat <- dat[!beyond, ]

  neg_t <- dat$time_months < 0
  add_log("negative_time_to_zero", sum(neg_t),
          paste(sprintf("t=%.3f", dat$time_months[neg_t]), collapse = "; "))
  dat$time_months[neg_t] <- 0

  out_range <- dat$surv < 0 | dat$surv > 1
  add_log("surv_clipped_to_unit_interval", sum(out_range),
          paste(sprintf("t=%.3f s=%.4f", dat$time_months[out_range],
                        dat$surv[out_range]), collapse = "; "))
  dat$surv <- pmin(pmax(dat$surv, 0), 1)

  ord <- order(dat$time_months, -dat$surv)
  add_log("reordered_by_time", sum(ord != seq_along(ord)), "")
  dat <- dat[ord, ]

  s_mono <- cummin(dat$surv)
  changed <- s_mono < dat$surv
  add_log("running_minimum", sum(changed),
          paste(sprintf("t=%.3f %.4f->%.4f", dat$time_months[changed],
                        dat$surv[changed], s_mono[changed]), collapse = "; "))
  dat$surv <- s_mono

  dup <- duplicated(dat[, c("time_months", "surv")])
  add_log("exact_duplicates_removed", sum(dup),
          paste(sprintf("t=%.3f s=%.4f", dat$time_months[dup], dat$surv[dup]), collapse = "; "))
  dat <- dat[!dup, ]

  has_origin <- any(dat$time_months == 0 & dat$surv == 1)
  if (!has_origin) {
    dat <- rbind(data.frame(time_months = 0, surv = 1), dat)
    add_log("origin_added", 1L, "(0, 1)")
  }

  rownames(dat) <- NULL
  list(curve = dat, log = log)
}
