#' Kaplan–Meier step data for plotting
#'
#' @param ipd IPD for one arm.
#' @return Data frame with `time` and `surv`, starting at (0, 1).
#' @noRd
km_steps <- function(ipd) {
  km <- survival::survfit(survival::Surv(time, status) ~ 1, data = ipd)
  data.frame(time = c(0, km$time), surv = c(1, km$surv))
}

#' Plot extrapolated survival for every candidate model
#'
#' One panel per model, each with the KM curve and the model's
#' extrapolation to the horizon. For OS, the curve with the
#' general-population mortality floor is shown dashed.
#'
#' @param fits Output of [fit_survival_models()].
#' @param ipd IPD for the arm.
#' @param grid Extrapolation grid (months).
#' @param title Plot title.
#' @param s_gpm General-population survival on `grid`, or `NULL`.
#' @return A ggplot object.
#' @export
plot_extrapolations <- function(fits, ipd, grid, title, s_gpm = NULL) {
  ok <- names(fits)[!vapply(fits, is.null, logical(1))]
  curves <- do.call(rbind, lapply(ok, function(m) {
    s <- predict_fit(fits[[m]], grid)
    d <- data.frame(model = m, time = grid, surv = s, curve = "Unconstrained")
    if (!is.null(s_gpm)) {
      d <- rbind(d, data.frame(model = m, time = grid, curve = "With general-population floor",
                               surv = apply_survival_constraints(s, s_gpm = s_gpm)$os))
    }
    d
  }))
  curves$model <- factor(curves$model, levels = ok)
  km <- km_steps(ipd)
  km <- km[rep(seq_len(nrow(km)), length(ok)), ]
  km$model <- factor(rep(ok, each = nrow(km) / length(ok)), levels = ok)
  ggplot2::ggplot(curves, ggplot2::aes(time / 12, surv)) +
    ggplot2::geom_step(data = km, ggplot2::aes(time / 12, surv), colour = "grey40",
                       linewidth = 0.4) +
    ggplot2::geom_line(ggplot2::aes(linetype = curve), colour = "#1b6ca8", linewidth = 0.5) +
    ggplot2::facet_wrap(~model, ncol = 3) +
    ggplot2::scale_linetype_manual(values = c("Unconstrained" = "solid",
                                              "With general-population floor" = "22"),
                                   name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(x * 100, "%")) +
    ggplot2::labs(x = "Years since randomisation", y = "Survival", title = title,
                  caption = "Grey: Kaplan–Meier of reconstructed IPD.") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = if (is.null(s_gpm)) "none" else "top")
}

#' Plot smoothed empirical hazard against model hazards
#'
#' The empirical hazard is a kernel-smoothed estimate (`muhaz`) over the
#' observed follow-up, up to the time when 10% of the arm remains at risk.
#' Model hazards are drawn over the same range.
#'
#' @param fits Output of [fit_survival_models()].
#' @param ipd IPD for the arm.
#' @param title Plot title.
#' @return A ggplot object.
#' @export
plot_hazards <- function(fits, ipd, title) {
  n <- nrow(ipd)
  t_sorted <- sort(ipd$time, decreasing = TRUE)
  t_max <- t_sorted[ceiling(0.1 * n)]
  mh <- muhaz::muhaz(ipd$time, ipd$status, min.time = 0, max.time = t_max)
  emp <- data.frame(time = mh$est.grid, haz = mh$haz.est)
  ok <- names(fits)[!vapply(fits, is.null, logical(1))]
  tt <- seq(0.25, t_max, length.out = 200)
  mod <- do.call(rbind, lapply(ok, function(m) {
    data.frame(model = m, time = tt, haz = predict_fit(fits[[m]], tt, "hazard"))
  }))
  mod$model <- factor(mod$model, levels = ok)
  emp <- emp[rep(seq_len(nrow(emp)), length(ok)), ]
  emp$model <- factor(rep(ok, each = nrow(emp) / length(ok)), levels = ok)
  ggplot2::ggplot(mod, ggplot2::aes(time, haz)) +
    ggplot2::geom_line(data = emp, colour = "grey40", linewidth = 0.5) +
    ggplot2::geom_line(colour = "#c8553d", linewidth = 0.5) +
    ggplot2::facet_wrap(~model, ncol = 3) +
    ggplot2::coord_cartesian(ylim = c(0, max(emp$haz) * 1.5)) +
    ggplot2::labs(x = "Months since randomisation", y = "Hazard (per month)", title = title,
                  caption = paste("Grey: kernel-smoothed hazard of reconstructed IPD (muhaz),",
                                  "to the time 10% remain at risk. Red: model hazard.")) +
    ggplot2::theme_minimal(base_size = 10)
}

#' Log-cumulative hazard plot for both arms
#'
#' Plots `log(-log(S_KM))` against `log(time)`. Parallel lines support
#' proportional hazards. Straight lines support a Weibull model, with a
#' slope of 1 for an exponential.
#'
#' @param ipd Stacked IPD with `arm`.
#' @param title Plot title.
#' @return A ggplot object.
#' @export
plot_log_cumhaz <- function(ipd, title) {
  d <- do.call(rbind, lapply(c("pembro", "chemo"), function(a) {
    k <- km_steps(ipd[ipd$arm == a, ])
    k <- k[k$time > 0 & k$surv > 0 & k$surv < 1, ]
    data.frame(arm = a, log_t = log(k$time), lch = log(-log(k$surv)))
  }))
  ggplot2::ggplot(d, ggplot2::aes(log_t, lch, colour = arm)) +
    ggplot2::geom_step(linewidth = 0.5) +
    ggplot2::scale_colour_manual(values = c(pembro = "#1b6ca8", chemo = "#c8553d"),
                                 labels = c(pembro = "Pembrolizumab", chemo = "Chemotherapy"),
                                 name = NULL) +
    ggplot2::labs(x = "log(time in months)", y = "log(cumulative hazard)", title = title) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "top")
}
