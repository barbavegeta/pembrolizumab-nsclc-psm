#' Cost-effectiveness plane
#'
#' @param psa Output of [run_psa()].
#' @param det Deterministic incremental results (`inc_qaly`, `inc_cost`).
#' @param thresholds Willingness-to-pay lines to draw.
#' @return A ggplot object.
#' @export
plot_ce_plane <- function(psa, det, thresholds = c(20000, 30000)) {
  lines <- data.frame(k = thresholds, label = paste0("£", formatC(thresholds, big.mark = ",",
                                                                   format = "d"), "/QALY"))
  ggplot2::ggplot(psa, ggplot2::aes(inc_qaly, inc_cost)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey70") +
    ggplot2::geom_point(alpha = 0.25, size = 0.8, colour = "#1b6ca8") +
    ggplot2::geom_abline(data = lines, ggplot2::aes(slope = k, intercept = 0, linetype = label),
                         colour = "grey30") +
    ggplot2::annotate("point", x = det[["inc_qaly"]], y = det[["inc_cost"]], colour = "#c8553d",
                      size = 3, shape = 18) +
    ggplot2::scale_y_continuous(labels = function(x) paste0("£", formatC(x, big.mark = ",", format = "d"))) +
    ggplot2::labs(x = "Incremental QALYs", y = "Incremental cost", linetype = NULL,
                  caption = "Points: PSA iterations. Red diamond: deterministic base case.") +
    ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "top")
}

#' Cost-effectiveness acceptability curve plot
#'
#' @param ce Output of [ceac()].
#' @return A ggplot object.
#' @export
plot_ceac <- function(ce) {
  ggplot2::ggplot(ce, ggplot2::aes(wtp, p_ce)) +
    ggplot2::geom_line(colour = "#1b6ca8", linewidth = 0.7) +
    ggplot2::geom_vline(xintercept = c(20000, 30000), linetype = "22", colour = "grey50") +
    ggplot2::scale_x_continuous(labels = function(x) paste0("£", x / 1000, "k")) +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(x * 100, "%")) +
    ggplot2::labs(x = "Willingness to pay per QALY",
                  y = "Probability pembrolizumab is cost-effective") +
    ggplot2::theme_minimal(base_size = 11)
}

#' Tornado plot of the one-way sensitivity analysis
#'
#' @param owsa Output of [run_owsa()].
#' @param base_icer Base-case ICER.
#' @param n_show Number of most influential inputs to show.
#' @return A ggplot object.
#' @export
plot_tornado <- function(owsa, base_icer, n_show = 15) {
  owsa$range <- abs(owsa$icer_high - owsa$icer_low)
  d <- utils::head(owsa[order(-owsa$range), ], n_show)
  d$label <- factor(d$description, levels = rev(d$description))
  long <- rbind(data.frame(label = d$label, end = "Low input value", icer = d$icer_low),
                data.frame(label = d$label, end = "High input value", icer = d$icer_high))
  ggplot2::ggplot(long) +
    ggplot2::geom_segment(ggplot2::aes(x = base_icer, xend = icer, y = label, yend = label,
                                       colour = end), linewidth = 4) +
    ggplot2::geom_vline(xintercept = base_icer, colour = "grey30") +
    ggplot2::scale_colour_manual(values = c("Low input value" = "#1b6ca8",
                                            "High input value" = "#c8553d"), name = NULL) +
    ggplot2::scale_x_continuous(labels = function(x) paste0("£", formatC(x, big.mark = ",", format = "d"))) +
    ggplot2::labs(x = "ICER (£ per QALY)", y = NULL) +
    ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(legend.position = "top")
}

#' Validation against TA531
#'
#' Compares modelled undiscounted mean survival with the 36.0 and 22.3
#' months reported in TA531 §3.18, for the base case and for the TA531-like
#' extrapolation (KM to 33 weeks, then exponential).
#'
#' @param base Base-case output of [summarise_results()].
#' @param scen Output of [run_scenarios()].
#' @return A data frame.
#' @export
ta531_validation <- function(base, scen) {
  ta531 <- c(pembro = 36.0, chemo = 22.3)
  pw <- scen[scen$scenario == "OS: KM to 33 weeks + exponential (TA531)", ]
  data.frame(
    analysis = c("TA531 company model (TA531 §3.18, July 2017 data)",
                 "This model, base case (lognormal, 2020 data)",
                 "This model, TA531-like extrapolation (KM to 33 weeks + exponential, 2020 data)"),
    mean_os_months_pembro = c(ta531[["pembro"]], 12 * base$totals["pembro", "ly_undisc"],
                              12 * pw$ly_undisc_pembro),
    mean_os_months_chemo = c(ta531[["chemo"]], 12 * base$totals["chemo", "ly_undisc"],
                             12 * pw$ly_undisc_chemo)
  )
}

#' Base-case modelled curves against the Kaplan–Meier data
#'
#' @param curves Output of [build_curves()].
#' @param t Grid (months) on which `curves` are defined.
#' @param ipd Named list (`os`, `pfs`) of stacked IPD.
#' @param horizon_months Plot limit.
#' @return A ggplot object.
#' @export
plot_base_curves <- function(curves, t, ipd, horizon_months = 240) {
  keep <- t <= horizon_months
  mod <- do.call(rbind, lapply(c("pembro", "chemo"), function(a) {
    rbind(data.frame(arm = a, endpoint = "Overall survival", time = t[keep], surv = curves[[a]]$os[keep]),
          data.frame(arm = a, endpoint = "Progression-free survival", time = t[keep],
                     surv = curves[[a]]$pfs[keep]))
  }))
  km <- do.call(rbind, lapply(c("os", "pfs"), function(ep) {
    do.call(rbind, lapply(c("pembro", "chemo"), function(a) {
      k <- km_steps(ipd[[ep]][ipd[[ep]]$arm == a, ])
      data.frame(arm = a, endpoint = if (ep == "os") "Overall survival" else "Progression-free survival",
                 time = k$time, surv = k$surv)
    }))
  }))
  cols <- c(pembro = "#1b6ca8", chemo = "#c8553d")
  labs <- c(pembro = "Pembrolizumab", chemo = "Chemotherapy")
  ggplot2::ggplot(mod, ggplot2::aes(time / 12, surv, colour = arm)) +
    ggplot2::geom_step(data = km, linewidth = 0.4, alpha = 0.7) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::facet_wrap(~endpoint) +
    ggplot2::scale_colour_manual(values = cols, labels = labs, name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(x * 100, "%")) +
    ggplot2::labs(x = "Years since randomisation", y = "Proportion",
                  caption = "Steps: Kaplan–Meier of reconstructed IPD. Lines: base-case model (lognormal, constrained).") +
    ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "top")
}
