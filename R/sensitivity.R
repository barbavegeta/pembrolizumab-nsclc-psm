#' Uncertainty distributions for scalar inputs
#'
#' Reads `psa_dist` and `psa_se` from `inputs.csv`. `psa_se` is either an
#' absolute standard error, a percentage of the mean (e.g. `"20%"`), or a
#' sample size for a proportion (e.g. `"n=151"`).
#'
#' @param inputs Output of [read_inputs()].
#' @return A data frame with `parameter`, `dist` (`"beta"` or `"gamma"`),
#'   `mean` and `se`, for parameters that vary.
#' @export
param_distributions <- function(inputs) {
  x <- inputs[inputs$psa_dist %in% c("beta", "gamma"), ]
  mean <- as.numeric(x$value)
  se <- vapply(seq_len(nrow(x)), function(i) {
    s <- x$psa_se[i]
    if (grepl("%$", s)) return(mean[i] * as.numeric(sub("%", "", s)) / 100)
    if (grepl("^n=", s)) {
      n <- as.numeric(sub("n=", "", s))
      return(sqrt(mean[i] * (1 - mean[i]) / n))
    }
    as.numeric(s)
  }, numeric(1))
  data.frame(parameter = x$parameter, dist = x$psa_dist, mean = mean, se = se)
}

#' Method-of-moments distribution parameters
#'
#' @param dist `"beta"` or `"gamma"`.
#' @param mean,se Mean and standard error.
#' @return For beta, `c(shape1, shape2)`; for gamma, `c(shape, rate)`.
#' @export
mom_params <- function(dist, mean, se) {
  if (dist == "beta") {
    k <- mean * (1 - mean) / se^2 - 1
    c(mean * k, (1 - mean) * k)
  } else {
    c(mean^2 / se^2, mean / se^2)
  }
}

#' Quantile of an input's uncertainty distribution
#'
#' @param dist,mean,se As in [mom_params()].
#' @param q Probability.
#' @return The quantile.
#' @export
dist_quantile <- function(dist, mean, se, q) {
  a <- mom_params(dist, mean, se)
  if (dist == "beta") stats::qbeta(q, a[1], a[2]) else stats::qgamma(q, a[1], a[2])
}

#' Draw PSA parameter sets
#'
#' Scalar inputs: beta for utilities and proportions, gamma for costs.
#' Adverse-event rates and subsequent-therapy proportions: beta from the
#' trial counts. Survival: multivariate normal draws of each chosen model's
#' parameters from its variance–covariance matrix
#' (`flexsurv::normboot.flexsurvreg`).
#'
#' @param p Deterministic parameters ([model_params()]).
#' @param data Output of [load_model_data()].
#' @param fits Fits by curve (named list of [fit_survival_models()] outputs).
#' @param spec Curve specification ([curve_spec()]) giving the models.
#' @param n Number of draws.
#' @param seed Random seed.
#' @return A list with `params` (a list of `n` parameter lists) and
#'   `surv_pars` (a list of `n` named lists of survival parameter vectors).
#' @export
draw_psa <- function(p, data, fits, spec, n, seed) {
  set.seed(seed)
  dists <- param_distributions(data$inputs)
  scalar <- vapply(seq_len(nrow(dists)), function(i) {
    a <- mom_params(dists$dist[i], dists$mean[i], dists$se[i])
    if (dists$dist[i] == "beta") stats::rbeta(n, a[1], a[2]) else stats::rgamma(n, a[1], a[2])
  }, numeric(n))
  scalar <- matrix(scalar, nrow = n, dimnames = list(NULL, dists$parameter))

  ae <- data$adverse_events
  rb <- function(rate, size) {
    if (rate <= 0) return(rep(0, n))
    stats::rbeta(n, rate * size, (1 - rate) * size)
  }
  ae_p <- vapply(seq_len(nrow(ae)), function(i) rb(ae$rate_pembro[i], ae$n_pembro[i]), numeric(n))
  ae_c <- vapply(seq_len(nrow(ae)), function(i) rb(ae$rate_chemo[i], ae$n_chemo[i]), numeric(n))
  sq <- data$subsequent
  sq_p <- vapply(seq_len(nrow(sq)), function(i) stats::rbeta(n, sq$n_pembro[i], sq$N_pembro[i] - sq$n_pembro[i]), numeric(n))
  sq_c <- vapply(seq_len(nrow(sq)), function(i) stats::rbeta(n, sq$n_chemo[i], sq$N_chemo[i] - sq$n_chemo[i]), numeric(n))

  surv <- list()
  for (key in names(fits)) {
    model <- if (startsWith(key, "os")) spec$os_model else spec$pfs_model
    fit <- fits[[key]][[model]]
    d <- flexsurv::normboot.flexsurvreg(fit, B = n, raw = TRUE, transform = FALSE)
    d <- matrix(d, nrow = n, dimnames = list(NULL, rownames(fit$res)))
    surv[[key]] <- d
  }

  params <- lapply(seq_len(n), function(i) {
    q <- p
    for (nm in colnames(scalar)) q[[nm]] <- unname(scalar[i, nm])
    q$ae_rate_pembro[] <- ae_p[i, ]
    q$ae_rate_chemo[] <- ae_c[i, ]
    q$subsq_prop_pembro[] <- sq_p[i, ]
    q$subsq_prop_chemo[] <- sq_c[i, ]
    q
  })
  surv_pars <- lapply(seq_len(n), function(i) lapply(surv, function(m) m[i, ]))
  list(params = params, surv_pars = surv_pars)
}

#' Run the model for one parameter set and return the headline results
#'
#' @param p Parameters.
#' @param settings Model settings.
#' @param ctx A list with `fits`, `ipd`, `t_ext`, `s_gpm`, `time` and `norms`.
#' @param surv_pars Optional survival parameter draws (see [build_curves()]).
#' @return The `incremental` vector from [summarise_results()], plus each
#'   arm's discounted costs and QALYs.
#' @export
run_once <- function(p, settings, ctx, surv_pars = NULL) {
  curves <- build_curves(settings$curves, ctx$fits, ctx$ipd, ctx$t_ext, ctx$s_gpm, surv_pars)
  res <- run_model(p, curves, settings, ctx$time, ctx$t_ext, ctx$norms)
  s <- summarise_results(res)
  c(s$incremental, cost_pembro = s$totals["pembro", "cost_disc"],
    cost_chemo = s$totals["chemo", "cost_disc"], qaly_pembro = s$totals["pembro", "qaly_disc"],
    qaly_chemo = s$totals["chemo", "qaly_disc"], ly_undisc_pembro = s$totals["pembro", "ly_undisc"],
    ly_undisc_chemo = s$totals["chemo", "ly_undisc"])
}

#' One-way sensitivity analysis
#'
#' Each varying input is set to the 2.5th and 97.5th percentiles of its
#' uncertainty distribution in turn. Grouped inputs (all AE unit costs,
#' all subsequent-therapy proportions per arm) are varied together.
#'
#' @param p,settings,ctx As in [run_once()].
#' @param data Output of [load_model_data()].
#' @return A data frame with `parameter`, `description`, `low_value`,
#'   `high_value`, `icer_low`, `icer_high`, `nmb30_low` and `nmb30_high`.
#' @export
run_owsa <- function(p, settings, ctx, data) {
  dists <- param_distributions(data$inputs)
  desc <- stats::setNames(data$inputs$description, data$inputs$parameter)
  cases <- lapply(seq_len(nrow(dists)), function(i) {
    lo <- dist_quantile(dists$dist[i], dists$mean[i], dists$se[i], 0.025)
    hi <- dist_quantile(dists$dist[i], dists$mean[i], dists$se[i], 0.975)
    nm <- dists$parameter[i]
    list(parameter = nm, description = desc[[nm]], low_value = lo, high_value = hi,
         set = function(q, v) { q[[nm]] <- v; q })
  })
  cases[[length(cases) + 1]] <- list(
    parameter = "ae_costs", description = "All adverse-event unit costs (±20%, gamma 95% interval)",
    low_value = dist_quantile("gamma", 1, 0.2, 0.025), high_value = dist_quantile("gamma", 1, 0.2, 0.975),
    set = function(q, v) { q$ae_cost_201415 <- q$ae_cost_201415 * v; q })
  for (a in c("pembro", "chemo")) {
    local({
      arm <- a
      cases[[length(cases) + 1]] <<- list(
        parameter = paste0("subsequent_", arm),
        description = sprintf("Subsequent-therapy use, %s arm (each therapy at its 2.5th/97.5th beta percentile)",
                              if (arm == "pembro") "pembrolizumab" else "chemotherapy"),
        low_value = NA, high_value = NA,
        set = function(q, v) {
          s <- data$subsequent
          n <- if (arm == "pembro") s$n_pembro else s$n_chemo
          N <- if (arm == "pembro") s$N_pembro else s$N_chemo
          prop <- stats::qbeta(v, n, N - n)
          q[[paste0("subsq_prop_", arm)]][] <- prop
          q
        }, quantiles = TRUE)
    })
  }
  do.call(rbind, lapply(cases, function(cs) {
    vlo <- if (isTRUE(cs$quantiles)) 0.025 else cs$low_value
    vhi <- if (isTRUE(cs$quantiles)) 0.975 else cs$high_value
    r_lo <- run_once(cs$set(p, vlo), settings, ctx)
    r_hi <- run_once(cs$set(p, vhi), settings, ctx)
    data.frame(parameter = cs$parameter, description = cs$description,
               low_value = cs$low_value, high_value = cs$high_value,
               icer_low = r_lo[["icer"]], icer_high = r_hi[["icer"]],
               nmb30_low = r_lo[["nmb_30000"]], nmb30_high = r_hi[["nmb_30000"]])
  }))
}

#' Probabilistic sensitivity analysis
#'
#' @param draws Output of [draw_psa()].
#' @param settings,ctx As in [run_once()].
#' @return A data frame with one row per iteration.
#' @export
run_psa <- function(draws, settings, ctx) {
  out <- lapply(seq_along(draws$params), function(i) {
    run_once(draws$params[[i]], settings, ctx, draws$surv_pars[[i]])
  })
  cbind(iteration = seq_along(out), as.data.frame(do.call(rbind, out)))
}

#' Cost-effectiveness acceptability curve
#'
#' @param psa Output of [run_psa()].
#' @param wtp Willingness-to-pay values (GBP per QALY).
#' @return A data frame with `wtp` and `p_ce` (probability pembrolizumab is
#'   cost-effective, i.e. incremental net monetary benefit > 0).
#' @export
ceac <- function(psa, wtp = seq(0, 100000, by = 1000)) {
  data.frame(wtp = wtp, p_ce = vapply(wtp, function(k) mean(k * psa$inc_qaly - psa$inc_cost > 0),
                                      numeric(1)))
}

#' Pembrolizumab price discount needed to reach a threshold ICER
#'
#' The discount applies to every pembrolizumab administration in the model,
#' including subsequent anti-PD-(L)1 therapy in the chemotherapy arm, which
#' is costed as pembrolizumab.
#'
#' @param p,settings,ctx As in [run_once()].
#' @param threshold ICER threshold (GBP per QALY).
#' @return The discount (0–1), 0 if the list-price ICER is already below the
#'   threshold, or `NA` if no discount reaches it.
#' @export
threshold_discount <- function(p, settings, ctx, threshold) {
  f <- function(d) { q <- p; q$price_discount_pembro <- d; run_once(q, settings, ctx)[["icer"]] - threshold }
  f0 <- f(0)
  if (f0 <= 0) return(0)
  if (f(1) > 0) return(NA_real_)
  stats::uniroot(f, c(0, 1), tol = 1e-6)$root
}
