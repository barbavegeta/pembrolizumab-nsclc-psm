#' Load all model input files
#'
#' @param dir Directory holding the input files (`data-raw`).
#' @return A list with `inputs`, `regimens`, `adverse_events`,
#'   `subsequent`, `nhscii`, `norms` and `lifetable`.
#' @export
load_model_data <- function(dir = "data-raw") {
  inputs <- read_inputs(file.path(dir, "inputs.csv"))
  rd <- function(f) utils::read.csv(file.path(dir, f), stringsAsFactors = FALSE)
  list(inputs = inputs, regimens = rd("regimens.csv"), adverse_events = rd("adverse_events.csv"),
       subsequent = rd("subsequent_therapy.csv"), nhscii = rd(input_value(inputs, "nhscii", FALSE)),
       norms = rd(input_value(inputs, "utility_norms", FALSE)),
       lifetable = rd(input_value(inputs, "lifetable", FALSE)))
}

#' Deterministic parameter set
#'
#' Collects every numeric input into one named list. Tables (regimen mix,
#' AE rates, subsequent-therapy use) are carried as vectors so that the
#' PSA can resample them.
#'
#' @param data Output of [load_model_data()].
#' @return A named list of parameters.
#' @export
model_params <- function(data) {
  inp <- data$inputs
  num <- suppressWarnings(as.numeric(inp$value))
  p <- as.list(stats::setNames(num, inp$parameter))
  p <- p[!is.na(num)]
  p$inflation_201415 <- prod(1 + data$nhscii$pct_increase_overall / 100)
  p$regimen_share <- stats::setNames(data$regimens$share_keynote024, data$regimens$regimen)
  p$ae_rate_pembro <- stats::setNames(data$adverse_events$rate_pembro, data$adverse_events$ae)
  p$ae_rate_chemo <- stats::setNames(data$adverse_events$rate_chemo, data$adverse_events$ae)
  p$ae_cost_201415 <- stats::setNames(data$adverse_events$unit_cost_201415, data$adverse_events$ae)
  s <- data$subsequent
  p$subsq_prop_pembro <- stats::setNames(s$n_pembro / s$N_pembro, s$therapy)
  p$subsq_prop_chemo <- stats::setNames(s$n_chemo / s$N_chemo, s$therapy)
  p$price_discount_pembro <- 0
  p
}

#' Model settings (scenario switches)
#'
#' @param utility `"ttd"` (time to death, base case) or `"progression"`.
#' @param cap_ttd360 Cap the 360+ day utility at the population norm.
#' @param tot `"pfs_capped"` (PFS with the stopping rule, base case) or
#'   `"pfs_uncapped"` (treatment to progression).
#' @param curves Output of [curve_spec()].
#' @return A list.
#' @export
model_settings <- function(utility = "ttd", cap_ttd360 = TRUE, tot = "pfs_capped",
                           curves = curve_spec()) {
  list(utility = match.arg(utility, c("ttd", "progression")), cap_ttd360 = cap_ttd360,
       tot = match.arg(tot, c("pfs_capped", "pfs_uncapped")), curves = curves)
}

#' Unit costs derived from the parameters
#'
#' @param p Parameter list from [model_params()].
#' @return A named list of unit costs (2024/25 GBP): `pembro_dose`
#'   (acquisition plus administration, after any price discount),
#'   `pembro_drug_dose`, `chemo_cycle`, `chemo_cycle_drug`, `pem_maint_dose`,
#'   `docetaxel_cycle`, `subsequent_course` (named by therapy), `ae_pembro`,
#'   `ae_chemo`, `pf_week`, `pd_week` and `terminal`.
#' @export
unit_costs <- function(p) {
  bsa <- p$bsa
  drug_mg <- function(d) p[[paste0("dose_mg_m2_", d)]] * bsa * p[[paste0("price_per_mg_", d)]]
  admin <- function(items) sum(vapply(strsplit(items, "\\+")[[1]], function(i) p[[i]], numeric(1)))

  pembro_drug <- p$pembro_dose_mg / 100 * p$price_pembro_100mg * (1 - p$price_discount_pembro) *
    p$dose_intensity_pembro
  pembro_dose <- pembro_drug + p$admin_sb12z_op

  reg <- data.frame(
    regimen = c("gem_carb", "gem_cis", "pac_carb", "pem_cis", "pem_carb"),
    platinum = c("carboplatin", "cisplatin", "carboplatin", "cisplatin", "carboplatin"),
    partner = c("gemcitabine", "gemcitabine", "paclitaxel", "pemetrexed", "pemetrexed"),
    partner_doses = c(2, 2, 1, 1, 1),
    admin = c("admin_sb14z_op+admin_sb15z_op", "admin_sb14z_dc+admin_sb15z_op",
              "admin_sb14z_op", "admin_sb14z_dc", "admin_sb14z_op")
  )
  share <- p$regimen_share[reg$regimen] / sum(p$regimen_share[reg$regimen])
  drug_cycle <- vapply(seq_len(nrow(reg)), function(i) {
    (drug_mg(reg$platinum[i]) + reg$partner_doses[i] * drug_mg(reg$partner[i])) * p$dose_intensity_chemo
  }, numeric(1))
  admin_cycle <- vapply(reg$admin, admin, numeric(1))
  chemo_cycle_drug <- sum(share * drug_cycle)
  chemo_cycle <- sum(share * (drug_cycle + admin_cycle))

  pem_maint <- drug_mg("pemetrexed") * p$dose_intensity_chemo + p$admin_sb12z_op
  docetaxel_cycle <- drug_mg("docetaxel") + p$admin_sb12z_op

  course <- stats::setNames(
    unname(c(p$subsequent_io_weeks / p$pembro_interval_weeks * pembro_dose,
             p$subsequent_doublet_cycles * chemo_cycle,
             p$subsequent_docetaxel_cycles * docetaxel_cycle)),
    c("Immunotherapy (anti-PD-(L)1)", "Platinum doublet chemotherapy", "Single-agent chemotherapy"))
  infl <- p$inflation_201415
  list(pembro_dose = pembro_dose, pembro_drug_dose = pembro_drug, chemo_cycle = chemo_cycle,
       chemo_cycle_drug = chemo_cycle_drug, pem_maint_dose = pem_maint,
       docetaxel_cycle = docetaxel_cycle, subsequent_course = course,
       ae_pembro = sum(p$ae_rate_pembro * p$ae_cost_201415) * infl,
       ae_chemo = sum(p$ae_rate_chemo * p$ae_cost_201415) * infl,
       pf_week = p$cost_pf_week_201415 * infl, pd_week = p$cost_pd_week_201415 * infl,
       terminal = p$cost_terminal_201415 * infl)
}

#' Proportion of the cohort dosed in each cycle
#'
#' A dose is given at the start of each cycle whose start time is in the
#' dosing schedule, to the proportion `s` on treatment at that time.
#'
#' @param s Proportion on treatment at cycle boundaries (length `n + 1`).
#' @param t_weeks Cycle boundaries in weeks.
#' @param start_week First dose (weeks).
#' @param interval_weeks Dosing interval (weeks).
#' @param max_doses Maximum number of doses (`Inf` for none).
#' @return A vector of length `n`: the proportion dosed in each cycle.
#' @export
dosing_stream <- function(s, t_weeks, start_week, interval_weeks, max_doses = Inf) {
  starts <- utils::head(t_weeks, -1)
  k <- (starts - start_week) / interval_weeks
  dose <- starts >= start_week - 1e-9 & abs(k - round(k)) < 1e-9 & round(k) < max_doses
  ifelse(dose, utils::head(s, -1), 0)
}

#' Run the model for both arms
#'
#' @param p Parameters from [model_params()].
#' @param curves Output of [build_curves()] on the extended grid `t_ext`.
#' @param settings Output of [model_settings()].
#' @param time Output of [model_time()].
#' @param t_ext Extended grid (months): `time$t_months` followed by at least
#'   360 days more, for time-to-death utilities.
#' @param norms EQ-5D norms table.
#' @return A list with `pembro` and `chemo` ([run_psm()] outputs) and
#'   `unit_costs`.
#' @export
run_model <- function(p, curves, settings, time, t_ext, norms) {
  n <- time$n_cycles
  idx <- seq_len(n + 1)
  cycle_weeks <- time$cycle_years * 365.25 / 7
  uc <- unit_costs(p)
  t_mid <- (utils::head(time$t_months, -1) + utils::tail(time$t_months, -1)) / 2
  age_mid <- p$start_age + t_mid / 12
  norm0 <- norm_utility(norms, p$start_age, p$prop_male)
  age_adj <- norm_utility(norms, age_mid, p$prop_male) / norm0

  arm_run <- function(a) {
    os <- curves[[a]]$os[idx]
    pfs <- curves[[a]]$pfs[idx]

    if (settings$utility == "ttd") {
      u <- p[c("u_ttd_0", "u_ttd_30", "u_ttd_180", "u_ttd_360")]
      if (settings$cap_ttd360) u$u_ttd_360 <- min(u$u_ttd_360, norm0)
      u_alive <- ttd_utility(curves[[a]]$os, t_ext, t_mid, u) * age_adj
      u_pf <- u_pd <- u_alive
    } else {
      u_pf <- p$u_pf * age_adj
      u_pd <- p$u_pd * age_adj
    }

    deaths <- -diff(os)
    pfs_events <- -diff(pfs)
    streams <- list()
    if (a == "pembro") {
      max_d <- if (settings$tot == "pfs_capped") p$pembro_max_doses else Inf
      dosed <- dosing_stream(pfs, time$t_weeks, 0, p$pembro_interval_weeks, max_d)
      streams$drug_1l <- dosed * uc$pembro_drug_dose
      streams$admin_1l <- dosed * p$admin_sb12z_op
      streams$maintenance <- rep(0, n)
      streams$ae <- c(uc$ae_pembro, rep(0, n - 1))
      subsq <- sum(p$subsq_prop_pembro * uc$subsequent_course[names(p$subsq_prop_pembro)])
    } else {
      dosed <- dosing_stream(pfs, time$t_weeks, 0, p$chemo_interval_weeks, p$chemo_max_cycles)
      streams$drug_1l <- dosed * uc$chemo_cycle_drug
      streams$admin_1l <- dosed * (uc$chemo_cycle - uc$chemo_cycle_drug)
      s_start <- stats::approx(time$t_weeks, pfs, xout = p$pem_maint_start_week)$y
      on_maint <- p$pem_maint_prop * pmin(pfs / s_start, 1)
      streams$maintenance <- dosing_stream(on_maint, time$t_weeks, p$pem_maint_start_week,
                                           p$chemo_interval_weeks) * uc$pem_maint_dose
      streams$ae <- c(uc$ae_chemo, rep(0, n - 1))
      subsq <- sum(p$subsq_prop_chemo * uc$subsequent_course[names(p$subsq_prop_chemo)])
    }
    streams$subsequent <- pfs_events * subsq
    streams$terminal <- deaths * uc$terminal

    run_psm(pfs, os, time$cycle_years, u_pf = u_pf, u_pd = u_pd,
            c_pf = uc$pf_week * cycle_weeks, c_pd = uc$pd_week * cycle_weeks,
            cost_streams = streams, rate_costs = p$discount_costs, rate_qalys = p$discount_qalys)
  }
  list(pembro = arm_run("pembro"), chemo = arm_run("chemo"), unit_costs = uc)
}

#' Cost-effectiveness summary of a model run
#'
#' @param res Output of [run_model()].
#' @param thresholds Willingness-to-pay thresholds (GBP per QALY).
#' @return A list with `totals` (one row per arm) and `incremental`
#'   (costs, LYs, QALYs, ICER and net monetary benefit at each threshold).
#' @export
summarise_results <- function(res, thresholds = c(20000, 30000)) {
  stopifnot(identical(names(res$pembro$totals), names(res$chemo$totals)))
  tot <- rbind(pembro = res$pembro$totals, chemo = res$chemo$totals)
  d <- tot["pembro", ] - tot["chemo", ]
  inc <- c(inc_cost = d[["cost_disc"]], inc_ly = d[["ly_disc"]], inc_ly_undisc = d[["ly_undisc"]],
           inc_qaly = d[["qaly_disc"]], icer = d[["cost_disc"]] / d[["qaly_disc"]])
  for (k in thresholds) inc[[paste0("nmb_", k)]] <- k * d[["qaly_disc"]] - d[["cost_disc"]]
  list(totals = tot, incremental = inc)
}
