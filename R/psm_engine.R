#' Model time grid
#'
#' @param horizon_years Time horizon in years.
#' @param cycle_weeks Cycle length in weeks.
#' @return A list with `n_cycles`, `cycle_years`, `t_weeks` (cycle
#'   boundaries from 0), `t_months` and `t_years`, each of length
#'   `n_cycles + 1`.
#' @export
model_time <- function(horizon_years, cycle_weeks = 1) {
  cycle_years <- cycle_weeks * 7 / 365.25
  n <- round(horizon_years / cycle_years)
  t_weeks <- seq(0, n) * cycle_weeks
  list(n_cycles = n, cycle_years = cycle_years, t_weeks = t_weeks,
       t_months = weeks_to_months(t_weeks), t_years = t_weeks * 7 / 365.25)
}

#' Partitioned survival state membership
#'
#' @param s_pfs,s_os PFS and OS at the cycle boundaries (same length,
#'   starting at 1). PFS must not exceed OS (apply
#'   [apply_survival_constraints()] first).
#' @return A data frame with `pf`, `pd` and `dead` at each cycle boundary.
#'   The three columns sum to 1 in every row.
#' @export
psm_states <- function(s_pfs, s_os) {
  stopifnot(length(s_pfs) == length(s_os), all(s_pfs <= s_os + 1e-12))
  data.frame(pf = s_pfs, pd = s_os - s_pfs, dead = 1 - s_os)
}

#' Half-cycle correction by the trapezoidal (life-table) method
#'
#' @param x Values at cycle boundaries (length `n + 1`).
#' @return Mean value within each cycle, `(x[j] + x[j + 1]) / 2`
#'   (length `n`).
#' @export
half_cycle_trapezoid <- function(x) (utils::head(x, -1) + utils::tail(x, -1)) / 2

#' Discount factors at cycle midpoints
#'
#' @param n_cycles Number of cycles.
#' @param cycle_years Cycle length in years.
#' @param rate Annual discount rate.
#' @return `(1 + rate)^-t` with `t` the midpoint of each cycle in years.
#' @export
discount_midpoint <- function(n_cycles, cycle_years, rate) {
  (1 + rate)^-((seq_len(n_cycles) - 0.5) * cycle_years)
}

#' Run the partitioned survival engine for one arm
#'
#' Generic engine: it knows nothing about the decision problem. It takes
#' constrained PFS and OS at the cycle boundaries, the per-cycle utility of
#' each living state, per-cycle state costs and any other per-cycle cost
#' streams already expressed per member of the starting cohort. It returns
#' the trace and the discounted and undiscounted totals.
#'
#' State occupancy in each cycle uses the trapezoidal half-cycle correction
#' ([half_cycle_trapezoid()]). Costs and QALYs are discounted at cycle
#' midpoints ([discount_midpoint()]).
#'
#' @param s_pfs,s_os Survival at cycle boundaries (length `n_cycles + 1`).
#' @param cycle_years Cycle length in years.
#' @param u_pf,u_pd Utility in each state per cycle (length `n_cycles` or 1).
#' @param c_pf,c_pd Cost per cycle of being in each state (length
#'   `n_cycles` or 1), applied to half-cycle-corrected occupancy.
#' @param cost_streams Named list of cohort-level costs per cycle (each of
#'   length `n_cycles`), e.g. drug acquisition or one-off event costs.
#' @param rate_costs,rate_qalys Annual discount rates.
#' @return A list with `trace` (one row per cycle) and `totals` (a named
#'   numeric vector of undiscounted and discounted life years, QALYs and
#'   costs by category).
#' @export
run_psm <- function(s_pfs, s_os, cycle_years, u_pf, u_pd, c_pf = 0, c_pd = 0,
                    cost_streams = list(), rate_costs = 0.035, rate_qalys = 0.035) {
  states <- psm_states(s_pfs, s_os)
  n <- nrow(states) - 1
  occ_pf <- half_cycle_trapezoid(states$pf)
  occ_pd <- half_cycle_trapezoid(states$pd)
  rec <- function(x) if (length(x) == 1) rep(x, n) else { stopifnot(length(x) == n); x }
  u_pf <- rec(u_pf); u_pd <- rec(u_pd); c_pf <- rec(c_pf); c_pd <- rec(c_pd)

  ly <- (occ_pf + occ_pd) * cycle_years
  qaly <- (occ_pf * u_pf + occ_pd * u_pd) * cycle_years
  streams <- c(list(state_pf = occ_pf * c_pf, state_pd = occ_pd * c_pd),
               lapply(cost_streams, rec))
  df_c <- discount_midpoint(n, cycle_years, rate_costs)
  df_q <- discount_midpoint(n, cycle_years, rate_qalys)

  trace <- data.frame(cycle = seq_len(n), t_start_years = (seq_len(n) - 1) * cycle_years,
                      pf_start = utils::head(states$pf, -1), pd_start = utils::head(states$pd, -1),
                      dead_start = utils::head(states$dead, -1), occ_pf = occ_pf, occ_pd = occ_pd,
                      u_pf = u_pf, u_pd = u_pd, ly = ly, qaly = qaly,
                      df_costs = df_c, df_qalys = df_q)
  for (nm in names(streams)) trace[[paste0("cost_", nm)]] <- streams[[nm]]
  cost_total <- Reduce(`+`, streams)
  trace$cost_total <- cost_total

  totals <- c(ly_undisc = sum(ly), ly_disc = sum(ly * df_q),
              qaly_undisc = sum(qaly), qaly_disc = sum(qaly * df_q),
              cost_undisc = sum(cost_total), cost_disc = sum(cost_total * df_c))
  by_cat <- vapply(streams, function(s) sum(s * df_c), numeric(1))
  names(by_cat) <- paste0("cost_disc_", names(by_cat))
  list(trace = trace, totals = c(totals, by_cat))
}
