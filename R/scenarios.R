#' Scenario analyses mirroring TA531
#'
#' Each scenario changes one element of the base case: the OS or PFS
#' distribution, the TA531 piecewise approach (KM to 23, 33 or 43 weeks,
#' then an exponential tail), treatment-effect waning at 3 or 5 years,
#' treatment to progression without the stopping rule, the utility
#' approach, and the discount rate.
#'
#' @param base Base-case settings ([model_settings()]).
#' @return A named list. Each element is a list with `group`, `settings` and
#'   `modify` (a function applied to the parameter list).
#' @export
scenario_list <- function(base) {
  sc <- list()
  add <- function(name, group, settings = base, modify = identity) {
    sc[[name]] <<- list(group = group, settings = settings, modify = modify)
  }
  with_curves <- function(...) {
    s <- base
    args <- utils::modifyList(base$curves, list(...))
    s$curves <- do.call(curve_spec, args)
    s
  }
  for (w in c(23, 33, 43)) {
    add(sprintf("OS: KM to %d weeks + exponential (TA531)", w), "KM cut-off",
        with_curves(os_model = "Piecewise", piecewise_weeks = w))
  }
  for (m in c("Log-logistic", "Generalised gamma", "RP spline (1 knot)", "Gompertz", "Weibull",
              "Exponential")) {
    add(paste("OS distribution:", m), "Distribution", with_curves(os_model = m))
  }
  for (m in c("Log-logistic", "RP spline (1 knot)")) {
    add(paste("PFS distribution:", m), "Distribution", with_curves(pfs_model = m))
  }
  add("Waning: HR = 1 from 3 years", "Waning", with_curves(waning_years = 3, waning_hr = 1))
  add("Waning: HR = 1 from 5 years", "Waning", with_curves(waning_years = 5, waning_hr = 1))
  s <- base; s$tot <- "pfs_uncapped"
  add("No 2-year stopping rule (treat to progression)", "Stopping rule", s)
  s <- base; s$utility <- "progression"
  add("Progression-based utilities", "Utilities", s)
  s <- base; s$cap_ttd360 <- FALSE
  add("Time-to-death utilities, 360+ days not capped at norm", "Utilities", s)
  add("Discount rate 1.5% (costs and QALYs)", "Discounting", base,
      function(q) { q$discount_costs <- 0.015; q$discount_qalys <- 0.015; q })
  sc
}

#' Run all scenarios
#'
#' @param p Base-case parameters.
#' @param scenarios Output of [scenario_list()].
#' @param ctx Context as in [run_once()].
#' @return A data frame with one row per scenario.
#' @export
run_scenarios <- function(p, scenarios, ctx) {
  do.call(rbind, lapply(names(scenarios), function(nm) {
    s <- scenarios[[nm]]
    r <- run_once(s$modify(p), s$settings, ctx)
    data.frame(scenario = nm, group = s$group, inc_cost = r[["inc_cost"]],
               inc_qaly = r[["inc_qaly"]], icer = r[["icer"]], nmb_20000 = r[["nmb_20000"]],
               nmb_30000 = r[["nmb_30000"]], ly_undisc_pembro = r[["ly_undisc_pembro"]],
               ly_undisc_chemo = r[["ly_undisc_chemo"]])
  }))
}
