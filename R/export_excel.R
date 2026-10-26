#' Export the base-case trace to Excel with transparent formulas
#'
#' Writes one sheet per arm with the per-cycle trace as values (state
#' membership, half-cycle-corrected occupancy, utilities, cost streams and
#' discount factors). Total and discounted costs, QALYs and life years per
#' cycle are Excel formulas. A summary sheet sums them with formulas and
#' computes the incremental results and the ICER.
#'
#' @param run Output of [run_model()].
#' @param path Output `.xlsx` path.
#' @return `path`.
#' @export
export_trace_excel <- function(run, path) {
  sheets <- list()
  refs <- list()
  for (a in c("pembro", "chemo")) {
    sheet <- if (a == "pembro") "Pembrolizumab" else "Chemotherapy"
    tr <- run[[a]]$trace
    cost_cols <- setdiff(grep("^cost_", names(tr), value = TRUE), "cost_total")
    d <- tr[, c("cycle", "t_start_years", "pf_start", "pd_start", "dead_start", "occ_pf", "occ_pd",
                "u_pf", "u_pd", "ly", "qaly", cost_cols, "df_costs", "df_qalys")]
    n <- nrow(d)
    rows <- seq_len(n) + 1
    L <- function(nm) .xl_col(match(nm, names(d)))
    total_col <- .xl_col(ncol(d) + 1)
    d$cost_total <- writexl::xl_formula(sprintf("=SUM(%s%d:%s%d)", L(cost_cols[1]), rows,
                                                L(cost_cols[length(cost_cols)]), rows))
    d$cost_disc <- writexl::xl_formula(sprintf("=%s%d*%s%d", total_col, rows, L("df_costs"), rows))
    d$qaly_disc <- writexl::xl_formula(sprintf("=%s%d*%s%d", L("qaly"), rows, L("df_qalys"), rows))
    d$ly_disc <- writexl::xl_formula(sprintf("=%s%d*%s%d", L("ly"), rows, L("df_qalys"), rows))
    rng <- function(nm) sprintf("SUM('%s'!%s2:%s%d)", sheet, .xl_col(match(nm, names(d))),
                                .xl_col(match(nm, names(d))), n + 1)
    refs[[a]] <- c(cost = rng("cost_disc"), qaly = rng("qaly_disc"), ly = rng("ly_disc"),
                   cost_undisc = rng("cost_total"), qaly_undisc = rng("qaly"), ly_undisc = rng("ly"))
    sheets[[sheet]] <- d
  }
  out <- c("cost", "qaly", "ly", "cost_undisc", "qaly_undisc", "ly_undisc")
  lab <- c("Discounted costs (GBP)", "Discounted QALYs", "Discounted life years",
           "Undiscounted costs (GBP)", "Undiscounted QALYs", "Undiscounted life years")
  sheets[["Summary"]] <- data.frame(
    Outcome = c(lab, "ICER (GBP per QALY)"),
    Pembrolizumab = writexl::xl_formula(c(paste0("=", refs$pembro[out]), "=\"\"")),
    Chemotherapy = writexl::xl_formula(c(paste0("=", refs$chemo[out]), "=\"\"")),
    Incremental = writexl::xl_formula(c(sprintf("=B%d-C%d", 2:7, 2:7), "=D2/D3"))
  )
  sheets[["Notes"]] <- data.frame(Note = c(
    "Values in the arm sheets come from the R pipeline (targets::tar_make()).",
    "Occupancy (occ_pf, occ_pd) uses the trapezoidal half-cycle correction.",
    "Discount factors (df_costs, df_qalys) are evaluated at cycle midpoints.",
    "cost_total, cost_disc, qaly_disc and ly_disc are formulas; the Summary sheet sums them.",
    "Costs are 2024/25 GBP at list prices."))
  sheets <- sheets[c("Summary", "Pembrolizumab", "Chemotherapy", "Notes")]
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(sheets, path)
  path
}

#' Excel column letters
#' @noRd
.xl_col <- function(i) {
  vapply(i, function(k) {
    s <- ""
    while (k > 0) { r <- (k - 1) %% 26; s <- paste0(LETTERS[r + 1], s); k <- (k - 1) %/% 26 }
    s
  }, character(1))
}
