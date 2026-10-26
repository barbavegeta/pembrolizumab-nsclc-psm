# Pipeline: digitised KM coordinates -> reconstructed IPD -> validation ->
# survival model fitting and selection evidence.
# Run with targets::tar_make().
library(targets)
library(tarchetypes)

tar_source("R")  # functions go into the global environment, so targets tracks changes to them
tar_option_set(packages = c("survival", "ggplot2", "flexsurv"))

dig <- function(f) file.path("data-raw", "digitised", f)
psa_n_override <- suppressWarnings(as.integer(Sys.getenv("PSM_PSA_N", NA)))

curve_targets <- lapply(c("os_pembro", "os_chemo", "pfs_pembro", "pfs_chemo"), function(id) {
  list(
    tar_target_raw(paste0("file_", id), dig(paste0(id, ".csv")), format = "file"),
    tar_target_raw(paste0("clean_", id), substitute(
      clean_digitised_curve(read_digitised_curve(F), x_max = axis_max(axis_limits, EP), label = ID),
      list(F = as.symbol(paste0("file_", id)), EP = sub("_.*", "", id), ID = id)
    ))
  )
})

list(
  curve_targets,
  tar_target(file_nar_os, dig("nar_os.csv"), format = "file"),
  tar_target(file_nar_pfs, dig("nar_pfs.csv"), format = "file"),
  tar_target(file_axis_limits, dig("axis_limits.csv"), format = "file"),
  tar_target(file_recon_inputs, "data-raw/reconstruction_inputs.csv", format = "file"),
  tar_target(nar_os, read_nar(file_nar_os)),
  tar_target(nar_pfs, read_nar(file_nar_pfs)),
  tar_target(axis_limits, read.csv(file_axis_limits)),
  tar_target(recon_inputs, read.csv(file_recon_inputs)),

  tar_target(ipd_os, rbind(
    reconstruct_from_inputs(clean_os_pembro$curve, nar_os, "pembro", recon_inputs, "os"),
    reconstruct_from_inputs(clean_os_chemo$curve, nar_os, "chemo", recon_inputs, "os"))),
  tar_target(ipd_pfs, rbind(
    reconstruct_from_inputs(clean_pfs_pembro$curve, nar_pfs, "pembro", recon_inputs, "pfs"),
    reconstruct_from_inputs(clean_pfs_chemo$curve, nar_pfs, "chemo", recon_inputs, "pfs"))),

  tar_target(validation, list(
    os  = validate_reconstruction(ipd_os, nar_os, published_targets()$os,
                                  acceptance_thresholds(), "os"),
    pfs = validate_reconstruction(ipd_pfs, nar_pfs, published_targets()$pfs,
                                  acceptance_thresholds(), "pfs"))),

  tar_target(plot_os_file, save_plot(plot_reconstruction(
    ipd_os, stack_curves(clean_os_pembro, clean_os_chemo), "Overall survival"),
    "output/validation/overlay_os.png"), format = "file"),
  tar_target(plot_pfs_file, save_plot(plot_reconstruction(
    ipd_pfs, stack_curves(clean_pfs_pembro, clean_pfs_chemo), "Progression-free survival"),
    "output/validation/overlay_pfs.png"), format = "file"),

  tar_target(reconstruction_report, write_reconstruction_report(
    validation,
    cleaning_log = rbind(clean_os_pembro$log, clean_os_chemo$log,
                         clean_pfs_pembro$log, clean_pfs_chemo$log),
    methods = describe_methods(list(os = ipd_os, pfs = ipd_pfs)),
    plot_files = c(os = basename(plot_os_file), pfs = basename(plot_pfs_file)),
    path = "output/validation/reconstruction.md"), format = "file"),

  tar_target(reconstruction_accepted, {
    reconstruction_report
    assert_reconstruction(validation)
  }),

  tar_target(ipd_os_file, write_ipd(ipd_os, "output/ipd/ipd_os.csv", reconstruction_accepted),
             format = "file"),
  tar_target(ipd_pfs_file, write_ipd(ipd_pfs, "output/ipd/ipd_pfs.csv", reconstruction_accepted),
             format = "file"),

  # Step 2: survival extrapolation ------------------------------------------
  tar_target(file_inputs, "data-raw/inputs.csv", format = "file"),
  tar_target(inputs, read_inputs(file_inputs)),
  tar_target(file_lifetable, file.path("data-raw", input_value(inputs, "lifetable", FALSE)),
             format = "file"),
  tar_target(lifetable, read.csv(file_lifetable)),
  tar_target(population, list(
    start_age = input_value(inputs, "start_age"),
    prop_male = input_value(inputs, "prop_male"),
    lifetable_edition = inputs$edition_year[inputs$parameter == "lifetable"])),
  tar_target(grid, extrapolation_grid(20)),
  tar_target(s_gpm, gpm_survival(lifetable, population$start_age, population$prop_male, grid)),

  tar_target(surv_ipd, {
    stopifnot(isTRUE(reconstruction_accepted))
    list(os = ipd_os, pfs = ipd_pfs)
  }),
  tar_target(curve_id, c("os_pembro", "os_chemo", "pfs_pembro", "pfs_chemo")),
  tar_target(curve_ipd, {
    ep <- sub("_.*", "", curve_id); a <- sub(".*_", "", curve_id)
    d <- surv_ipd[[ep]]
    list(id = curve_id, endpoint = ep, arm = a, ipd = d[d$arm == a, c("time", "status")])
  }, pattern = map(curve_id), iteration = "list"),
  tar_target(surv_fits, list(id = curve_ipd$id, fits = fit_survival_models(curve_ipd$ipd)),
             pattern = map(curve_ipd), iteration = "list"),
  tar_target(selection, model_selection_table(surv_fits$fits, curve_ipd$ipd, curve_ipd$endpoint,
                                              curve_ipd$arm, s_gpm, grid),
             pattern = map(surv_fits, curve_ipd), iteration = "list"),
  tar_target(piecewise, piecewise_table(curve_ipd$ipd, curve_ipd$endpoint, curve_ipd$arm,
                                        c(23, 33, 43), s_gpm = s_gpm, grid = grid),
             pattern = map(curve_ipd), iteration = "list"),
  tar_target(ph, list(os = ph_test(ipd_os), pfs = ph_test(ipd_pfs))),

  tar_target(fig_extrap, save_plot(plot_extrapolations(
      surv_fits$fits, curve_ipd$ipd, grid, curve_title(curve_ipd$id),
      if (curve_ipd$endpoint == "os") s_gpm else NULL),
    sprintf("output/survival/extrap_%s.png", curve_ipd$id), height = 7),
    pattern = map(surv_fits, curve_ipd), format = "file"),
  tar_target(fig_haz, save_plot(plot_hazards(surv_fits$fits, curve_ipd$ipd,
                                             curve_title(curve_ipd$id)),
    sprintf("output/survival/haz_%s.png", curve_ipd$id), height = 7),
    pattern = map(surv_fits, curve_ipd), format = "file"),
  tar_target(fig_lch, c(
    lch_os = save_plot(plot_log_cumhaz(ipd_os, "Overall survival"), "output/survival/lch_os.png"),
    lch_pfs = save_plot(plot_log_cumhaz(ipd_pfs, "Progression-free survival"),
                        "output/survival/lch_pfs.png")), format = "file"),

  tar_target(selection_csv, {
    p <- "output/survival/model_selection.csv"
    write.csv(dplyr::bind_rows(selection), p, row.names = FALSE); p
  }, format = "file"),
  tar_target(survival_report, write_survival_report(
    dplyr::bind_rows(selection), dplyr::bind_rows(piecewise), ph, population,
    figures = survival_figure_index(fig_extrap, fig_haz, fig_lch),
    path = "output/survival/model_selection.md"), format = "file"),

  # Step 3: partitioned survival model --------------------------------------
  tar_target(model_files, c(file.path("data-raw", c(
      "inputs.csv", "regimens.csv", "adverse_events.csv", "subsequent_therapy.csv", "nhscii.csv",
      "eq5d_norms_dsu2022.csv", "lifetables_england_2022_2024.csv", "emit_2025_26.csv",
      "ncc_2024_25.csv")), "data-raw/inputs.csv"), format = "file"),
  tar_target(model_data, { model_files; load_model_data("data-raw") }),
  tar_target(params, model_params(model_data)),
  tar_target(mtime, model_time(params$horizon_years, params$cycle_weeks)),
  tar_target(t_ext, c(mtime$t_months, max(mtime$t_months) + weeks_to_months(1:53))),
  tar_target(s_gpm_ext, gpm_survival(model_data$lifetable, params$start_age, params$prop_male, t_ext)),
  tar_target(fits_by_curve, stats::setNames(lapply(surv_fits, `[[`, "fits"),
                                            vapply(surv_fits, `[[`, "", "id"))),
  tar_target(ipd_list, list(os = ipd_os, pfs = ipd_pfs)),
  tar_target(base_settings, model_settings(
    curves = curve_spec(os_model = input_value(model_data$inputs, "os_model", FALSE),
                        pfs_model = input_value(model_data$inputs, "pfs_model", FALSE)))),
  tar_target(base_curves, build_curves(base_settings$curves, fits_by_curve, ipd_list, t_ext,
                                       s_gpm_ext)),
  tar_target(base_case, {
    run <- run_model(params, base_curves, base_settings, mtime, t_ext, model_data$norms)
    list(run = run, summary = summarise_results(run))
  }),
  tar_target(trace_files, write_traces(base_case$run, "output/model"), format = "file"),
  tar_target(trace_excel, export_trace_excel(base_case$run, "output/model/base_case_trace.xlsx"),
             format = "file"),

  # Step 4: analyses ---------------------------------------------------------
  tar_target(ctx, list(fits = fits_by_curve, ipd = ipd_list, t_ext = t_ext, s_gpm = s_gpm_ext,
                       time = mtime, norms = model_data$norms)),
  # PSM_PSA_N (fast CI profile) is read when the pipeline is defined, so the
  # target changes only when the requested number of iterations changes.
  tar_target(psa_n, {
    n <- !!psa_n_override
    if (is.na(n)) as.integer(params$psa_iterations) else n
  }),
  tar_target(owsa, run_owsa(params, base_settings, ctx, model_data)),
  tar_target(scenarios, run_scenarios(params, scenario_list(base_settings), ctx)),
  tar_target(thresholds, data.frame(
    threshold = c(20000, 30000),
    discount = vapply(c(20000, 30000), function(k) threshold_discount(params, base_settings, ctx, k),
                      numeric(1)))),
  tar_target(psa_draws, draw_psa(params, model_data, fits_by_curve, base_settings$curves,
                                 psa_n, params$psa_seed)),
  tar_target(psa, run_psa(psa_draws, base_settings, ctx)),
  tar_target(psa_ceac, ceac(psa)),
  tar_target(validation_ta531, ta531_validation(base_case$summary, scenarios)),
  tar_target(diagnostics, model_diagnostics(base_curves, params, mtime, model_data$norms,
                                            fits_by_curve, owsa)),
  tar_target(file_facts, "data-raw/published_facts.csv", format = "file"),
  tar_target(facts, {
    f <- read.csv(file_facts)
    stats::setNames(as.list(f$value), f$key)
  }),

  tar_target(results_files, {
    dir <- "output/results"
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    out <- c(
      save_csv(as.data.frame(base_case$summary$totals), file.path(dir, "base_case_totals.csv"),
               row.names = TRUE),
      save_csv(as.data.frame(t(base_case$summary$incremental)), file.path(dir, "base_case_incremental.csv")),
      save_csv(owsa, file.path(dir, "owsa.csv")),
      save_csv(scenarios, file.path(dir, "scenarios.csv")),
      save_csv(thresholds, file.path(dir, "threshold_discount.csv")),
      save_csv(psa, file.path(dir, "psa.csv")),
      save_csv(psa_ceac, file.path(dir, "ceac.csv")),
      save_csv(validation_ta531, file.path(dir, "validation_ta531.csv")),
      save_plot(plot_ce_plane(psa, base_case$summary$incremental), file.path(dir, "ce_plane.png")),
      save_plot(plot_ceac(psa_ceac), file.path(dir, "ceac.png")),
      save_plot(plot_tornado(owsa, base_case$summary$incremental[["icer"]]),
                file.path(dir, "tornado.png"), height = 5.5),
      save_plot(plot_base_curves(base_curves, t_ext, ipd_list), file.path(dir, "base_curves.png"))
    )
    out
  }, format = "file"),

  # Report (skipped when the Quarto CLI is unavailable, e.g. the fast test profile)
  if (nzchar(Sys.which("quarto")) && !identical(Sys.getenv("PSM_SKIP_REPORT"), "true")) {
    tar_quarto(report, path = "report/report.qmd", quiet = FALSE)
  }
)
