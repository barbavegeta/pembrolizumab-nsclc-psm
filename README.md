# pembrolizumab-nsclc-psm

A reproducible R implementation of a three-state partitioned survival
cost-effectiveness model of first-line **pembrolizumab versus platinum-based
chemotherapy** for metastatic non-small-cell lung cancer (NSCLC) with PD-L1
tumour proportion score ≥ 50%.

The model uses the five-year KEYNOTE-024 analysis (Reck et al. 2021) and
follows the structure of NICE technology appraisal TA531. Every input is
traced to a public source, survival extrapolation follows NICE DSU TSD 14,
19 and 21, and `targets::tar_make()` rebuilds everything from the digitised
data.

> **This is not NICE's model, and its results are not NICE's results.**
> TA531's decision-making ICERs include confidential discounts. This model
> uses list prices and gives a price-threshold analysis instead. The report
> explains why its results differ from the appraisal.

## What the pipeline does

1. **Reconstructs individual patient data** from the Kaplan–Meier curves in
   Figure 2 of Reck et al. (2021), using the Guyot et al. (2012) algorithm
   (`IPDfromKM`). The reconstruction is checked against the published
   medians, hazard ratios, landmark survival and numbers at risk, with
   acceptance thresholds set in advance; the pipeline stops if any check
   fails (`output/validation/reconstruction.md`).
2. **Fits and compares survival models** for OS and PFS in each arm:
   exponential, Weibull, Gompertz, log-logistic, lognormal, generalised
   gamma, and Royston–Parmar splines with 1–3 knots. It reports AIC/BIC,
   smoothed hazards, log-cumulative-hazard plots, extrapolated 5-, 10- and
   20-year survival, and the TA531 piecewise approach
   (`output/survival/model_selection.md`). OS is never allowed to fall below
   general-population mortality, and PFS never exceeds OS.
3. **Runs the partitioned survival model**: weekly cycles, trapezoidal
   half-cycle correction, a 20-year horizon, an NHS and PSS perspective, and
   3.5% discounting. Pembrolizumab has a 2-year stopping rule. Utilities are
   based on time to death, with age adjustment. Costs cover drug acquisition,
   administration, monitoring, adverse events, terminal care, and subsequent
   therapy including immunotherapy after crossover.
4. **Analyses**: base case with net monetary benefit, one-way sensitivity
   analysis (tornado), probabilistic sensitivity analysis (1,000 iterations,
   fixed seed; CE plane and CEAC), TA531-style scenarios, a pembrolizumab
   price-threshold analysis, and validation against TA531's modelled survival.
5. **Renders the report** `report/report.qmd` to HTML. Every number in it is
   read from the pipeline outputs.

## Repository layout

| Path | Contents |
|---|---|
| `R/` | Functions, documented with roxygen. The engine (`psm_engine.R`, `constraints.R`, `curves.R`) is generic and separate from the TA531-specific inputs (`model.R`). |
| `_targets.R` | The pipeline. |
| `data-raw/` | Digitised coordinates, the extraction script, and every model input with its source (`inputs.csv`, `regimens.csv`, `adverse_events.csv`, `subsequent_therapy.csv`, plus the prepared price, unit-cost, life-table and utility-norm files). |
| `data-raw/prepare_*.R` | Scripts that download each public source, check its md5 against the edition used, and write the CSVs. |
| `docs/source_checks.md` | The exact sentence or table cell behind every fact used, with its location. |
| `DECISIONS.md` | Every methodological decision: who made it, and why. |
| `tests/testthat/` | Unit and integration tests. |
| `output/` | Pipeline outputs: validation, survival evidence, traces and results. |
| `report/report.qmd`, `report/report.html` | The technical report and its rendered, self-contained HTML. |
| `output/model/base_case_trace.xlsx` | Base-case trace for both arms, with the cost and QALY sums and the ICER as Excel formulas. |
| `app/app.R` | A small Shiny app for running scenarios with the same model functions (`shiny::runApp("app")` after `tar_make()`). |

## Running it

Requirements: R 4.3, the system library `libglpk40` (needed by `igraph`, a
dependency of `targets`), and the Quarto CLI to render the report.

```r
renv::restore()        # packages pinned in renv.lock (Posit Package Manager snapshot 2024-06-03)
targets::tar_make()    # full pipeline, including the 1,000-iteration PSA and the report
testthat::test_local() # tests
```

The fast profile used in continuous integration runs 50 PSA iterations and
skips the report:

```sh
PSM_PSA_N=50 PSM_SKIP_REPORT=true Rscript -e 'targets::tar_make()'
```

The GitHub Actions workflow (`.github/workflows/check.yml`) runs the fast
pipeline and the tests on every push.

To rebuild the inputs from their public sources:

```sh
Rscript data-raw/prepare_lifetables.R
Rscript data-raw/prepare_unit_costs.R
python data-raw/extract_km_vectors.py PMC8280089.1.pdf data-raw/digitised   # see data-raw/digitised/README.md
```

## Sources

* Reck M, Rodríguez-Abreu D, Robinson AG, et al. Five-year outcomes with
  pembrolizumab versus chemotherapy for metastatic non–small-cell lung cancer
  with PD-L1 tumor proportion score ≥ 50%. *J Clin Oncol*
  2021;39(21):2339–2349. doi:10.1200/JCO.21.00174.
  [PMC8280089](https://pmc.ncbi.nlm.nih.gov/articles/PMC8280089). Figure 2A
  (OS) and 2B (PFS); Appendix Table A1 (subsequent therapies).
* NICE TA531, [committee discussion](https://www.nice.org.uk/guidance/ta531/chapter/3-Committee-discussion)
  and committee papers ([2](https://www.nice.org.uk/guidance/ta531/documents/committee-papers-2),
  [3](https://nice.org.uk/guidance/ta531/documents/committee-papers-3)): model
  structure, baseline characteristics, utilities, regimen mix, resource use and
  the costs that are not redacted.
* NICE DSU Technical Support Documents 14 (survival analysis), 19
  (partitioned survival analysis) and 21 (flexible survival models).
* NICE health technology evaluations: the manual (PMG36), reference case.
* Guyot P, Ades AE, Ouwens MJ, Welton NJ. *BMC Med Res Methodol* 2012;12:9.
* Office for National Statistics. National life tables: England, 2022 to 2024
  (released 10 December 2025).
* eMIT national database, 1 July 2025 to 30 June 2026.
* NHS England. National Cost Collection 2024/25 (published 6 November 2025).
* Jones KC, Weatherly H, Birch S, et al. Unit Costs of Health and Social Care
  2025 Manual (NHS Cost Inflation Index).
* dm+d / BNF: Keytruda 100 mg/4 ml, NHS indicative price.
* Hernández Alava M, Pudney S, Wailoo A. Estimating EQ-5D by age and sex for
  the UK. NICE DSU report, 2022.
* EMA. Keytruda summary of product characteristics, section 4.2.

## Licensing

The code is MIT-licensed (`LICENSE`).

Reck et al. (2021) is published under CC BY-NC-ND 4.0. This repository does
**not** reproduce the article's figures or PDF. It stores only numerical
coordinates extracted from Figure 2 and the printed numbers at risk, with
attribution in `data-raw/digitised/README.md`. ONS, eMIT and NHS England data
are used under the Open Government Licence v3.0.
