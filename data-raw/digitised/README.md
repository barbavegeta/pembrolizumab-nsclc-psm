# Digitised Kaplan–Meier data: KEYNOTE-024 five-year analysis

## Source and attribution

The coordinates in this folder were extracted from Figure 2 of:

> Reck M, Rodríguez-Abreu D, Robinson AG, et al. Five-year outcomes with
> pembrolizumab versus chemotherapy for metastatic non–small-cell lung cancer
> with PD-L1 tumor proportion score ≥ 50%. *J Clin Oncol* 2021;39(21):2339–2349.
> doi:10.1200/JCO.21.00174. Open access at
> <https://pmc.ncbi.nlm.nih.gov/articles/PMC8280089>.

- Figure 2A: overall survival (OS), intention-to-treat population.
- Figure 2B: progression-free survival (PFS), investigator review.
- Data cut-off: 1 June 2020.

The article is published under the Creative Commons Attribution
Non-Commercial No Derivatives 4.0 licence
(<https://creativecommons.org/licenses/by-nc-nd/4.0/>). The repository does
not reproduce the figure, any crop of it, any image derived from it, or the
article PDF. It holds only the numerical coordinates of the plotted curves
and the numbers at risk printed beneath them.

## Method

Figure 2 is embedded in the article PDF as vector artwork. Each arm's
Kaplan–Meier curve is a single stroked path in its own colour (pembrolizumab
blue, RGB 0/0.463/0.659; chemotherapy red, RGB 0.671/0.086/0.173). The
coordinates are therefore read directly from the path geometry rather than
traced from an image. Each vertex of the path is a corner of a KM step, and
every step is captured.

`../extract_km_vectors.py` does the extraction:

1. Checks the PDF against the PMC Open Access copy (md5
   `4b2d15722a5c6f1228124ce57c1c78f2`).
2. Calibrates each panel separately. It fits a least-squares line through
   the tick marks drawn on each axis: one tick per month on the x-axis, and
   ticks from 10% to 100% on the y-axis. The maximum residual is 0.011
   months for x and below 0.0001 for survival.
3. Writes the vertices of each arm's step path, in drawing order, as
   `time_months,surv`.
4. Censoring tick marks are separate paths. They are not used, because the
   Guyot algorithm infers censoring from the numbers at risk.

The numbers at risk were read from the PDF's text layer and match the
printed tables.

To reproduce:

```
curl -o PMC8280089.1.pdf https://pmc-oa-opendata.s3.amazonaws.com/PMC8280089.1/PMC8280089.1.pdf
python data-raw/extract_km_vectors.py PMC8280089.1.pdf data-raw/digitised
```

This requires Python 3 with PyMuPDF (2.0) and numpy. Keep the downloaded
PDF out of the repository.

## Files

| File | Content | Columns |
|---|---|---|
| `os_pembro.csv` | Fig 2A, pembrolizumab curve | `time_months,surv` |
| `os_chemo.csv` | Fig 2A, chemotherapy curve | `time_months,surv` |
| `pfs_pembro.csv` | Fig 2B, pembrolizumab curve | `time_months,surv` |
| `pfs_chemo.csv` | Fig 2B, chemotherapy curve | `time_months,surv` |
| `nar_os.csv` | Fig 2A numbers at risk, every 6 months from 0 to 72 | `time_months,pembro,chemo` |
| `nar_pfs.csv` | Fig 2B numbers at risk, every 6 months from 0 to 66 | `time_months,pembro,chemo` |
| `axis_limits.csv` | Last labelled x-axis tick of each panel (OS 72, PFS 66) | `endpoint,x_max_months` |

- `time_months` is time in months.
- `surv` is survival as a proportion between 0 and 1.
- The curve files are inputs to `R/digitise_clean.R`. Every cleaning step it
  applies is logged in `output/validation/reconstruction.md`.
