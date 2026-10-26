# Source checks

Each fact used for validation or model structure is checked against the source
document and quoted exactly. These facts are **not** model inputs, unless an
entry in `data-raw/inputs.csv` cites them separately.

Sources checked on 2026-10-23:

- **[Reck 2021]** Reck M, Rodríguez-Abreu D, Robinson AG, et al. Five-year
  outcomes with pembrolizumab versus chemotherapy for metastatic non–small-cell
  lung cancer with PD-L1 tumor proportion score ≥ 50%. *J Clin Oncol*
  2021;39(21):2339–2349. doi:10.1200/JCO.21.00174. PMC8280089. Data cut-off
  1 June 2020.
- **[TA531]** NICE. Pembrolizumab for untreated PD-L1-positive metastatic
  non-small-cell lung cancer. Technology appraisal guidance TA531, section 3
  (Committee discussion).
  <https://www.nice.org.uk/guidance/ta531/chapter/3-Committee-discussion>.
  Paragraph numbers below are the guidance paragraph numbers.

Status key: ✅ confirmed as stated in the brief; ⚠️ confirmed, with a nuance
noted.

## KEYNOTE-024 five-year analysis (Reck 2021)

| # | Fact | Status | Exact quote | Location |
|---|---|---|---|---|
| R1 | 154 patients randomised to pembrolizumab and 151 to chemotherapy | ✅ | "Three hundred five patients were randomly assigned: 154 to pembrolizumab and 151 to chemotherapy." | Abstract, Results |
| R2 | Median OS 26.3 vs 13.4 months; HR 0.62 (95% CI 0.48–0.81) | ✅ | "Median OS was 26.3 months (95% CI, 18.3 to 40.4 months) in the pembrolizumab group and 13.4 months (95% CI, 9.4 to 18.3 months) in the chemotherapy group (HR, 0.62; 95% CI, 0.48 to 0.81; Fig 2A)." | Results, Efficacy |
| R3 | 5-year OS 31.9% vs 16.3% | ✅ | "Kaplan-Meier estimates of OS at 5 years were 31.9% (95% CI, 24.5 to 39.5) in the pembrolizumab group and 16.3% (95% CI, 10.6 to 23.0) in the chemotherapy group." | Results, Efficacy |
| R4 | Median PFS 7.7 vs 5.5 months; HR 0.50 (95% CI 0.39–0.65) | ✅ | "Median PFS was 7.7 months (95% CI, 6.1 to 10.2 months) for patients in the pembrolizumab group versus 5.5 months (95% CI, 4.2 to 6.2 months) for patients in the chemotherapy group (HR, 0.50; 95% CI, 0.39 to 0.65; Fig 2B)." | Results, Efficacy |
| R5 | 3- and 5-year PFS (landmarks used under D-003) | ✅ | "Kaplan-Meier estimates of PFS at 3 and 5 years were 22.8% (95% CI, 16.3 to 29.9) and 12.8% (95% CI, 7.4 to 19.8), respectively, in the pembrolizumab group and 4.1% (95% CI, 1.3 to 9.4) and not reached (no evaluable patients) in the chemotherapy group." | Results, Efficacy |
| R6 | Effective crossover rate 66.0% | ✅ | "Among patients initially assigned to chemotherapy, 99 received anti‒PD-1 or PD-L1 (PD-[L]1) therapy (83 crossed over to pembrolizumab on-study, and 16 additional patients received anti‒PD-(L)1 therapy outside the study), representing a 66.0% effective crossover rate." | Results, Patients |
| R7 | OS event totals, same data cut (used under D-004) | ✅ | "At data cutoff, 103 patients (66.9%) in the pembrolizumab group and 123 patients (81.5%) in the chemotherapy group had died." | Results, Efficacy |
| R8 | PFS event totals, same data cut (used under D-004) | ✅ | "At data cutoff, 126 patients (81.8%) in the pembrolizumab group and 141 patients (93.4%) in the chemotherapy group had progressed per investigator assessment or died." | Results, Efficacy |
| R9 | Data cut-off date | ✅ | "Median time from randomization to data cutoff (June 1, 2020) was 59.9 months (range, 55.1-68.4 months)." | Results, Patients |
| R10 | The published HRs are stratified (relevant to D-005) | ✅ | "A stratified Cox model with Efron's method of tie handling was used to determine HRs and 95% CIs; randomization stratification factors were applied to the Cox models." | Methods, Statistical analysis |
| R11 | Stratification factors | ✅ | "Random assignment was stratified by Eastern Cooperative Oncology Group performance status (0 or 1), histology (squamous or nonsquamous), and region (East Asian or non–East Asian enrollment)." | Methods |
| R12 | Trial dosing and 35-cycle maximum | ✅ | "Patients were randomly assigned (1:1) to receive 200 mg of intravenous pembrolizumab once every 3 weeks for 35 cycles or investigator's choice of platinum chemotherapy (carboplatin/cisplatin) plus pemetrexed or gemcitabine, or carboplatin plus paclitaxel for 4 to 6 cycles." | Methods |
| R13 | Licence | ✅ | "Creative Commons Attribution Non-Commercial No Derivatives 4.0 License: https://creativecommons.org/licenses/by-nc-nd/4.0/" | Article footer |
| R14 | What the figures show (for digitising) | ✅ | "Kaplan-Meier estimates of (A) OS and (B) PFS in the pembrolizumab group and the chemotherapy group." … "Tick marks indicate censoring times." | Figure 2 caption |

## NICE TA531 committee discussion

| # | Fact | Status | Exact quote | Paragraph |
|---|---|---|---|---|
| N1 | 20-year time horizon | ✅ | "The ERG highlighted that although the final overall survival analysis from the company includes an additional 6 months of data, this only accounts for 10% of the 20‑year time horizon." | 3.10 |
| N2 | 2-year (35-cycle) maximum treatment duration | ✅ | "The committee was aware that the maximum possible treatment duration with pembrolizumab in KEYNOTE‑024 was 2 years (35 cycles)." | 3.6 |
| N3 | 2-year stopping rule in the model | ✅ | "The committee concluded that implementing a 2‑year stopping rule in the model was appropriate." | 3.8 |
| N4 | KM to 33 weeks, then exponential; 23 and 43 weeks as alternatives | ⚠️ | "To estimate overall survival, the company's base-case analysis used 33‑week Kaplan–Meier data from KEYNOTE‑024. After 33 weeks, the company fitted separate exponential models to the data. The company's base-case incremental cost-effectiveness ratio (ICER) using this approach was £30,244 per quality-adjusted life year (QALY) gained. In scenario analyses the company explored alternative time points to extrapolate the trial data; 23 and 43 weeks of Kaplan–Meier data." **Nuance:** the guidance describes this approach for **OS** only. It does not state the PFS extrapolation method. | 3.9 |
| N5 | ICERs for the 23- and 43-week cut-offs | ✅ | "Both these analyses increased the company's base-case ICER; to £31,321 per QALY gained for the 23‑week time point and £33,829 per QALY gained for the 43‑week time point." | 3.9 |
| N6 | Waning: OS HR set to 1 at 3 and 5 years | ✅ | "The committee noted the company's scenario analyses, which explored stopping the treatment effect by setting the overall survival hazard ratio at 1 at different time points in the model. This increased the company's base-case ICER from £30,244 per QALY gained to £44,483 and £36,156 per QALY gained for the 3- and 5‑year scenarios." | 3.11 |
| N7 | Data cut used in TA531 | ✅ | "…the tail of the pembrolizumab Kaplan–Meier overall survival curve based on the latest KEYNOTE‑024 data (July 2017)." | 3.11 |
| N8 | Utilities based on time to death | ⚠️ | "The utility values for pembrolizumab and standard care were pooled (adjusted for age) and divided into 4 groups based on time to death (from less than 30 days to at least 360 days)." **Nuance:** the committee preferred the company's values with the ≥360-days group set to the UK population norm: "…preferred to use values from KEYNOTE‑024 combined with the ERG's approach of setting the utility value for 360 days to death to the UK population norm." (3.14). The values themselves are not given in the guidance. | 3.13, 3.14 |
| N9 | Company base-case ICER £30,244 per QALY, including a 2-year stopping rule | ✅ | "It was aware that the company's base-case ICER was £30,244 per QALY gained (including a 2‑year stopping rule)." | 3.16 |
| N10 | Committee range about £30,000 to under £50,000 | ✅ | "However, it agreed that the ICERs, ranging from £30,000 per QALY gained to less than £50,000 per QALY gained, on which it was basing its decision were associated with uncertainty, particularly for overall survival and the duration of treatment effect." | 3.16 |
| N11 | ICERs include confidential discounts | ⚠️ | "It noted that the company's model included the agreed commercial access agreement discount for pembrolizumab and an estimated discount for pemetrexed." (3.7); "The decision-making ICERs cannot be presented because they include the commercial access agreement reduction for pemetrexed." (3.16). **Nuance:** the discounts are for pembrolizumab **and** pemetrexed. | 3.7, 3.16 |
| N12 | Modelled survival 36.0 vs 22.3 months | ⚠️ | "The committee heard that the average number of months of life gained with pembrolizumab, as estimated by the company's economic model, is 36.0 months compared with 22.3 months for standard care." **Nuance:** this is mean modelled survival per arm. The guidance does not say whether it is discounted, but values used for the end-of-life criteria are conventionally undiscounted. | 3.18 |
| N13 | 200 mg vs 2 mg/kg dosing issue | ⚠️ | "The ERG highlighted that in the economic model, the company assumed the fixed 200 mg dose of pembrolizumab for people in the standard care group (who go on to have pembrolizumab after disease progression), which is not in line with the marketing authorisation or how it is given in NHS clinical practice." **Nuance:** the issue concerns **subsequent** pembrolizumab in the chemotherapy arm, not first-line dosing. | 3.15 |
| N14 | No crossover adjustment | ✅ | "The committee agreed that no adjustment was needed to account for treatment switching in the trial." | 3.5 |
| N15 | Follow-up in the data used by TA531 | ✅ | "The median duration of follow‑up was 25.2 months; 14.9% of people remained on pembrolizumab compared with 1.3% on standard care." | 3.5 |

## NICE TA531 committee papers

Source: NICE TA531 committee papers
(<https://www.nice.org.uk/guidance/ta531/documents/committee-papers-3>),
downloaded 2026-10-23. Page numbers are PDF page numbers.

| # | Fact | Status | Exact quote | Location |
|---|---|---|---|---|
| C1 | Model starting age and sex mix (used under D-012) | ✅ | "Table 52. Baseline characteristics of patients included in the model … Average age 65 - KEYNOTE-024 CSR; Proportion male 64.6% - KEYNOTE-024 CSR" | Company submission, Table 52, p245 |
| C2 | The same values in the ERG report | ✅ | "Table 35 Model baseline patient characteristics … Age 65 years KEYNOTE-024 trial; Proportion of male patients 64.6% KEYNOTE-024 trial" | ERG report, Table 35, p495 |
| C3 | Three-state partitioned survival structure | ✅ | "The model consisted of three health states: pre-progression, post-progression and death (see Figure 26)." | Company submission §5.2.2, p245 |
| C4 | Weekly cycle with half-cycle correction | ✅ | "A cycle length of one week was considered sufficient to reflect the patterns of treatment administration and the transitions to disease progression and death. In line with previous submissions, a half-cycle correction was implemented to mitigate bias." | Company submission §5.2.2, p245 |
| C5 | Pembrolizumab dosing in the company model | ✅ | "Pembrolizumab was implemented in the model as per the anticipated licensed dosing regimen, i.e. administered intravenously at a fixed dose of 200mg over 30 minutes every 3 weeks." | ERG report §5.3.5, p495 |
| C6 | Time-to-death and progression-based utilities (used under D-019) | ✅ | Table 62: "≥360 … 0.808 0.021 (0.767, 0.850); [180, 360) … 0.712 0.025; [30, 180) … 0.598 0.026; <30 … 0.48 0.075" (pooled). Table 63: "Progression-Free … 0.778 0.008 … Progressive … 0.668 0.02" (pooled) | Company submission Tables 62–63, p269 |
| C7 | Chemotherapy regimen mix | ⚠️ | Table 69 (KEYNOTE-024, all): "Gem + Car 13.3% … Gem + Cis 7.3% … Pac + Car 11.3% … Pemx + Cis 44.0% … Pemx + Car 24.0%". **Nuance:** Table 54 (p251) gives pemetrexed/carboplatin 44% and pemetrexed/cisplatin 24%, the reverse. Table 69 is used (D-020). | Company submission Table 69, p284 |
| C8 | Pemetrexed maintenance | ✅ | "In the chemotherapy group, one patient withdrew consent before receiving planned study treatment (Figure 6), and 46 patients received pemetrexed maintenance therapy (Figure 6)." | Company submission, p155 |
| C9 | Maximum chemotherapy cycles and dose intensity | ✅ | "A maximum treatment duration of 18 weeks (i.e., 6 cycles for the platinum-based therapies administrated every 3 weeks) was used…"; "on average, 99.21% of patients on pembrolizumab and 97.05% of patients on overall platinum-based chemotherapy received their planned doses." | Company submission, p286 |
| C10 | Weekly monitoring costs and terminal care | ✅ | "The estimated per week monitoring and disease management costs were £76.75 and £125.87 respectively for the PFS and PPS periods." (p290–291); "The estimated one-off terminal costs were £4,735.73" (p292) | Company submission §5.5.6 |
| C11 | Subsequent therapy durations | ✅ | "The duration of the second line treatment for docetaxel is assumed to be 3 cycles (i.e., 9 weeks) and 8.7 cycles (i.e., 26.1 weeks) for pembrolizumab based on data observed in KEYNOTE-024." | Company submission §5.3.6, p265 and p289 |
| C12 | Crossover proportion used by the company (2017 data) | ✅ | "For patients in the SOC arm, patients are assumed to receive pembrolizumab based on the proportion of patients who crossed over in KEYNOTE-024 (43.7%), with the rest of the patients assumed to receive docetaxel." | Company submission §5.3.6, p265 |

## Other sources

| # | Fact | Status | Exact quote | Location |
|---|---|---|---|---|
| M1 | Reference-case discount rate | ✅ | "For the reference case, costs and health effects should be discounted at the same rate of 3.5% per year." | NICE health technology evaluations manual (PMG36), §4.5.1 |
| S1 | Current pembrolizumab dosing | ✅ | The recommended dose of KEYTRUDA in adults is either 200 mg every 3 weeks or 400 mg every 6 weeks (paraphrase of section 4.2, as returned by a web search of the EMA document; not quoted directly). | Keytruda SmPC (EMA product information), section 4.2 |
| R15 | Subsequent therapies by arm (used under D-022) | ✅ | Table A1: "With one or more subsequent oncologic treatment 77 (50.0) / 103 (68.2)"; "Immunotherapy 20 (13.0) / 99 (65.6)"; "Platinum doublet chemotherapy with or without bevacizumab 57 (37.0) / 8 (5.3)"; "Single-agent chemotherapy with or without antiangiogenic agent 31 (20.1) / 43 (28.5)" | Reck 2021, Appendix Table A1 |

## Not yet checked

None outstanding for Step 3.
