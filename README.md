# SEER ampullary morphology statistical analysis code

R code for *Detailed morphology coding preferentially captures surgical patients with ampullary adenocarcinoma in SEER*.

The primary analysis quantifies the association of recorded surgery with detailed morphology coding in SEER 17, 2018–2022. SEER 8 provides historical coverage and supporting analyses. Exploratory competing-risk analyses examine the applicability of registry subtype evidence. This distribution includes analysis and figure generators, variable definitions, and a table-to-output map. It contains no individual-level SEER data, fitted patient objects, manuscript, or result datasets.

## Data and software
Obtain independently authorized labeled UTF-8, tab-separated, 60-column SEER Research Case Listing exports from the November 2025 submission: SEER 8 (1975–2023) and SEER 17 (2000–2023, selecting 2018–2022). Select malignant C24.1 tumors and known age 20+, without histology, microscopic-confirmation, first-primary, or survival restrictions. See [input fields](seer8/docs/INPUT.md) and [variable algorithms](docs/VARIABLES.md).

Requires R 4.6.1, survival 3.8-6, cmprsk 2.2-12, and the base splines/parallel packages. Figure scripts require Cairo TIFF support and Arial or a system substitute. Keep input data and generated outputs outside the repository.

## Execution
Run from this directory, replacing the private example paths:

```sh
Rscript analysis/analyses.R /private/SEER17.txt /private/results/seer17 seer17
Rscript analysis/confirmation_sensitivity.R /private/SEER17.txt /private/results/confirmation
Rscript seer8/run_analysis.R /private/SEER8.txt /private/results/seer8 full 1,2,3,4,5,6,7,8,9,11,12
Rscript analysis/analyses.R /private/SEER8.txt /private/results/seer8_coding
Rscript figures/figure1.R /private/results/seer17 /private/results/confirmation /private/results/figures
Rscript figures/figure2.R /private/SEER8.txt /private/SEER17.txt /private/results/seer17 /private/results/figures
Rscript seer8/tests/check_results.R /private/results/seer8
python tools/verify_manifest.py
```

SEER8 module 10 optionally generates diagnostic plots; the two scripts in figures/ generate the manuscript figures at 600 dpi. Numerical bootstrap analyses use 1,000 patient-clustered resamples and fixed seeds recorded in each script. The historical runner also supports smoke mode (five resamples) for software checks; smoke intervals must never be used for inference. Full bootstrap execution can be computationally intensive. Warnings, support diagnostics, and failed fits are recorded by the analyses.

[MANUSCRIPT_MAP.csv](docs/MANUSCRIPT_MAP.csv) maps Tables 1–3, S1–S19, and Figures 1–2 to their numerical generators. Table 1 contains code definitions and S8 contains reporting recommendations; neither is a fitted statistical result. The scripts generate numerical sources, not Word document layout.

## Interpretation
Observed coding coverage and the surgery proportion within coding groups use different denominators. The primary standardized association averages predictions over the 2,535 complete-case reference records, not all 2,675 reference records, and is not a causal surgery effect. SEER 17 overlaps geographically with SEER 8. Registry labels are not pathology-validated. Exploratory survival models with different overlap weights target different populations; their contrast cannot isolate surgery adjustment, and wide intervals do not establish subtype equivalence. The confirmation sensitivity adds only 15 reference records and does not exclude broader diagnostic selection.

## Citation and license
Use CITATION.cff and the version-specific Zenodo archive associated with the GitHub release. MIT license. SEER access and use remain subject to the applicable data-use agreement.
