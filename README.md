# SEER ampullary morphology analysis

R statistical code for the analysis of detailed morphology coding, clinical
selection, and cancer-death cumulative incidence in ampullary adenocarcinoma.
This is a **code-only** repository. It contains no manuscript, patient records,
analysis dataset, fitted patient-level objects, table values, or figure files.

## Authorized input

Obtain your own SEER authorization and accept the applicable data-use agreement:
<https://seer.cancer.gov/data/>. No data access is granted by this repository.
Use SEER*Stat 9.1.0.2: Incidence - SEER Research Data, 8 Registries,
Nov 2025 Sub (1975-2023), released April 2026. Select malignant C24.1,
known age 20+, diagnosis years 1975-2023; do not restrict histology,
diagnostic confirmation, first primary, or survival completeness. Export labels,
not numeric representations, as UTF-8 tab-separated text with the original column
names and header. The original broad export has 60 columns; see `docs/INPUT.md`.
The earlier reference-only export is not required; code reconstructs its definition.
Exact original-cohort count assertions intentionally detect a different input extract.
Later data releases may not reproduce the published extraction.

## Software and execution

Tested analysis environment: R 4.6.1, with survival, cmprsk, splines, and parallel.
Package versions from the original analyses are in `docs/SOFTWARE.md`.
Nothing installs R or packages, accesses a network, or uploads results.

From the repository directory, with **both input and output outside it**:

```sh
Rscript run_analysis.R /private/authorized_export.txt /private/analysis_results full
```

On Windows, quote paths containing spaces. Full mode preserves original random
seeds, 1,000 patient-clustered bootstrap attempts for risk/probability analyses,
and 300 attempts for the older nested descriptive contrast. Full execution may
take considerable time. Supporting composition models use four local R workers.

```sh
Rscript run_analysis.R /private/authorized_export.txt /private/smoke_results smoke
```

Smoke mode uses five resamples and is **not** a publication-analysis mode. Repeated
execution overwrites only generated files in the chosen output directory; use a
new output directory to retain an earlier run. To run selected existing modules:

```sh
Rscript run_analysis.R /private/authorized_export.txt /private/analysis_results full 6,8,9,10
```

Module 1 must have already created the private derived RDS. Module 5 uses module
4's estimator definitions; module 6 also loads those definitions without running
module 4's analysis loop. Module 10 needs completed modules 6 and 8.

After a full run, check the locally generated primary intervals and independent
estimator checks (these tests read private outputs but do not publish them):

```sh
Rscript tests/check_results.R /private/analysis_results
python tests/check_distribution.py
```

The Python check uses only the standard library and is optional for statistical
execution; it audits the code distribution for prohibited files and sensitive
patterns. It does not install anything.

## Modules and analysis coverage

| Module | Role | Manuscript result mapping |
|---|---|---|
| 1 | Definitions, eligibility, descriptive cohorts, sequential cause-specific Cox models | Table 1; secondary supplementary results |
| 2 | Overlap diagnostics and interval-specific secondary hazards | Supplementary methods/results |
| 3 | Original full-window overlap-weighted risks and bootstrap diagnostics | Tables S1-S2 |
| 4 | Comparator, era-specific grade, full censoring-model diagnostics, influence analysis | Table S3; supplementary diagnostics |
| 5 | Reduced censoring models for reported IPCW sensitivity | Table S3 |
| 6 | Supported three-year primary risks, exclude-8255, strict surgery, era comparisons | Table 3; Tables S1-S2 |
| 7 | Supporting hypothetical composition standardization and descriptive subcode risks | Tables S4-S5 |
| 8 | Stage/surgery/morphology audit and modern coding probabilities | Table 4; Tables S6-S7 |
| 9 | Independent formula/predict check of coding probabilities | Computational verification |
| 10 | TIFF plot of coding coverage, probabilities and primary risk contrasts | Figure 1 |
| 11 | Follow-up, missingness and exclusion summaries | Results and Table 2 context |
| 12 | Modern cohort characteristics, era denominators and flow | Table 2 |

Table 2 clinical characteristics are computed by modules 11-12; no document-building
code or manuscript tables are included.
Table S8 is a narrative reporting application, not a statistical model output.
All analysis scopes were exploratory; the code archive is not preregistration.

## Interpretation

Detailed codes are 8144/8163/8255; NOS is 8140. Code 8255 does not uniformly
establish a mixed intestinal-pancreatobiliary phenotype. U.S. coding introduced
8163 for 2018 diagnoses and previously directed pancreatobiliary-type cases to
8255. Summary Stage distant is not equivalent to AJCC M1. Analyses target
different populations; their contrasts are not sequential causal decomposition.
Conditional coding probabilities are not a validated prediction model.
Surgical restriction does not establish subtype equivalence or surgery benefit.
Bootstrap assesses internal uncertainty, not independent external validation.
Full-window late-era horizons, nonsurgical contrasts and hypothetical composition
standardization have sparse support; preserve warnings and failures, not just CIs.
The unadjusted 8144/8163 comparison also has unequal diagnosis-era composition.

## Privacy and outputs

The runner refuses input/output inside this repository. Generated RDS files include
restricted patient records and identifiers needed only for local clustered resampling.
Aggregate crosswalks can contain small cells. **Do not redistribute any generated
output without a separate SEER disclosure review.** `.gitignore` is a second layer,
not permission to commit restricted files. Session reports and error logs may contain
local paths; they also remain outside this repository.

The public distribution is assembled using an explicit file allowlist. See
`tests/check_distribution.py`. Only source code and technical documentation are
distributed. Author emails, institutions, user-directory paths and patient values
are excluded. Git and citation metadata use the selected repository account only.

## Citation and license

Use the versioned GitHub release and the Zenodo DOI once it has been issued.
No DOI is asserted before successful archival. Citation metadata are supplied in
`CITATION.cff`. MIT applies to this code only; it grants no rights to SEER data.
