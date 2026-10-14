# Validation boundaries

The archive is a portability refactor of the analysis scripts, not a new analysis.
The import accepts the all-malignant export directly rather than requiring a
second duplicate reference export. Inputs and outputs are user-supplied private
paths. Original full-mode seeds, model specifications and estimators are retained.
No manuscript-writing, author-metadata or rendering scripts are included.

The R scripts assert the original extraction counts, stage/surgery consistency,
model-matrix rank, valid risk calculations and independent weighted
Aalen-Johansen agreement with survival::survfit. Module 9 independently checks
the coding-probability model through formula-based glm and predict.

Local release validation uses the authorized data outside this repository to
compare reproduced outputs against the frozen original outputs. Neither those
data nor results are deposited here. A five-resample smoke run only verifies
execution; it is not a substitute for full bootstrap intervals.

## Version 1.0.0 release checks

On 2026-10-15 all twelve modules completed locally in full mode with R 4.6.1.
Twenty designated CSV outputs matched the original frozen outputs exactly after
CSV parsing. The comparisons covered descriptive and secondary hazards,
original and supported-window overlap risks, grade/comparator sensitivity,
reduced IPCW, composition standardization, descriptive subcode risks, and the
coding audit/probability model, including saved bootstrap draws where applicable.
Independent AJ/model checks and the distribution scanner passed. The runner also
refused a test attempting to write analysis output inside the repository.
These are computational reproducibility checks, not external scientific validation.

Known limitations: boundary fits and failed resamples are retained and counted;
unsupported nonsurgical contrasts must not be used as scientific evidence;
the composition model relies on sparse-cell extrapolation; the subtype comparison
mixes diagnosis eras; none is external validation or a causal estimate.
