# Variable algorithms and analysis populations

The exact historical fields and grade algorithms are preserved in
../seer8/docs/INPUT.md. Both inputs are tumor-record Case Listing exports;
Patient ID is used only locally for clustered resampling. Patient+record keys
must be unique. No first-primary restriction is applied.

| Variable | Exact field or rule | Missing handling |
|---|---|---|
| Site/behavior | Primary Site - labeled = C24.1-Ampulla of Vater; Behavior recode for analysis = Malignant | Enforced at import |
| Diagnosis year | Year of diagnosis | Primary SEER17 2018–2022; historical SEER8 1975–2023 |
| Age | Age recode with single ages and 85+ | Known age >=20; 85+ top-coded at 85; fixed natural spline, 3 df |
| Reference morphology | Histologic Type ICD-O-3: 8140,8144,8163,8210,8211,8255,8260:8263,8480,8481,8490 | Same set with/without microscopic confirmation |
| Confirmation | Diagnostic Confirmation: positive histology; positive cytology; histology with immunophenotyping/genetic studies; microscopic method unspecified | Required in primary; removed only in confirmation sensitivity |
| Detailed coding | Broad 8144/8163/8255; narrower 8144/8163 | 8255 remains in reference denominator when not an outcome |
| Comparator | 8140 (NOS); secondary all remaining reference codes | NOS is unspecified morphology, not a verified phenotype |
| Stage | Combined Summary Stage with Expanded Regional Codes (2004+) | Localized; 3 regional labels combined; distant; other/blank unknown |
| Surgery | RX Summ--Surg Prim Site (1998-2022): 00 No;10–90 Yes | Other/missing Unknown. Primary coding uses this definition in both datasets |
| Surgery sensitivity | SEER17 removal 30/40/60; SEER8 strict resection20–49 or60 | Different clinical restrictions; neither identifies a named operation or curative intent |
| Sex | Sex, as recorded | Male/female field |
| Race/origin | Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic) | NHAIAN and unknown combined as Other/unknown |
| Survival time | Survival months | Negative or>=9999 invalid;0 assigned0.5 months |
| Event | Vital status recode (study cutoff used), SEER cause-specific death classification | Alive0;cancer death1;other death2;unclassified missing |
| Coding model eligibility | Reference cohort, known surgery and stage | Stage sensitivities retain unknown as specified; no survival requirement |
| Survival eligibility | Valid time and classifiable event | Known stage/surgery additionally required for adjusted SEER17 comparisons |

## Models and estimands
Primary SEER17 coding: logistic outcome detailed coding; surgery + linear year +
stage + age ns3 + sex + race/origin. Standardization sets surgery for every
model record and averages fitted probabilities over that same cohort. Each
sensitivity standardizes over its own analysis population. Patient-clustered
percentile bootstrap1000; fixed spline basis; refit each sample. SEER8 2010–2022
uses diagnosis era2010–2017/2018–2022; modern SEER8 uses linear year.

Prognosis: overlap-weighted Aalen–Johansen cancer cumulative incidence with
competing other deaths; weighted Kaplan–Meier all-cause risk. Common SEER17
three-year window2018–2020; one-year2018–2022. SEER8 three-year2010–2019.
Bootstrap models are refitted. No multiplicity adjustment. Propensity weighting
targets a model-specific overlap population, not a fixed reference population.

## Interpreting the primary proportion and contrast
Coding coverage: detailed-code count / reference count within a surgery stratum.
Surgery composition: surgery count / all records within a coding group.
The primary standardization target is the 2,535-record known-stage, known-surgery
model cohort. The 24.5 percentage-point contrast describes its fitted coding
association, not all-reference coverage and not a causal effect of surgery.
The unrestricted-confirmation reference cohort adds only 15 records (2,690 versus
2,675), limiting what this sensitivity can establish about diagnostic selection.
The exploratory prognosis models with/without surgery use the same 1,392 records
but different overlap targets; their change cannot be assigned to surgery alone.
