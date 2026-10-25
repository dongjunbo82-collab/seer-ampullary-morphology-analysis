# Input schema and variable algorithms

Provide the original 60-column all-malignant labeled SEER*Stat export. Extra raw
columns are retained privately but no raw column values are in the distribution.
These columns are accessed by the code (literal field names, **not patient values**):

- Patient ID (used locally for clustering, never distributed)
- Record number recode (used locally for unique tumor-record keys)
- Year of diagnosis
- Age recode with single ages and 85+
- Sex
- Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)
- Primary Site - labeled
- Histologic Type ICD-O-3
- Behavior recode for analysis
- Diagnostic Confirmation
- SEER historic stage A (1973-2015)
- Combined Summary Stage with Expanded Regional Codes (2004+)
- RX Summ--Surg Prim Site (1998-2022)
- Survival months
- Vital status recode (study cutoff used)
- SEER cause-specific death classification
- Grade Recode (thru 2017)
- Grade Clinical (2018+)
- Grade Pathological (2018+)
- Derived AJCC M, 7th ed (2010-2015)
- Derived EOD 2018 M Recode (2018+)
- Derived EOD 2018 T Recode (2018+)

Reference morphology: 8140, 8144, 8163, 8210-8211, 8255, 8260-8263,
8480-8481, 8490. Expanded sensitivity adds 8160 and 8500 without a
microscopic-confirmation restriction. Reference cohort requires positive histology,
positive exfoliative cytology, histology plus immunophenotyping/genetic studies,
or microscopic confirmation with method unspecified; literal labels are in module 1.
Detailed = reference with 8144/8163/8255; narrow excludes 8255; NOS = 8140.

Age 85+ maps to 85. Race/origin American Indian/Alaska Native and unknown map to
other/unknown. Modern stage uses the combined field, not pooled AJCC editions.
Before 2004 the historical stage is descriptive only. Surgery 00 maps to no,
10-90 to recorded surgery, otherwise unknown; do not apply the exported field
to 2023. Strict resection sensitivity uses 20-49 or 60 without assuming cure.
Survival months below zero or >=9999 are invalid; zero maps to 0.5 months.
Status: 0 living/censored, 1 SEER-classified cancer death, 2 other-cause death.
Unclassifiable events are excluded only from survival analyses.

Official resources:
- https://seer.cancer.gov/data/
- https://seer.cancer.gov/icd-o-3/sitetype.icdo3.20220429.pdf
- https://seer.cancer.gov/archive/icd-o-3/errata.sitetype.20180117.pdf
- https://staging.seer.cancer.gov/eod_public/input/3.2/ampulla_vater/ss2018/
- https://seer.cancer.gov/archive/manuals/2022/AppendixC/Surgery_Codes_Other_Sites_2022.pdf
