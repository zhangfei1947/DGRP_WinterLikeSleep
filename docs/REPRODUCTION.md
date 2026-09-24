# Reproduction instructions

Run commands from the repository root unless a `cd` is shown. Extract the assets first. The original analysis used R 4.4.2, PLINK 1.9.0-b.7.11 and GCTA 1.95.3. R packages include data.table, damr, behavr, sleepr, ggplot2, ggrepel, patchwork and lme4; recorded session information accompanies this file. Install versions appropriate to your operating system. Executables are not redistributed.

## Figures and tables

Use the four commands in the main README. They use frozen processed inputs and do not rerun GWAS. All 11 manuscript figures can be regenerated. The generated figure source TSVs can be compared with `results/source_data/`; visual rasterization can differ with fonts or graphics-library versions.

## Raw DAM to minute and individual sleep measures

```sh
cd reproducibility
Rscript scripts/rebuild_minute_phenotypes.R ../data/raw_dam ../rerun/minute_rebuild '^Batch[[:alnum:]]+_DGRP$' true
```

The four principal output files are `fly_qc_all_batches.csv`, `daily_qc_all_batches.csv`, `window_metrics_all_batches.csv` and `individual_phenotypes_all_batches.csv`. The release validation rebuilt all 30 batches from the cropped files and compared these files against the original full-history analysis. Primary fly retention is read from the original GIGEM summaries; this is deliberate and preserves the published cohort.

## Profiles, PCA and line estimates

```sh
cd reproducibility
Rscript scripts/extract_profile_phenotypes.R ../data/raw_dam ../rerun/profile_phenotypes
Rscript scripts/analyze_profile_phenotypes.R ../rerun/profile_phenotypes
```

These commands use the supplied frozen minute-window metrics and retention decisions, reconstruct half-hour profiles from the raw records, and regenerate PCA and batch-adjusted line phenotypes. The extraction independently checks sleep and coverage agreement. To rebuild every stage consecutively, first compare the raw-minute output with the frozen inputs; then use the validated files in a working copy. The default extraction never silently replaces the supplied frozen cohort.

Time-window and PCA sensitivity scripts are also supplied:

```sh
cd reproducibility
Rscript scripts/readiness_baseline.R
Rscript scripts/readiness_pca.R
```

These write their corresponding sensitivity directories. Run in a working copy if you want to preserve the release files for checksum verification. Original primary estimates, alternate-QC results, and the final supplementary tables are included even where the exploratory development code is outside this curated package.

## GWAS

Set `GCTA_BIN` to an installed GCTA 1.95.3 executable; `GWAS_THREADS` defaults to four. The commands below preserve the supplied scaled GRM, numerical covariates and allele coding:

```sh
GCTA_BIN=/path/to/gcta64 Rscript code/run_primary_gwas.R primary
# Optional: one trait only
GCTA_BIN=/path/to/gcta64 Rscript code/run_primary_gwas.R primary NightSleepChange
# Optional alternative model; use one of the values below
GCTA_BIN=/path/to/gcta64 Rscript code/run_primary_gwas.R chromosome_LOCO
Rscript code/export_eligible_gwas.R
```

Other model names are `raw_batch`, `RINT`, `exclude_influential` and `legacy_residual`. New model outputs go into `rerun/gwas/`; the wrapper refuses to overwrite an existing result. The eligibility export validates and summarizes the supplied frozen scans and writes to `rerun/eligible_gwas/`. Apply the same masks to rerun scans for comparisons. Floating-point details can depend on the platform and GCTA build.

## Provenance and scope

The code was selected to reproduce the main data transformations and manuscript outputs without bundling every development script or abandoned analysis. Release copies change only raw-data default paths and, for the minute reconstruction, cap data.table at two threads. The profile-analysis provenance list omits a development-only test script that is not part of this release. New wrapper scripts expose model parameters and input paths. Recorded provenance paths were shortened to remove machine-specific home directories; scientific values were not edited. The original working directory and raw histories remain unchanged.

See `validation.json` for checks actually performed for this release. File checksums establish the identity of the shared files; they do not replace scientific validation or independent replication.
