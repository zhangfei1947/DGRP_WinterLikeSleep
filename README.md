# DGRP WinterLikeSleep

Data and reusable analysis code for **Natural variation in sleep and activity profiles during social isolation under short photoperiod in Drosophila melanogaster**, by Fei Zhang, Esther Doria, and Wanhe Li.

Release **v1.0.0** contains the experimental DAM recording windows, metadata and recorded QC decisions, processed phenotypes, genotype inputs and population-structure matrix, all 36 genome-wide association scans, and the source data for Figures 1–5, Figures S1–S6, Tables 1–2, and Supplementary Tables S1–S8. The primary analyses contain 176 eligible DGRP lines and 6,228 retained flies including CS controls. Association candidates are exploratory; no association passed the stated Bonferroni thresholds.

## Download

Clone this repository, then download the versioned archives from [release v1.0.0](https://github.com/zhangfei1947/DGRP_WinterLikeSleep/releases/tag/v1.0.0). Extract all archives in the repository root. The archives restore the directory structure expected by the code. Alternatively, with Python 3.12 or later:

```sh
python3 code/download_data.py
python3 code/verify_files.py
Rscript code/check_environment.R
```

`ASSETS.json` records download URLs, file sizes, and SHA-256 hashes. `MANIFEST.tsv` records the extracted files. GitHub's automatically generated source-code ZIP does **not** include the large data archives. For figure reproduction alone, obtain `reproduction_inputs.tar.gz` and the six `gwas_*.tar.gz` archives; DAM files and BED genotypes are needed for earlier stages.

## Contents

| Directory | Contents |
| --- | --- |
| `data/raw_dam/` | 30 batch folders, 286 cropped monitor files, loading metadata, original retention summaries and exclusion lists |
| `data/metadata/` | Batch and monitor inventories and exact cropping/checksum manifest |
| `reproducibility/publication_analysis/` | Frozen minute-derived measures, recorded QC, genotype BED/BIM/FAM, and the original scaled GRM |
| `reproducibility/profile_phenotype_discovery/` | Individual and line phenotypes, half-hour profiles, PCA models/loadings, and phenotype QC |
| `reproducibility/profile_gwas_call90_maf10/` | Model inputs, complete association scans, allele counts, candidate summaries, and sensitivity results |
| `reproducibility/manuscript_readiness/` | Time-window, baseline, PCA, and fly-QC sensitivity results used in the manuscript |
| `reproducibility/scripts/` | Selected original analysis/figure scripts, with portable raw-data defaults |
| `code/` | Entry points for downloading, verifying, cropping, checking dependencies, running GWAS, and exporting eligible statistics |
| `results/` | Machine-readable manuscript tables and figure source data |
| `docs/` | Reproduction instructions, data/GWAS definitions, software versions, and release validation |

## Reproduce the manuscript figures and tables

From the repository root after extracting data:

```sh
cd reproducibility
Rscript scripts/manuscript_profile_figures.R
Rscript scripts/check_manuscript_profile_figures.R
Rscript scripts/manuscript_gwas_figures.R
python3 scripts/manuscript_tables.py
```

Outputs are written to `reproducibility/manuscript_figures_tables/`. These commands rebuild all 11 figures and the manuscript tables from the provided frozen inputs. Python table rendering requires `reportlab` and `pypdf`. R dependencies and recorded versions are described in `docs/REPRODUCTION.md` and the session-information files.

Earlier stages, optional GWAS reruns, and the scope of the selected code are described in [Reproduction instructions](docs/REPRODUCTION.md). Raw-window and QC conventions are described in [Data definitions](docs/DATA.md); the association mask, coded alleles and GRM are described in [GWAS inputs and outputs](docs/GWAS.md).

## Reuse and attribution

Original analysis code and documentation in this repository are distributed under the [MIT license](LICENSE). The MIT license does not relicense third-party software or the published DGRP resource. Cite this study when using its recordings or phenotypes and cite Mackay et al. (2012; doi:10.1038/nature10811) and Huang et al. (2014; doi:10.1101/gr.171546.113) when using the DGRP genotype resource. Third-party resources retain their original terms. No manuscript DOI has been assigned in this release.

ChatGPT (OpenAI; versions 5.6 and 6.0) assisted analysis planning, statistical-method discussions, code development/debugging, execution, visualization, literature searches, interpretation, and manuscript preparation under author direction. The authors checked the code, results, cited references, and manuscript text and retain responsibility for the work.
