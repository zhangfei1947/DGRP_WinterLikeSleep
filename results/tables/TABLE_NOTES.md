# Table notes and source index

Two compact main tables accompany eight supplementary table groups. Main PDF tables use Helvetica, white backgrounds and three horizontal rules. TSV files retain full numeric precision and flat headers. No source data or analysis results are overwritten.

## Main tables

- Table 1: six primary phenotype definitions, units, eligibility and rank repeatability.
- Table 2: six primary GWAS scans. Detailed candidate rows remain supplementary.

Interpretation follows the confirmed 12L:12D group-housing to 8L:16D individual-housing design at unchanged 21.5 C. D6-D1 is a within-individual-housing comparison, not a net isolation effect. D2-D6 is an operational summary window. No other laboratory manuscript's raw data are included.

## Supplementary Table S1. Cohort and recording design

Counts are flies, distinct DGRP lines or monitor files as named. CS is retained in the cohort/design records but excluded from DGRP line counts. Analysis timestamps retain the source's UTC encoding; they are not verified local loading or lights-on clock times.

- [S1a: batch inventory](supplementary/TableS1a_batch_inventory.tsv) - 30 rows.
- [S1b: line batch retention](supplementary/TableS1b_line_batch_retention.tsv) - 286 rows.
- [S1c: monitor design](supplementary/TableS1c_monitor_design.tsv) - 286 rows.

## Supplementary Table S2. Phenotype definitions and discovery quality

The dictionary covers all 106 explored phenotype fields, not 106 primary endpoints. Unit, definition, aggregation and eligibility are explicit columns. Quality rho values are rank correlations; intervals are bootstrap bounds. Historical machine names such as State_D2_D6 are retained and must be interpreted using the current Table 1 notes.

- [S2a: phenotype dictionary](supplementary/TableS2a_phenotype_dictionary.tsv) - 106 rows.
- [S2b: phenotype quality](supplementary/TableS2b_phenotype_quality.tsv) - 106 rows.
- [S2c: minimum n sensitivity](supplementary/TableS2c_minimum_n_sensitivity.tsv) - 318 rows.

## Supplementary Table S3. Line-level phenotype estimates

The six-trait table contains 176 lines x six traits. raw and batchFE are in the trait units listed in Table 1. all-phenotype estimates join to S2a by phenotype for units. n_flies/total_flies are phenotype-valid fly counts, n_batches counts eligible recording batches. S3c retains all line-batch cells, including ineligible cells: descriptive eligible cells require n_flies >= 8 and a finite value. The original batch-adjustment fit may use finite cells with n_flies >= 6, including CS, before reporting eligible DGRP line estimates. Line-batch scalar values are medians; PC scores are means. Do not treat all repeated cells or flies as independent genotypes.

- [S3a: six primary line estimates](supplementary/TableS3a_six_primary_line_estimates.tsv) - 1,056 rows.
- [S3b: all batch adjusted line estimates](supplementary/TableS3b_all_batch_adjusted_line_estimates.tsv) - 17,517 rows.
- [S3c: all line batch estimates](supplementary/TableS3c_all_line_batch_estimates.tsv) - 30,316 rows.

## Supplementary Table S4. GWAS scans and nominal associations

variant in scan summaries denotes analysis variant, not SNP identity. There are six primary and 30 sensitivity scans. The nominal table retains all 311 records across these scans; filter variant=primary to obtain 50 primary variant-trait records (47 unique variants) with P < 1e-5. b and se are effect/standard-error per coded A1 allele in phenotype units, except RINT scans use the rank-inverse-normalized scale. Freq/MAF/call_rate are fractions, bp is genomic base position. No scan has a Bonferroni-significant hit.

- [S4a: all 36 scan summaries](supplementary/TableS4a_all_36_scan_summaries.tsv) - 36 rows.
- [S4b: all nominal associations](supplementary/TableS4b_all_nominal_associations.tsv) - 311 rows.
- [S4c: testing thresholds](supplementary/TableS4c_testing_thresholds.tsv) - 1 rows.

## Supplementary Table S5. Selected candidate regions and diagnostics

Ten physically separated lead candidates per trait give 60 variant-trait records; this is a 250 kb proximity screen, not LD clumping or evidence for 60 independent loci. Only 31 of these records satisfy P < 1e-5. Coordinates and allele-aware annotations retain source fields, warnings and all effects. An annotation identifies genomic context, not causal gene validation. P-values are dimensionless; effect estimates use trait units except RINT scans, which use rank-inverse-normalized units and are not directly magnitude-comparable to original-scale effects. Carrier counts are lines, not minor-allele copies.

- [S5a: 60 annotated candidates](supplementary/TableS5a_60_annotated_candidates.tsv) - 60 rows.
- [S5b: all annotation effects](supplementary/TableS5b_all_annotation_effects.tsv) - 382 rows.
- [S5c: candidate sensitivities](supplementary/TableS5c_candidate_sensitivities.tsv) - 360 rows.
- [S5d: leave one line out](supplementary/TableS5d_leave_one_line_out.tsv) - 10,560 rows.
- [S5e: missingness diagnostics](supplementary/TableS5e_missingness_diagnostics.tsv) - 60 rows.

## Supplementary Table S6. Sleep-window and baseline sensitivity

Sleep amounts and changes are minutes; rho values are Spearman correlations. Window-agreement comparisons use matched flies for each contrast pair. Independent-fly split percentile ranges summarize split-to-split variation, not confidence intervals. Calendar sensitivity omits the overlapping line 705 Batch3/Batch4a pair. Alternative windows do not isolate causal social effects or establish completed entrainment.

- [S6a: common cohort window agreement](supplementary/TableS6a_common_cohort_window_agreement.tsv) - 9 rows.
- [S6b: window repeatability](supplementary/TableS6b_window_repeatability.tsv) - 18 rows.
- [S6c: baseline change associations](supplementary/TableS6c_baseline_change_associations.tsv) - 6 rows.
- [S6d: independent fly splits](supplementary/TableS6d_independent_fly_splits.tsv) - 24 rows.
- [S6e: nonoverlapping calendar repeatability](supplementary/TableS6e_nonoverlapping_calendar_repeatability.tsv) - 6 rows.

## Supplementary Table S7. Fly-QC and candidate robustness

Schemes include the recorded cohort, current-rule simulation and relaxed QC alternatives. additional/lost are flies relative to the recorded cohort. rho is rank correlation; phenotype differences retain original units, while rms_difference_over_baseline_sd is standardized. Candidate robustness holds 60 selected records fixed and refits models; it is not a fresh genome-wide discovery scan. Effect ratios are relative to original estimates; P-values and nominal thresholds are as named.

- [S7a: cohort counts](supplementary/TableS7a_cohort_counts.tsv) - 7 rows.
- [S7b: phenotype concordance](supplementary/TableS7b_phenotype_concordance.tsv) - 72 rows.
- [S7c: repeatability](supplementary/TableS7c_repeatability.tsv) - 42 rows.
- [S7d: candidate robustness](supplementary/TableS7d_candidate_robustness.tsv) - 60 rows.
- [S7e: candidate comparisons](supplementary/TableS7e_candidate_comparisons.tsv) - 360 rows.

## Supplementary Table S8. PCA structure and stability

PCA uses 48 half-hour profile bins (zt is hours after lights-on). variance_fraction is a fraction; loadings/PC scores use their stored PCA scaling. Centered-shape PCA removes each profile mean but retains amplitude. Bootstrap percentile bounds have the resampling unit named in each row. Conditional-selection frequencies concern the original automatic discovery panel, not a probability that one of the six current primary traits is biologically true.

- [S8a: variance explained](supplementary/TableS8a_variance_explained.tsv) - 288 rows.
- [S8b: loadings](supplementary/TableS8b_loadings.tsv) - 864 rows.
- [S8c: scalar correlations](supplementary/TableS8c_scalar_correlations.tsv) - 1,494 rows.
- [S8d: bootstrap stability](supplementary/TableS8d_bootstrap_stability.tsv) - 52 rows.
- [S8e: PC1 loading intervals](supplementary/TableS8e_PC1_loading_intervals.tsv) - 192 rows.
- [S8f: conditional selection](supplementary/TableS8f_conditional_selection.tsv) - 51 rows.

## Provenance and checks

`TABLE_SOURCE_MANIFEST.tsv` records exact source paths, SHA-256 hashes, row counts and transformations. Supplementary sources are copied byte-for-byte. The two main TSV tables are keyed projections from existing estimates and contain unrounded numeric values. Display PDFs/Markdown round rho to three decimals and P-values to three significant digits.

All retained batch labels are preserved. Sample sizes differ by phenotype; the whole retained cohort is not the denominator of every analysis. Column names and historical machine labels are not silently renamed. Review the unit/definition columns in S2a together with Table 1 and these notes.
