#!/usr/bin/env python3
"""Assemble manuscript tables from frozen results; never modify inputs.

Run with the bundled Python runtime (reportlab and pypdf are required).
Supplementary TSVs are byte-for-byte source copies. No new inference is run.
"""
from pathlib import Path
import csv
import hashlib
import json
import shutil
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "manuscript_figures_tables" / "tables"
MAIN = OUT / "main"
SUPP = OUT / "supplementary"
for p in (MAIN, SUPP):
    p.mkdir(parents=True, exist_ok=True)


def read(relative):
    with (ROOT / relative).open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows, columns=None):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns or list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


manifest = []


def track(source, target, transform):
    path = ROOT / source
    manifest.append(dict(source=source, source_sha256=sha(path), source_rows=len(read(source)),
                         output=str(target.relative_to(OUT)), transformation=transform,
                         output_sha256=sha(target)))


TRAITS = "profile_gwas_call90_maf10/inputs/primary_traits.tsv"
QUALITY = "profile_phenotype_discovery/qa/phenotype_quality_rank.tsv"
GWAS = "profile_gwas_call90_maf10/results/gwas_summary.tsv"
LABELS = "manuscript_readiness/provenance/primary_trait_reporting_labels.tsv"
traits = read(TRAITS)
quality = {r["phenotype"]: r for r in read(QUALITY)}
labels = {r["trait"]: r for r in read(LABELS)}
gwas_all = read(GWAS)
gwas = {r["trait"]: r for r in gwas_all if r["variant"] == "primary"}
assert len(traits) == len(gwas) == 6
assert len(gwas_all) == 36
assert all(int(r["n_lines"]) == 176 and int(r["eligible_variants"]) == 1428812 for r in gwas.values())
assert all(int(r["per_trait_hits"]) == int(r["six_trait_hits"]) == 0 for r in gwas_all)

short_names = {
    "NightSleepChange": "Night sleep change", "MovingShapePC1": "Moving shape PC1",
    "NightPWake": "Night P(Wake)", "NightPDoze": "Night P(Doze)",
    "LDHarmonicR2": "LD harmonic R2", "TotalSleep": "Total sleep",
}
definitions = {
    "NightSleepChange": "Night sleep on D6 minus D1 (ZT8-24)",
    "MovingShapePC1": "PC1 of mean-centered 48-bin moving profiles",
    "NightPWake": "Corrected log odds of one-minute inactive-to-active transitions (ZT8-24)",
    "NightPDoze": "Corrected log odds of one-minute active-to-inactive transitions (ZT8-24)",
    "LDHarmonicR2": "Daily moving-profile variance explained by fixed 24 h and 12 h harmonics",
    "TotalSleep": "Sleep minutes per 24 h day",
}
units = {"NightSleepChange": "min", "MovingShapePC1": "PC score", "NightPWake": "log odds",
         "NightPDoze": "log odds", "LDHarmonicR2": "fraction", "TotalSleep": "min/day"}

table1 = []
display1 = []
for t in traits:
    tr = t["trait"]
    q = quality[t["phenotype"]]
    assert int(q["eligible_lines"]) == 176
    r = dict(trait=tr, phenotype=t["phenotype"], label=short_names[tr], definition=definitions[tr],
             window=labels[tr]["window"], units=units[tr], eligible_lines=int(q["eligible_lines"]),
             eligible_line_batches=int(q["eligible_line_batches"]), repeated_lines=int(q["repeated_lines"]),
             batch_pairs=int(q["batch_pairs"]), repeatability_rho=float(q["repeatability_rho"]),
             rho_ci_low=float(q["rho_ci_low"]), rho_ci_high=float(q["rho_ci_high"]))
    table1.append(r)
    display1.append([r["label"], r["definition"], r["window"], r["units"], str(r["eligible_lines"]),
                     f'{r["repeatability_rho"]:.3f}\n({r["rho_ci_low"]:.3f}, {r["rho_ci_high"]:.3f})'])
table1_notes = [
    "Male DGRP flies were group-housed for 2-5 days at 21.5 C under 12L:12D before individual DAM housing at 21.5 C under 8L:16D. D1 and D6 are both individually housed observations. Negative night sleep change means less sleep on D6.",
    "All six traits include 176 eligible DGRP lines. Eligibility requires at least one line-batch with eight valid flies. There are 256 eligible line-batches for five traits and 253 for P(Doze). Scalar line-batch values are medians; PC scores are means. Primary GWAS inputs are batch-adjusted line estimates.",
    "Repeatability is the symmetrized Spearman correlation between repeated line-batch estimates. Parentheses show 95% percentile intervals from 200 genotype-cluster bootstrap draws. Five traits use 67 repeated lines (95 batch pairs); P(Doze) uses 65 lines (91 pairs). These are rank-repeatability estimates, not heritability.",
    "D2-D6 denotes an operational summary window, not demonstrated steady state or completed entrainment. P(Wake)/P(Doze) use one-minute activity states rather than the five-minute sleep definition. Centering retains amplitude differences; harmonic fit under LD does not test endogenous DD rhythmicity.",
]
table1_headers = ["Phenotype", "Definition", "Window", "Unit", "Lines", "Repeatability\n(95% interval)"]

table2 = []
display2 = []
for t in traits:
    tr = t["trait"]
    r = gwas[tr]
    table2.append(dict(trait=tr, phenotype_label=short_names[tr], n_lines=int(r["n_lines"]),
                       eligible_variants=int(r["eligible_variants"]), valid_tests=int(r["valid_tests"]),
                       lambda_gc=float(r["lambda_gc"]), minimum_p=float(r["minimum_p"]),
                       lead_variant=r["lead_variant"], lead_minor_carriers=int(r["lead_minor_carriers"]),
                       nominal_hits=int(r["nominal_hits"]), per_trait_hits=int(r["per_trait_hits"]),
                       six_trait_hits=int(r["six_trait_hits"]),
                       trait_threshold=float(r["trait_threshold"]), family_threshold=float(r["family_threshold"])))
    display2.append([short_names[tr], f'{float(r["minimum_p"]):.2e}', r["lead_variant"],
                     r["lead_minor_carriers"], r["nominal_hits"], f'{float(r["lambda_gc"]):.3f}', r["six_trait_hits"]])
table2_notes = [
    "Primary GCTA mixed-model scans used 176 DGRP lines and 1,428,812 eligible variants per trait: genotype call rate >=90% and MAF >=10% in the analysis cohort. Models used batch-adjusted line phenotypes, the existing genomic relationship matrix, and jointly fitted Wolbachia/inversion covariates.",
    "Lead variant is the lowest-P variant for that trait. Minor carriers counts lines carrying at least one minor allele at that variant. Nominal hits are variant-trait records with P < 1e-5, not independent loci: 50 records across six traits represent 47 unique variants.",
    "Bonferroni thresholds are 3.499e-8 per trait and 5.832e-9 across six traits. No variant passes either threshold. Family hits reports the six-trait threshold. Lambda GC is genomic inflation. These scans generate association candidates, not validated causal genes.",
]
table2_headers = ["Phenotype", "Minimum P", "Lead variant", "Minor\ncarriers", "Nominal\nhits", "Lambda\nGC", "Family\nhits"]


def pdf_table(path, title, headers, rows, widths, notes):
    body = ParagraphStyle("body", fontName="Helvetica", fontSize=10.5, leading=13, textColor=colors.HexColor("#202020"))
    header = ParagraphStyle("header", parent=body, fontName="Helvetica-Bold")
    small = ParagraphStyle("note", parent=body, fontSize=9.5, leading=12)
    heading = ParagraphStyle("heading", parent=body, fontName="Helvetica-Bold", fontSize=13, leading=17, spaceAfter=14)
    def p(v, style):
        return Paragraph(escape(str(v)).replace("\n", "<br/>"), style)
    cells = [[p(v, header) for v in headers]] + [[p(v, body) for v in row] for row in rows]
    tab = Table(cells, colWidths=widths, repeatRows=1, hAlign="LEFT")
    tab.setStyle(TableStyle([
        ("LINEABOVE", (0, 0), (-1, 0), 0.9, colors.black),
        ("LINEBELOW", (0, 0), (-1, 0), 0.6, colors.black),
        ("LINEBELOW", (0, -1), (-1, -1), 0.9, colors.black),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 4), ("RIGHTPADDING", (0, 0), (-1, -1), 7),
        ("TOPPADDING", (0, 0), (-1, -1), 8), ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
    ]))
    doc = SimpleDocTemplate(str(path), pagesize=(792, 612), rightMargin=36, leftMargin=36,
                            topMargin=31, bottomMargin=28, title=title, author="")
    story = [Paragraph(escape(title), heading), tab, Spacer(1, 14)]
    for note in notes:
        story.extend([Paragraph(escape(note), small), Spacer(1, 5)])
    doc.build(story)
    pages = PdfReader(path).pages
    assert len(pages) == 1, (path, len(pages))
    text = pages[0].extract_text()
    assert all(short_names[t["trait"]] in text for t in traits)


def md_table(path, title, headers, rows, notes):
    lines = [f"# {title}", "", "| " + " | ".join(h.replace("\n", " ") for h in headers) + " |",
             "| " + " | ".join("---" for _ in headers) + " |"]
    lines.extend("| " + " | ".join(str(v).replace("\n", " ") for v in row) + " |" for row in rows)
    lines.append("")
    for note in notes:
        lines.extend([note, ""])
    path.write_text("\n".join(lines))


for num, title, headers, display, raw, widths, notes in [
    (1, "Table 1. Primary sleep and activity phenotypes", table1_headers, display1, table1, [121, 236, 72, 64, 43, 184], table1_notes),
    (2, "Table 2. Primary GWAS summary", table2_headers, display2, table2, [148, 95, 145, 85, 82, 83, 82], table2_notes),
]:
    stem = MAIN / f"Table{num}_{'primary_phenotypes' if num == 1 else 'GWAS_summary'}"
    write_tsv(stem.with_suffix(".tsv"), raw)
    md_table(stem.with_suffix(".md"), title, headers, display, notes)
    pdf_table(stem.with_suffix(".pdf"), title, headers, display, widths, notes)
    for source in ([TRAITS, QUALITY, LABELS] if num == 1 else [TRAITS, GWAS]):
        track(source, stem.with_suffix(".tsv"), "Keyed selection of six primary traits; preserve numeric precision; concise display labels")

groups = [
    ("S1", "Cohort and recording design", [
        ("batch_inventory", "manuscript_readiness/provenance/batch_inventory.tsv"),
        ("line_batch_retention", "profile_phenotype_discovery/qa/cohort_retention_by_line_batch.tsv"),
        ("monitor_design", "manuscript_readiness/provenance/line_batch_monitor_design.tsv"),
    ], "Counts are flies, distinct DGRP lines or monitor files as named. CS is retained in the cohort/design records but excluded from DGRP line counts. Analysis timestamps retain the source's UTC encoding; they are not verified local loading or lights-on clock times."),
    ("S2", "Phenotype definitions and discovery quality", [
        ("phenotype_dictionary", "profile_phenotype_discovery/phenotype_dictionary.tsv"),
        ("phenotype_quality", QUALITY),
        ("minimum_n_sensitivity", "profile_phenotype_discovery/qa/minimum_n_sensitivity.tsv"),
    ], "The dictionary covers all 106 explored phenotype fields, not 106 primary endpoints. Unit, definition, aggregation and eligibility are explicit columns. Quality rho values are rank correlations; intervals are bootstrap bounds. Historical machine names such as State_D2_D6 are retained and must be interpreted using the current Table 1 notes."),
    ("S3", "Line-level phenotype estimates", [
        ("six_primary_line_estimates", "profile_gwas_call90_maf10/inputs/line_phenotypes.tsv"),
        ("all_batch_adjusted_line_estimates", "profile_phenotype_discovery/data/line_phenotypes_batchFE.tsv"),
        ("all_line_batch_estimates", "profile_phenotype_discovery/data/line_batch_phenotypes.tsv"),
    ], "The six-trait table contains 176 lines x six traits. raw and batchFE are in the trait units listed in Table 1. all-phenotype estimates join to S2a by phenotype for units. n_flies/total_flies are phenotype-valid fly counts, n_batches counts eligible recording batches. S3c retains all line-batch cells, including ineligible cells: descriptive eligible cells require n_flies >= 8 and a finite value. The original batch-adjustment fit may use finite cells with n_flies >= 6, including CS, before reporting eligible DGRP line estimates. Line-batch scalar values are medians; PC scores are means. Do not treat all repeated cells or flies as independent genotypes."),
    ("S4", "GWAS scans and nominal associations", [
        ("all_36_scan_summaries", GWAS),
        ("all_nominal_associations", "profile_gwas_call90_maf10/results/all_nominal_hits.tsv"),
        ("testing_thresholds", "profile_gwas_call90_maf10/qa/testing_thresholds.tsv"),
    ], "variant in scan summaries denotes analysis variant, not SNP identity. There are six primary and 30 sensitivity scans. The nominal table retains all 311 records across these scans; filter variant=primary to obtain 50 primary variant-trait records (47 unique variants) with P < 1e-5. b and se are effect/standard-error per coded A1 allele in phenotype units, except RINT scans use the rank-inverse-normalized scale. Freq/MAF/call_rate are fractions, bp is genomic base position. No scan has a Bonferroni-significant hit."),
    ("S5", "Selected candidate regions and diagnostics", [
        ("60_annotated_candidates", "profile_gwas_call90_maf10/results/primary_leads_annotated.tsv"),
        ("all_annotation_effects", "profile_gwas_call90_maf10/results/lead_annotations_all_effects.tsv"),
        ("candidate_sensitivities", "profile_gwas_call90_maf10/results/lead_sensitivity.tsv"),
        ("leave_one_line_out", "profile_gwas_call90_maf10/results/lead_leave_one_line_out.tsv"),
        ("missingness_diagnostics", "profile_gwas_call90_maf10/results/lead_missingness_diagnostics.tsv"),
    ], "Ten physically separated lead candidates per trait give 60 variant-trait records; this is a 250 kb proximity screen, not LD clumping or evidence for 60 independent loci. Only 31 of these records satisfy P < 1e-5. Coordinates and allele-aware annotations retain source fields, warnings and all effects. An annotation identifies genomic context, not causal gene validation. P-values are dimensionless; effect estimates use trait units except RINT scans, which use rank-inverse-normalized units and are not directly magnitude-comparable to original-scale effects. Carrier counts are lines, not minor-allele copies."),
    ("S6", "Sleep-window and baseline sensitivity", [
        ("common_cohort_window_agreement", "manuscript_readiness/baseline/time_window_common_fly_cohort_agreement.tsv"),
        ("window_repeatability", "manuscript_readiness/baseline/time_window_repeatability.tsv"),
        ("baseline_change_associations", "manuscript_readiness/baseline/baseline_change_associations.tsv"),
        ("independent_fly_splits", "manuscript_readiness/baseline/independent_fly_split_summary.tsv"),
        ("nonoverlapping_calendar_repeatability", "manuscript_readiness/provenance/repeatability_nonoverlapping_calendar_sensitivity.tsv"),
    ], "Sleep amounts and changes are minutes; rho values are Spearman correlations. Window-agreement comparisons use matched flies for each contrast pair. Independent-fly split percentile ranges summarize split-to-split variation, not confidence intervals. Calendar sensitivity omits the overlapping line 705 Batch3/Batch4a pair. Alternative windows do not isolate causal social effects or establish completed entrainment."),
    ("S7", "Fly-QC and candidate robustness", [
        ("cohort_counts", "manuscript_readiness/qc/cohort_counts.tsv"),
        ("phenotype_concordance", "manuscript_readiness/qc/concordance.tsv"),
        ("repeatability", "manuscript_readiness/qc/repeatability.tsv"),
        ("candidate_robustness", "manuscript_readiness/candidate_qc/candidate_robustness_summary.tsv"),
        ("candidate_comparisons", "manuscript_readiness/candidate_qc/candidate_comparisons.tsv"),
    ], "Schemes include the recorded cohort, current-rule simulation and relaxed QC alternatives. additional/lost are flies relative to the recorded cohort. rho is rank correlation; phenotype differences retain original units, while rms_difference_over_baseline_sd is standardized. Candidate robustness holds 60 selected records fixed and refits models; it is not a fresh genome-wide discovery scan. Effect ratios are relative to original estimates; P-values and nominal thresholds are as named."),
    ("S8", "PCA structure and stability", [
        ("variance_explained", "profile_phenotype_discovery/pca/variance_explained.tsv"),
        ("loadings", "profile_phenotype_discovery/pca/loadings.tsv"),
        ("scalar_correlations", "profile_phenotype_discovery/pca/scalar_interpretation_correlations.tsv"),
        ("bootstrap_stability", "manuscript_readiness/pca/pca_bootstrap_summary.tsv"),
        ("PC1_loading_intervals", "manuscript_readiness/pca/PC1_loading_bootstrap_intervals.tsv"),
        ("conditional_selection", "manuscript_readiness/pca/conditional_selection_summary.tsv"),
    ], "PCA uses 48 half-hour profile bins (zt is hours after lights-on). variance_fraction is a fraction; loadings/PC scores use their stored PCA scaling. Centered-shape PCA removes each profile mean but retains amplitude. Bootstrap percentile bounds have the resampling unit named in each row. Conditional-selection frequencies concern the original automatic discovery panel, not a probability that one of the six current primary traits is biologically true."),
]

notes = ["# Table notes and source index", "", "Two compact main tables accompany eight supplementary table groups. Main PDF tables use Helvetica, white backgrounds and three horizontal rules. TSV files retain full numeric precision and flat headers. No source data or analysis results are overwritten.", "", "## Main tables", "", "- Table 1: six primary phenotype definitions, units, eligibility and rank repeatability.", "- Table 2: six primary GWAS scans. Detailed candidate rows remain supplementary.", "", "Interpretation follows the confirmed 12L:12D group-housing to 8L:16D individual-housing design at unchanged 21.5 C. D6-D1 is a within-individual-housing comparison, not a net isolation effect. D2-D6 is an operational summary window. No other laboratory manuscript's raw data are included.", ""]
supp_counts = {}
for number, title, files, note in groups:
    notes.extend([f"## Supplementary Table {number}. {title}", "", note, ""])
    for i, (name, source) in enumerate(files):
        code = number + chr(ord("a") + i)
        target = SUPP / f"Table{code}_{name}.tsv"
        shutil.copyfile(ROOT / source, target)
        track(source, target, "Byte-for-byte copy; no filtering or rounding")
        assert sha(ROOT / source) == sha(target)
        supp_counts[code] = len(read(source))
        notes.append(f"- [{code}: {name.replace('_', ' ')}](supplementary/{target.name}) - {supp_counts[code]:,} rows.")
    notes.append("")
notes.extend(["## Provenance and checks", "", "`TABLE_SOURCE_MANIFEST.tsv` records exact source paths, SHA-256 hashes, row counts and transformations. Supplementary sources are copied byte-for-byte. The two main TSV tables are keyed projections from existing estimates and contain unrounded numeric values. Display PDFs/Markdown round rho to three decimals and P-values to three significant digits.", "", "All retained batch labels are preserved. Sample sizes differ by phenotype; the whole retained cohort is not the denominator of every analysis. Column names and historical machine labels are not silently renamed. Review the unit/definition columns in S2a together with Table 1 and these notes.", ""])
(OUT / "TABLE_NOTES.md").write_text("\n".join(notes))
write_tsv(OUT / "TABLE_SOURCE_MANIFEST.tsv", manifest)

nominal = [r for r in read("profile_gwas_call90_maf10/results/all_nominal_hits.tsv") if r["variant"] == "primary"]
leads = read("profile_gwas_call90_maf10/results/primary_leads_annotated.tsv")
primary_lines = read("profile_gwas_call90_maf10/inputs/line_phenotypes.tsv")
assert len(nominal) == 50 and len({r["SNP"] for r in nominal}) == 47
assert len(leads) == 60 and sum(r["nominal_screen"] == "TRUE" for r in leads) == 31
assert len(primary_lines) == 1056
assert all(len({r["Genotype"] for r in primary_lines if r["trait"] == t["trait"]}) == 176 for t in traits)
checks = dict(main_tables=2, main_pdf_pages=[1, 1], supplementary_groups=len(groups),
              supplementary_files=sum(len(x[2]) for x in groups), supplementary_rows=supp_counts,
              source_copy_hashes_match=True, primary_traits=6, primary_lines_per_trait=176,
              primary_eligible_variants=1428812, all_scans=36, primary_nominal_records=len(nominal),
              unique_primary_nominal_variants=47, physically_screened_candidate_records=len(leads),
              nominal_records_in_candidate_subset=31, bonferroni_hits_all_scans=0,
              pdf_text_labels_verified=True, visual_review="Required after rendering; not asserted by this script")
(OUT / "TABLE_VALIDATION.json").write_text(json.dumps(checks, indent=2) + "\n")
print(json.dumps(checks, indent=2))
