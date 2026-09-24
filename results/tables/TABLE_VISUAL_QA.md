# Table visual verification

Both final one-page main PDFs were rendered with Poppler at 1600 pixels on the long edge and visually inspected on 2026-09-11. The tables contain all six traits, readable body text and notes, three horizontal rules, no vertical cell borders, and no clipped or overlapping text. The exported text was also checked for all six phenotype labels.

Numeric verification is recorded in `TABLE_VALIDATION.json`. Supplementary files are flat machine-readable TSVs, not typeset tables, and their byte-for-byte source integrity is recorded in `TABLE_SOURCE_MANIFEST.tsv`.

Poppler required a temporary local fontconfig file pointing to macOS system fonts and a writable temporary cache. This affected rendering only, not the PDF contents or analysis results.
