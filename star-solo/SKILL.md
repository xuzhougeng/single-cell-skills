---
name: star-solo
description: Run STARsolo for single-cell RNA-seq mapping, demultiplexing, UMI deduplication, and per-cell quantification. Use when building STAR references, running 10x droplet protocols (v2/v3/5'), configuring barcode/UMI geometry, matching CellRanger outputs, filtering cells (knee or EmptyDrops_CR), handling multi-gene reads, or ingesting BAM inputs for STARsolo.
---

# STARsolo workflow

Use this skill to assemble correct STARsolo commands and reference inputs. Prefer the protocol recipes and option glossary in `references/star-solo.md` instead of re-deriving parameters.

## Quick start

1) Decide protocol and barcode geometry.
   - Choose droplet vs SmartSeq.
   - For 10x droplet, confirm chemistry (v2 vs v3) and read order (cDNA read first).

2) Build the STAR genome index with the same FASTA/GTF used by CellRanger if matching outputs.
   - Consider `--genomeSAsparseD 3` for CellRanger-like index (slower, lower RAM).

3) Run STAR with STARsolo enabled.
   - Use `--soloType CB_UMI_Simple` (or `Droplet`) for 10x.
   - Provide `--soloCBwhitelist` and set UMI length for v3.

4) Apply cell filtering (optional).
   - Use `--soloCellFilter EmptyDrops_CR` for CellRanger-like filtering.
   - Use `--runMode soloCellFiltering` for filtering existing raw matrices.

5) Validate outputs.
   - Inspect `Solo.out/` matrices and `Summary.csv`.
   - Confirm counts with expected chemistry and reference.

## Common decisions

- Match CellRanger raw counts by using CellRanger FASTA/GTF and matching UMI/cell-barcode settings.
- Match CellRanger 4/5 by enabling adapter clipping and the CellRanger 3+ barcode/UMI options.
- Enable additional features (GeneFull, SJ, Velocyto) only when needed.
- Choose a multi-gene strategy only if you need recovery of multi-mapping reads.

## References

- Use `references/star-solo.md` for protocol-specific command templates, barcode geometry, and parameter glossary.
