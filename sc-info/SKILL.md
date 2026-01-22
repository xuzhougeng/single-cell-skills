---
name: sc-info
description: Summarize single-cell dataset information across formats (.rds/.qs Seurat or SummarizedExperiment, .h5ad AnnData, 10X .h5, 10x MTX). Use when users want a quick overview of cells/genes, QC summaries, assays, metadata fields, reductions/clusters, sparsity, and output as terminal table, JSON, or Markdown with optional file save.
---

# sc-info

Summarize single-cell dataset info with minimal loading and consistent outputs.

## Workflow

1) Identify input type by extension or directory contents.
2) Choose engine:
   - R: Seurat or SummarizedExperiment .rds/.qs
   - Python: AnnData .h5ad
   - 10X raw (h5 or MTX): either R or Python
3) Load data with the lightest reader available.
4) Collect summary fields (see "Outputs").
5) Print terminal table, and optionally emit JSON and Markdown to file.

## Quick start

R (Seurat/SummarizedExperiment/10X):

```bash
scripts/sc_info.R <input_path> [--json-out out.json] [--md-out out.md] [--table-out out.txt]
```

Python (AnnData/10X):

```bash
scripts/sc_info.py <input_path> [--json-out out.json] [--md-out out.md] [--table-out out.txt]
```

## Input detection

- **Seurat/SummarizedExperiment**: `.rds` or `.qs`
- **AnnData**: `.h5ad`
- **10X H5**: `.h5` with 10X layout
- **10X MTX**: directory containing `matrix.mtx` (or `.gz`) plus `barcodes.tsv` and `features.tsv` or `genes.tsv`

## Outputs

For Seurat/SummarizedExperiment/AnnData:
- **Counts**: n_cells, n_genes, n_nonzero, sparsity
- **QC summaries**: per-cell nCount/nFeature (or equivalent), percent.mt if present
- **Assays/layers**: list assays (Seurat) or layers (AnnData)
- **Metadata fields**: obs/meta column names
- **Reductions/clusters**: PCA/UMAP/tSNE availability and cluster labels (e.g., Seurat `seurat_clusters`, AnnData `leiden` or `louvain`)

For raw 10X inputs (h5 or MTX), only report:
- n_cells, n_genes, n_nonzero, sparsity
- Basic per-cell nCount/nFeature if cheap to compute
Do not infer reductions, clusters, or extra metadata.

## Output formats

- **Terminal table**: a compact summary for quick inspection
- **JSON**: structured output for downstream tooling
- **Markdown**: human-readable report with short sections

When saving to file, keep the same content as stdout.

## References

- R workflows and field mapping: `references/r.md`
- Python workflows and field mapping: `references/python.md`
