---
name: sc-format-convert
description: Convert single-cell data formats between Seurat (.rds/.qs), AnnData (.h5ad), and 10X (h5/MTX). Use when users ask to convert Seurat <-> h5ad, 10X -> Seurat, or 10X -> h5ad.
---

# sc-format-convert

Convert common single-cell file formats using lightweight scripts.

## Quick start

Seurat RDS -> h5ad:

```bash
scripts/seurat2h5ad.py input.rds -o output.h5ad
```

10X (h5 or MTX dir) -> h5ad:

```bash
scripts/tenx_to_h5ad.py <input_10x> -o output.h5ad
```

10X (h5 or MTX dir) -> Seurat RDS:

```bash
Rscript scripts/tenx_to_seurat.R <input_10x> -o output.rds
```

h5ad -> Seurat RDS:

```bash
Rscript scripts/h5ad_to_seurat.R input.h5ad -o output.rds
```

## Inputs

- **10X H5**: `.h5` (10X layout)
- **10X MTX**: directory containing `matrix.mtx(.gz)` and `barcodes.tsv(.gz)` plus `features.tsv(.gz)` or `genes.tsv(.gz)`

## Notes

## Environment / dependency check

All scripts in `sc-format-convert/scripts/` do a quick dependency check on startup.
If required R/Python packages are missing, the script will exit early and print an install hint.

## Dependencies

- **Python scripts** (`*.py`): `scanpy`, `anndata`, `numpy`, `pandas`, `scipy`
- **R scripts** (`*.R`): `Seurat`, `Matrix`
- **h5ad -> Seurat** (`scripts/h5ad_to_seurat.R`): additionally requires `schard` ([`cellgeni/schard`](https://github.com/cellgeni/schard); reticulate-free)
