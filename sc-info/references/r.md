# R workflows (Seurat / SummarizedExperiment / 10X)

Use R for `.rds`/`.qs` and optionally 10X raw inputs.

## Seurat (.rds/.qs)

```r
obj <- readRDS(path) # or qs::qread for .qs
class(obj)
```

Core fields:
- n_cells: `ncol(obj)`
- n_genes: `nrow(obj)`
- assays: `Assays(obj)` or `names(obj@assays)`
- metadata fields: `colnames(obj@meta.data)`
- reductions: `Reductions(obj)` or `names(obj@reductions)`
- clusters: `obj@meta.data$seurat_clusters` if present
- counts matrix: `GetAssayData(obj, slot = "counts")`

Sparsity:
```r
counts <- GetAssayData(obj, slot = "counts")
n_nonzero <- Matrix::nnzero(counts)
sparsity <- 1 - n_nonzero / (nrow(counts) * ncol(counts))
```

Basic per-cell QC:
```r
nCount <- Matrix::colSums(counts)
nFeature <- Matrix::colSums(counts > 0)
```

If percent.mt exists:
```r
if ("percent.mt" %in% colnames(obj@meta.data)) {
  mt <- obj@meta.data$percent.mt
}
```

## SummarizedExperiment (.rds/.qs)

```r
obj <- readRDS(path)
counts <- SummarizedExperiment::assay(obj)
n_cells <- ncol(counts)
n_genes <- nrow(counts)
metadata_fields <- colnames(SummarizedExperiment::colData(obj))
```

## 10X H5

```r
counts <- Seurat::Read10X_h5(path)
if (is.list(counts)) counts <- counts[[1]]
```

## 10X MTX

```r
counts <- Seurat::Read10X(path)
if (is.list(counts)) counts <- counts[[1]]
```

For raw 10X, only compute basic counts and sparsity.
