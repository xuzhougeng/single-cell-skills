# Python workflows (AnnData / 10X)

Use Python for `.h5ad` and optionally 10X raw inputs.

## AnnData (.h5ad)

```python
import anndata as ad
import numpy as np
import scipy.sparse as sp

adata = ad.read_h5ad(path, backed="r")
```

Core fields:
- n_cells: `adata.n_obs`
- n_genes: `adata.n_vars`
- layers: `list(adata.layers.keys())`
- metadata fields: `list(adata.obs.columns)`
- reductions: `list(adata.obsm.keys())`
- clusters: `obs` columns like `leiden` or `louvain` if present

Sparsity:
```python
X = adata.X
nnz = X.nnz if sp.issparse(X) else np.count_nonzero(X)
sparsity = 1 - nnz / (adata.n_obs * adata.n_vars)
```

Basic per-cell QC:
```python
if sp.issparse(X):
    n_count = np.asarray(X.sum(axis=1)).ravel()
    n_feature = np.asarray((X > 0).sum(axis=1)).ravel()
else:
    n_count = X.sum(axis=1)
    n_feature = (X > 0).sum(axis=1)
```

## 10X H5

```python
import scanpy as sc
adata = sc.read_10x_h5(path)
```

## 10X MTX

```python
adata = sc.read_10x_mtx(path, var_names="gene_symbols")
```

For raw 10X, only compute basic counts and sparsity.
