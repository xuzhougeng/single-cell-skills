#!/usr/bin/env python3

import argparse
import json
import os
import sys

import numpy as np
import scipy.sparse as sp


def detect_10x_mtx_dir(path):
    if not os.path.isdir(path):
        return False
    matrix_files = ["matrix.mtx", "matrix.mtx.gz"]
    barcode_files = ["barcodes.tsv", "barcodes.tsv.gz"]
    feature_files = ["features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"]
    has_matrix = any(os.path.exists(os.path.join(path, f)) for f in matrix_files)
    has_barcodes = any(os.path.exists(os.path.join(path, f)) for f in barcode_files)
    has_features = any(os.path.exists(os.path.join(path, f)) for f in feature_files)
    return has_matrix and has_barcodes and has_features


def summarize_numeric(x):
    if x is None or len(x) == 0:
        return None
    return {
        "min": float(np.min(x)),
        "median": float(np.median(x)),
        "mean": float(np.mean(x)),
        "max": float(np.max(x)),
    }


def pick_expression_matrix(adata):
    """
    Pick a reasonable expression matrix for QC/count summaries.

    Many .h5ad store counts in layers (e.g. layers['counts']) and leave adata.X=None.
    In backed mode, some matrices may be disk-backed wrappers; this function only
    selects the matrix, conversion is handled later.
    """
    if getattr(adata, "X", None) is not None:
        return adata.X, "X"

    if getattr(adata, "raw", None) is not None and getattr(adata.raw, "X", None) is not None:
        return adata.raw.X, "raw.X"

    layers = getattr(adata, "layers", None)
    if layers is not None and len(layers.keys()) > 0:
        # Common conventions in the wild
        for key in ("counts", "raw_counts", "count", "umi", "UMI"):
            if key in layers:
                return layers[key], f"layers[{key}]"
        key = next(iter(layers.keys()))
        return layers[key], f"layers[{key}]"

    raise ValueError(
        "No expression matrix found: adata.X is None and no adata.raw.X / adata.layers present."
    )


def _coerce_to_memory(X):
    """
    Convert disk-backed wrappers to in-memory array/sparse matrix when needed.
    """
    # anndata's sparse_dataset.SparseDataset exposes to_memory()
    if hasattr(X, "to_memory") and callable(getattr(X, "to_memory")):
        return X.to_memory()
    return X


def compute_counts(X):
    X = _coerce_to_memory(X)

    if sp.issparse(X):
        nnz = int(X.nnz)
        n_count = np.asarray(X.sum(axis=1)).ravel()
        n_feature = np.asarray((X > 0).sum(axis=1)).ravel()
        n_cells, n_genes = X.shape
    else:
        X = np.asarray(X)
        if X.ndim != 2:
            raise ValueError(
                f"Expression matrix must be 2D (cells x genes); got ndim={X.ndim}."
            )
        nnz = int(np.count_nonzero(X))
        n_count = np.asarray(X.sum(axis=1)).ravel()
        n_feature = np.asarray((X > 0).sum(axis=1)).ravel()
        n_cells, n_genes = X.shape

    denom = n_cells * n_genes
    sparsity = 1 - (nnz / denom if denom else 0.0)
    return nnz, sparsity, n_count, n_feature, n_cells, n_genes


def render_table(out):
    lines = [
        "========================================",
        "sc-info (Python)",
        "========================================",
        f"Input: {out['input_path']}",
        f"Type: {out['input_type']}",
        f"Cells: {out['counts']['n_cells']}",
        f"Genes: {out['counts']['n_genes']}",
        f"Nonzero: {out['counts']['n_nonzero']}",
        f"Sparsity: {out['counts']['sparsity']:.4f}",
    ]
    if out["assays"]:
        lines.append(f"Assays: {', '.join(out['assays'])}")
    if out["layers"]:
        lines.append(f"Layers: {', '.join(out['layers'])}")
    if out["metadata_fields"]:
        lines.append(f"Metadata fields: {len(out['metadata_fields'])}")
    if out["reductions"]:
        lines.append(f"Reductions: {', '.join(out['reductions'])}")
    if out["clusters"]:
        lines.append(f"Clusters: {', '.join(out['clusters'])}")
    lines.append("========================================")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Summarize single-cell dataset info (AnnData/10X)."
    )
    parser.add_argument("input_path")
    parser.add_argument("--json-out", dest="json_out")
    parser.add_argument("--md-out", dest="md_out")
    parser.add_argument("--table-out", dest="table_out")
    args = parser.parse_args()

    input_path = args.input_path
    if not os.path.exists(input_path):
        print(f"Input does not exist: {input_path}", file=sys.stderr)
        return 1

    ext = os.path.splitext(input_path)[1].lower().lstrip(".")
    adata = None
    input_type = None

    if ext == "h5ad":
        import anndata as ad

        adata = ad.read_h5ad(input_path, backed="r")
        input_type = "h5ad"
    elif ext == "h5":
        import scanpy as sc

        adata = sc.read_10x_h5(input_path)
        input_type = "10x_h5"
    elif detect_10x_mtx_dir(input_path):
        import scanpy as sc

        adata = sc.read_10x_mtx(input_path, var_names="gene_symbols")
        input_type = "10x_mtx"
    elif ext in ("rds", "qs"):
        print("Use the R script for .rds/.qs inputs.", file=sys.stderr)
        return 1
    else:
        print(f"Unsupported input type: {ext}", file=sys.stderr)
        return 1

    X, matrix_source = pick_expression_matrix(adata)
    nnz, sparsity, n_count, n_feature, n_cells, n_genes = compute_counts(X)

    qc = {
        "n_count": summarize_numeric(n_count),
        "n_feature": summarize_numeric(n_feature),
    }

    metadata_fields = list(adata.obs.columns)
    layers = list(adata.layers.keys())
    reductions = list(adata.obsm.keys())
    clusters = []
    for key in ("leiden", "louvain"):
        if key in adata.obs.columns:
            clusters.append(key)
    assays = []

    percent_mt = None
    for key in ("percent.mt", "pct_counts_mt", "pct_counts_mito", "percent_mito"):
        if key in adata.obs.columns:
            percent_mt = summarize_numeric(adata.obs[key].to_numpy())
            break
    if percent_mt is not None:
        qc["percent_mt"] = percent_mt

    output = {
        "input_path": input_path,
        "input_type": input_type,
        "engine": "Python",
        "matrix_source": matrix_source,
        "counts": {
            "n_cells": int(n_cells),
            "n_genes": int(n_genes),
            "n_nonzero": int(nnz),
            "sparsity": float(sparsity),
        },
        "qc": qc,
        "assays": assays,
        "layers": layers,
        "metadata_fields": metadata_fields,
        "reductions": reductions,
        "clusters": clusters,
    }

    table_text = render_table(output)
    print(table_text)

    if args.table_out:
        with open(args.table_out, "w", encoding="utf-8") as f:
            f.write(table_text + "\n")
    if args.json_out:
        with open(args.json_out, "w", encoding="utf-8") as f:
            json.dump(output, f, indent=2)
    if args.md_out:
        md = []
        md.append("# sc-info")
        md.append("")
        md.append(f"- Input: {output['input_path']}")
        md.append(f"- Type: {output['input_type']}")
        md.append(f"- Engine: {output['engine']}")
        md.append("")
        md.append("## Counts")
        md.append("")
        md.append(f"- Cells: {output['counts']['n_cells']}")
        md.append(f"- Genes: {output['counts']['n_genes']}")
        md.append(f"- Nonzero: {output['counts']['n_nonzero']}")
        md.append(f"- Sparsity: {output['counts']['sparsity']:.4f}")
        if layers:
            md.append("")
            md.append("## Layers")
            md.append("")
            md.append(f"- {', '.join(layers)}")
        if metadata_fields:
            md.append("")
            md.append("## Metadata fields")
            md.append("")
            md.append(f"- {', '.join(metadata_fields)}")
        if reductions:
            md.append("")
            md.append("## Reductions")
            md.append("")
            md.append(f"- {', '.join(reductions)}")
        if clusters:
            md.append("")
            md.append("## Clusters")
            md.append("")
            md.append(f"- {', '.join(clusters)}")
        with open(args.md_out, "w", encoding="utf-8") as f:
            f.write("\n".join(md) + "\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
