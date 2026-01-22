#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
10X (h5 or MTX) to h5ad converter.
"""

import argparse
import os
import sys

import scanpy as sc


def _is_mtx_dir(path):
    mtx = os.path.join(path, "matrix.mtx")
    mtx_gz = os.path.join(path, "matrix.mtx.gz")
    barcodes = os.path.join(path, "barcodes.tsv")
    barcodes_gz = os.path.join(path, "barcodes.tsv.gz")
    features = os.path.join(path, "features.tsv")
    features_gz = os.path.join(path, "features.tsv.gz")
    genes = os.path.join(path, "genes.tsv")
    genes_gz = os.path.join(path, "genes.tsv.gz")

    return (
        (os.path.exists(mtx) or os.path.exists(mtx_gz))
        and (os.path.exists(barcodes) or os.path.exists(barcodes_gz))
        and (
            os.path.exists(features)
            or os.path.exists(features_gz)
            or os.path.exists(genes)
            or os.path.exists(genes_gz)
        )
    )


def _infer_output_path(input_path, output_path):
    if output_path:
        return output_path
    if os.path.isdir(input_path):
        base = os.path.basename(os.path.abspath(input_path))
        return f"{base}.h5ad"
    base, _ = os.path.splitext(input_path)
    return f"{base}.h5ad"


def main():
    ap = argparse.ArgumentParser(
        description="Convert 10X h5 or MTX directory to h5ad.",
    )
    ap.add_argument("input", help="10X .h5 file or MTX directory")
    ap.add_argument("-o", "--output", help="Output h5ad path")
    ap.add_argument(
        "--var-names",
        choices=["gene_symbols", "gene_ids"],
        default="gene_symbols",
        help="Which column to use as gene names for MTX input (default: gene_symbols)",
    )
    ap.add_argument(
        "--make-genes-unique",
        action="store_true",
        help="Make gene names unique",
    )
    ap.add_argument(
        "--genome",
        default=None,
        help="Genome name for 10X h5 (optional)",
    )
    ap.add_argument(
        "--include-all-assays",
        action="store_true",
        help="Include non-GEX assays for 10X h5 (sets gex_only=False)",
    )
    args = ap.parse_args()

    input_path = args.input
    if not os.path.exists(input_path):
        sys.stderr.write(f"[error] Input not found: {input_path}\n")
        sys.exit(1)

    if os.path.isfile(input_path) and input_path.endswith(".h5") and not input_path.endswith(".h5ad"):
        adata = sc.read_10x_h5(
            input_path,
            genome=args.genome,
            gex_only=not args.include_all_assays,
        )
    else:
        mtx_dir = input_path
        if os.path.isfile(input_path) and input_path.endswith((".mtx", ".mtx.gz")):
            mtx_dir = os.path.dirname(input_path)
        if not os.path.isdir(mtx_dir) or not _is_mtx_dir(mtx_dir):
            sys.stderr.write("[error] Input is not a 10X h5 file or MTX directory.\n")
            sys.exit(2)
        adata = sc.read_10x_mtx(
            mtx_dir,
            var_names=args.var_names,
            make_unique=args.make_genes_unique,
        )

    if args.make_genes_unique and adata.var_names.is_unique is False:
        adata.var_names_make_unique()

    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    out_path = _infer_output_path(input_path, args.output)
    adata.write(out_path, compression="lzf")
    print(f"[ok] Wrote: {out_path}")
    print(f"[ok] Shape: cells x genes = {adata.n_obs} x {adata.n_vars}")


if __name__ == "__main__":
    main()
