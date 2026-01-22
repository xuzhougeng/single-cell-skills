#!/usr/bin/env python3

"""Run the Metacell divide-and-conquer pipeline on an AnnData object."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import scanpy as sc

try:
    import metacells as mc
except ImportError as exc:
    raise SystemExit(
        "The 'metacells' package is required. Install it with `pip install metacells`."
    ) from exc


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Construct metacells from a single-cell AnnData object using the "
            "divide-and-conquer pipeline."
        )
    )
    parser.add_argument(
        "adata",
        help="Input AnnData (.h5ad) file containing single-cell data.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output path for the metacell AnnData (defaults to <input>.metacells.h5ad).",
    )
    parser.add_argument(
        "--target-metacell-size",
        type=int,
        default=96,
        help=(
            "Approximate number of cells per metacell. "
            "Use a smaller value (e.g. 72) when per-cell quality is high. "
            "Default: 96."
        ),
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=123456,
        help="Random seed forwarded to the metacell pipeline. Default: 123456.",
    )
    parser.add_argument(
        "--metacell-name",
        default="preliminary.metacells",
        help=(
            "Bookkeeping name stored in AnnData.uns and used by metacells when "
            "collecting metacells. Default: preliminary.metacells."
        ),
    )
    parser.add_argument(
        "--group-key",
        default="group",
        help=(
            "Observation column to propagate to metacells and compute fractions for. "
            "Set to empty string to skip this step. Default: group."
        ),
    )
    parser.add_argument(
        "--max-parallel-piles",
        type=int,
        help="Override the guessed number of parallel piles.",
    )
    parser.add_argument(
        "--lateral-gene",
        dest="lateral_genes",
        action="append",
        default=[],
        help="Explicit gene marked as lateral (repeat flag for multiple genes).",
    )
    parser.add_argument(
        "--lateral-gene-pattern",
        dest="lateral_patterns",
        action="append",
        default=[],
        help="Regular expression marking lateral genes (repeat as needed).",
    )
    parser.add_argument(
        "--noisy-gene",
        dest="noisy_genes",
        action="append",
        default=[],
        help="Explicit gene marked as noisy (repeat flag for multiple genes).",
    )
    parser.add_argument(
        "--output-cells",
        help=(
            "Optional output path to write the original single-cell AnnData with the "
            "cell-to-metacell assignments (obs['metacell_name']). Defaults to "
            "<input>_with_metacells.h5ad next to the input."
        ),
    )
    parser.add_argument(
        "--no-output-cells",
        action="store_true",
        help="Do not write the single-cell AnnData with metacell assignments.",
    )
    return parser.parse_args()


def _resolve_output_path(input_path: Path, explicit_output: str | None) -> Path:
    if explicit_output:
        return Path(explicit_output).resolve()

    if input_path.name.endswith(".h5ad"):
        output_name = input_path.name[:-5] + ".metacells.h5ad"
    else:
        output_name = input_path.name + ".metacells.h5ad"
    return (input_path.parent / output_name).resolve()


def _resolve_output_cells_path(input_path: Path, explicit_output: str | None) -> Path:
    if explicit_output:
        return Path(explicit_output).resolve()
    return (input_path.parent / f"{input_path.stem}_with_metacells.h5ad").resolve()


def _log(message: str) -> None:
    sys.stderr.write(message + "\n")


def main() -> None:
    args = _parse_args()

    adata_path = Path(args.adata).resolve()
    if not adata_path.exists():
        raise SystemExit(f"Input AnnData not found: {adata_path}")

    output_path = _resolve_output_path(adata_path, args.output)

    _log(f"[info] Reading cells from {adata_path}")
    adata = sc.read_h5ad(adata_path)

    # Ensure data is float type for metacells pipeline
    import numpy as np
    if adata.X.dtype not in (np.float32, np.float64):
        _log(f"[info] Converting data from {adata.X.dtype} to float32")
        adata.X = adata.X.astype(np.float32)

    mc.pl.mark_lateral_genes(
        adata,
        lateral_gene_names=args.lateral_genes,
        lateral_gene_patterns=args.lateral_patterns,
    )
    if args.noisy_genes:
        mc.pl.mark_noisy_genes(adata, noisy_gene_names=args.noisy_genes)
    else:
        mc.pl.mark_noisy_genes(adata, noisy_gene_names=[])

    if args.max_parallel_piles:
        max_parallel = args.max_parallel_piles
    else:
        max_parallel = mc.pl.guess_max_parallel_piles(adata)
    mc.pl.set_max_parallel_piles(max_parallel)
    _log(f"[info] Using up to {max_parallel} parallel piles")

    _log(
        "[info] Running divide-and-conquer pipeline with target size "
        f"{args.target_metacell_size}"
    )
    with mc.ut.progress_bar():
        mc.pl.divide_and_conquer_pipeline(
            adata,
            target_metacell_size=args.target_metacell_size,
            random_seed=args.random_seed,
        )

    metacells = mc.pl.collect_metacells(
        adata,
        name=args.metacell_name,
        random_seed=args.random_seed,
    )
    _log(
        f"[info] Built {metacells.n_obs} metacells covering {metacells.n_vars} genes"
    )

    group_key = (args.group_key or "").strip()
    if group_key:
        if group_key not in adata.obs.columns:
            _log(f"[warn] '{group_key}' not found in adata.obs; skipping group transfer")
        else:
            mc.tl.convey_obs_to_group(
                adata=adata,
                gdata=metacells,
                property_name=group_key,
                to_property_name=group_key,
                method=mc.ut.most_frequent,
            )
            mc.tl.convey_obs_fractions_to_group(
                adata=adata,
                gdata=metacells,
                property_name=group_key,
                to_property_name=group_key,
            )
            _log(f"[info] Propagated '{group_key}' and its fractions to metacells")

    metacells.write(output_path, compression="lzf")
    _log(f"[ok] Wrote metacells to {output_path}")

    # Save cell-to-metacell assignments back to the single-cell AnnData.
    # collect_metacells() should populate obs['metacell_name'] (unassigned => 'Outliers').
    if "metacell_name" not in adata.obs.columns:
        raise SystemExit(
            "[error] Missing adata.obs['metacell_name'] after metacell collection. "
            "This indicates the pipeline did not produce cell-to-metacell assignments."
        )

    if not args.no_output_cells:
        output_cells = _resolve_output_cells_path(adata_path, args.output_cells)
        adata.write(output_cells, compression="lzf")
        _log(f"[ok] Wrote cells with metacell assignments to {output_cells}")


if __name__ == "__main__":
    main()
