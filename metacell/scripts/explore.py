#!/usr/bin/env python3

"""Explore single-cell data characteristics for all species."""

import argparse
import scanpy as sc
import pandas as pd
from pathlib import Path
import sys
import re

def explore_h5ad(filepath):
    """Explore basic characteristics of a single h5ad file."""
    print(f"\n{'='*60}")
    print(f"Processing: {filepath.name}")
    print('='*60)

    try:
        adata = sc.read_h5ad(filepath)

        species_code = filepath.stem.replace('_with_metacells', '')
        n_cells = adata.n_obs
        n_genes = adata.n_vars

        # Check for celltype column
        celltype_cols = [col for col in adata.obs.columns if 'type' in col.lower() or 'anno' in col.lower()]
        celltype_key = celltype_cols[0] if celltype_cols else None

        # Check gene prefixes for organellar genes
        gene_names = adata.var_names.tolist()

        # Common (imperfect) heuristics for organellar gene naming across species.
        # We keep this lightweight and report detected "prefix families" for quick inspection.
        mito_prefixes = set()
        chloro_prefixes = set()
        leading_alpha = re.compile(r"^([A-Za-z]+(?:-[A-Za-z]+)?)")

        def _prefix_family(gene: str) -> str:
            m = leading_alpha.match(gene)
            if m:
                return m.group(1)
            # Fallback: first 5 chars is usually enough to spot conventions.
            return gene[:5]

        for gene in gene_names:
            g = str(gene)
            g_up = g.upper()
            # Mitochondria
            if g_up.startswith(("ATMG", "ATM", "MT-", "MT_", "MT.")):
                mito_prefixes.add(_prefix_family(g))
            # Chloroplast / plastid
            if g_up.startswith(("ATCG", "ATC", "PT-", "PT_", "PT.", "CHLORO")):
                chloro_prefixes.add(_prefix_family(g))

        # Determine target metacell size based on cell count
        if n_cells < 5000:
            target_size = 18
        elif n_cells < 10000:
            target_size = 36
        elif n_cells < 20000:
            target_size = 72
        else:
            target_size = 108

        result = {
            'species': species_code,
            'n_cells': n_cells,
            'n_genes': n_genes,
            'celltype_key': celltype_key,
            'target_size': target_size,
            'mito_prefix': ','.join(sorted(mito_prefixes)) if mito_prefixes else 'None',
            'chloro_prefix': ','.join(sorted(chloro_prefixes)) if chloro_prefixes else 'None'
        }

        print(f"  Cells: {n_cells:,}")
        print(f"  Genes: {n_genes:,}")
        print(f"  Celltype column: {celltype_key}")
        print(f"  Target metacell size: {target_size}")
        print(f"  Mitochondrial prefixes: {result['mito_prefix']}")
        print(f"  Chloroplast prefixes: {result['chloro_prefix']}")

        if celltype_key:
            n_types = adata.obs[celltype_key].nunique()
            print(f"  Number of cell types: {n_types}")
            print(f"  Cell types: {', '.join(adata.obs[celltype_key].unique()[:10])}")
            if n_types > 10:
                print(f"    (showing first 10 of {n_types})")

        return result

    except Exception as e:
        print(f"  ERROR: {e}")
        return None

def _parse_args():
    parser = argparse.ArgumentParser(
        description="Explore basic characteristics of one or more .h5ad files."
    )
    parser.add_argument(
        "path",
        nargs="?",
        default="h5ad",
        help="Path to a .h5ad file or a directory containing .h5ad files (default: h5ad).",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help=(
            "Where to write the summary CSV. Use '-' to write to stdout. "
            "Default: metacell/species_characteristics.csv next to this skill."
        ),
    )
    return parser.parse_args()


def _collect_h5ad_files(path: Path):
    if path.is_file():
        if path.suffix != ".h5ad":
            print(f"Error: {path} is not a .h5ad file")
            sys.exit(1)
        return [path]

    if not path.exists():
        print(f"Error: {path} directory not found")
        sys.exit(1)

    return sorted([
        f for f in path.glob("*.h5ad")
        if not f.name.startswith(('_', '.')) and
        'with_metacells' not in f.name and
        'metacells_by' not in f.name and
        not f.name.endswith('.csv')
    ])

def _default_output_path() -> Path:
    # Place next to the metacell skill folder, regardless of current working directory.
    # metacell/scripts/explore.py -> metacell/species_characteristics.csv
    return (Path(__file__).resolve().parents[1] / "species_characteristics.csv")


def main():
    args = _parse_args()
    path = Path(args.path)

    h5ad_files = _collect_h5ad_files(path)

    print(f"Found {len(h5ad_files)} h5ad files to process")

    results = []
    for filepath in h5ad_files:
        result = explore_h5ad(filepath)
        if result:
            results.append(result)

    # Save summary
    if results:
        df = pd.DataFrame(results)
        if args.output == "-":
            df.to_csv(sys.stdout, index=False)
            output_file = "<stdout>"
        else:
            output_file = Path(args.output).resolve() if args.output else _default_output_path()
            output_file.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(output_file, index=False)
        print(f"\n{'='*60}")
        print(f"Summary saved to: {output_file}")
        print('='*60)
        print(df.to_string(index=False))

if __name__ == "__main__":
    main()
