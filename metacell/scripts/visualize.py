#!/usr/bin/env python3

"""Simple clean metacell visualization - scatter plot style like the paper."""

import argparse
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from pathlib import Path
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgba
import scipy.sparse as sp
from scipy.spatial import cKDTree


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("cells", help="Single-cell h5ad")
    parser.add_argument("metacells", help="Metacell h5ad")
    parser.add_argument("-o", "--output", default="metacell_clean.pdf")
    parser.add_argument("--title", default="Metacell Visualization")
    parser.add_argument("--celltype-key", default="celltype_anno")
    parser.add_argument(
        "--colors",
        default="",
        help=(
            "Optional explicit color mapping for cell types. "
            "Accepts a string like \"'SCN/Me'='#E9967A', 'Ph SE'='#56b417'\". "
            "Keys must match values in --celltype-key. "
            "Unspecified types fall back to the default palette."
        ),
    )
    parser.add_argument("--figsize", type=float, nargs=2, default=[10, 8])
    parser.add_argument("--node-size", type=float, default=100)
    parser.add_argument(
        "--node-alpha",
        type=float,
        default=1.0,
        help="Alpha for metacell (top-layer) points. Range: 0~1 (default: 1.0).",
    )
    parser.add_argument(
        "--node-linewidth",
        type=float,
        default=1.6,
        help="Outline width for metacell points (default: 1.6).",
    )
    parser.add_argument(
        "--outline-darken",
        type=float,
        default=0.7,
        help="Darken factor for outline color (default: 0.7).",
    )
    parser.add_argument(
        "--draw-nn",
        action="store_true",
        help="If set, draw nearest-neighbor edges between metacell points (based on metacell 2D coords).",
    )
    parser.add_argument(
        "--nn-k",
        type=int,
        default=1,
        help="Number of nearest neighbors per metacell to connect (default: 1).",
    )
    parser.add_argument(
        "--nn-alpha",
        type=float,
        default=0.15,
        help="Alpha for NN edges (default: 0.15).",
    )
    parser.add_argument(
        "--nn-linewidth",
        type=float,
        default=0.6,
        help="Line width for NN edges (default: 0.6).",
    )
    return parser.parse_args()

def parse_color_mapping(s: str) -> dict[str, str]:
    """
    Parse a color mapping string like:
      "'SCN/Me'='#E9967A', 'Ph SE'='#56b417'"
    into {"SCN/Me": "#E9967A", "Ph SE": "#56b417"}.
    """
    if not s or not s.strip():
        return {}

    # Accept common wrappers like c(...) from R if user passes it.
    txt = s.strip()
    if txt.startswith("c(") and txt.endswith(")"):
        txt = txt[2:-1].strip()

    # Match 'key'='#hex' or "key"="#hex" with flexible spacing.
    pairs = re.findall(
        r"""['"]([^'"]+)['"]\s*=\s*['"](#(?:[0-9A-Fa-f]{6}|[0-9A-Fa-f]{8}))['"]""",
        txt,
    )
    if not pairs:
        raise SystemExit(
            "[error] Failed to parse --colors. "
            "Expected pairs like \"'Type'='#RRGGBB'\" separated by commas."
        )
    return {k: v for k, v in pairs}


def _is_likely_counts(adata) -> bool:
    # If scanpy log1p metadata exists, treat as already logged.
    if isinstance(getattr(adata, "uns", None), dict) and adata.uns.get("log1p") is not None:
        return False
    X = adata.X
    rng = np.random.default_rng(0)

    if sp.issparse(X):
        data = X.data
        if data.size == 0:
            return False
        sample = data if data.size <= 10_000 else rng.choice(data, size=10_000, replace=False)
    else:
        arr = np.asarray(X)
        flat = arr.ravel()
        if flat.size == 0:
            return False
        sample = flat if flat.size <= 10_000 else rng.choice(flat, size=10_000, replace=False)

    # Heuristic: counts are non-negative and often have large-ish maxima.
    try:
        if np.issubdtype(sample.dtype, np.integer):
            return True
    except Exception:
        pass
    if np.nanmin(sample) < 0:
        return False
    return np.nanmax(sample) > 50 and np.nanmedian(sample) >= 1

def compute_umap_if_needed(adata):
    """Compute UMAP if not present."""
    if 'X_umap' in adata.obsm:
        print("[info] Using existing UMAP")
        return adata

    print("[info] Computing UMAP")
    if _is_likely_counts(adata):
        print("[info] Input looks like counts; applying normalize_total + log1p for embedding")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
    sc.pp.pca(adata, n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata, random_state=42)

    return adata


def _darken_rgba(rgba, factor: float = 0.7):
    """Darken an RGBA color by scaling RGB channels."""
    r, g, b, a = rgba
    return (r * factor, g * factor, b * factor, a)


def draw_nearest_neighbor_edges(
    ax,
    pos: np.ndarray,
    *,
    k: int = 1,
    color: str = "k",
    alpha: float = 0.15,
    linewidth: float = 0.6,
    zorder: int = 2,
):
    """
    Draw k-NN edges for 2D coordinates.
    pos: (n, 2)
    """
    if pos.shape[0] < 2:
        return
    if k < 1:
        return

    k_eff = min(k, pos.shape[0] - 1)
    tree = cKDTree(pos)
    # include self => query k_eff + 1
    _, idxs = tree.query(pos, k=k_eff + 1)
    nbrs = np.atleast_2d(idxs)[:, 1:]  # drop self

    i = np.arange(pos.shape[0])[:, None]
    a = np.minimum(i, nbrs).ravel()
    b = np.maximum(i, nbrs).ravel()
    edges = np.unique(np.stack([a, b], axis=1), axis=0)

    segments = np.stack([pos[edges[:, 0]], pos[edges[:, 1]]], axis=1)
    lc = LineCollection(segments, colors=color, linewidths=linewidth, alpha=alpha, zorder=zorder)
    ax.add_collection(lc)


def project_metacells(adata_cells, adata_meta):
    """Project metacells to UMAP."""
    print("[info] Projecting metacells to UMAP")

    cell_umap = adata_cells.obsm['X_umap']
    if 'metacell_name' not in adata_cells.obs.columns:
        raise SystemExit(
            "[error] cells AnnData is missing obs['metacell_name'] required for projection. "
            "Run metacell/scripts/pipeline.py to generate <input>_with_metacells.h5ad."
        )
    cell_mc = adata_cells.obs['metacell_name'].astype(str).values

    mc_names = adata_meta.obs_names
    mc_umap = np.full((len(mc_names), 2), np.nan, dtype=float)

    for i, mc_name in enumerate(mc_names):
        mask = cell_mc == mc_name
        if mask.sum() > 0:
            mc_umap[i] = cell_umap[mask].mean(axis=0)

    adata_meta.obsm['X_umap'] = mc_umap
    return adata_meta


def plot_clean(
    adata_cells,
    adata_meta,
    output,
    title,
    celltype_key,
    color_mapping: dict[str, str],
    figsize,
    node_size,
    node_alpha: float,
    node_linewidth: float,
    outline_darken: float,
    draw_nn: bool,
    nn_k: int,
    nn_alpha: float,
    nn_linewidth: float,
):
    """Clean scatter plot visualization."""
    print("[info] Creating clean visualization")

    cell_umap = adata_cells.obsm['X_umap']
    mc_umap = adata_meta.obsm['X_umap']
    mc_finite = np.isfinite(mc_umap).all(axis=1)

    # Get cell types
    if celltype_key in adata_meta.obs.columns:
        mc_types = adata_meta.obs[celltype_key].astype(str).values
    else:
        mc_types = np.array(['Unknown'] * adata_meta.n_obs)

    # Colors
    unique_types = sorted(set(mc_types))
    n_types = len(unique_types)
    type_colors: dict[str, object] = {}
    if color_mapping:
        # Use explicit mapping when available; fill the rest with a palette.
        missing = [t for t in unique_types if t not in color_mapping]
        palette = sns.color_palette('tab20' if len(missing) > 10 else 'tab10', len(missing))
        for t in unique_types:
            if t in color_mapping:
                type_colors[t] = color_mapping[t]
        for t, c in zip(missing, palette):
            type_colors[t] = c
    else:
        palette = sns.color_palette('tab20' if n_types > 10 else 'tab10', n_types)
        type_colors = dict(zip(unique_types, palette))

    # Plot
    fig, ax = plt.subplots(figsize=figsize)

    # Background cells (gray, larger)
    ax.scatter(cell_umap[:, 0], cell_umap[:, 1],
              s=3, c='lightgray', alpha=0.3,
              rasterized=True, zorder=1)

    # Optional: nearest-neighbor edges between metacells (on the current 2D coords)
    if draw_nn:
        pos = mc_umap[mc_finite]
        draw_nearest_neighbor_edges(
            ax,
            pos,
            k=nn_k,
            alpha=nn_alpha,
            linewidth=nn_linewidth,
            zorder=2,
        )

    # Metacells as scatter points
    for ct in unique_types:
        mask = (mc_types == ct) & mc_finite
        if mask.sum() > 0:
            face = to_rgba(type_colors[ct], alpha=node_alpha)
            edge = _darken_rgba(face, factor=outline_darken)
            ax.scatter(mc_umap[mask, 0], mc_umap[mask, 1],
                      s=node_size, c=[face],
                      edgecolors=[edge], linewidths=node_linewidth,
                      label=ct, zorder=3)

    # Legend
    ax.legend(bbox_to_anchor=(1.02, 0.5), loc='center left',
             frameon=True, fontsize=9)

    # Labels
    ax.set_xlabel('UMAP 1', fontsize=11)
    ax.set_ylabel('UMAP 2', fontsize=11)

    n_cells = len(adata_cells)
    n_mc = adata_meta.n_obs
    ax.set_title(f"{title}\n{n_cells:,} cells, {n_mc} metacells",
                fontsize=12, fontweight='bold')

    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"[ok] Saved to {output}")


def main():
    args = parse_args()
    color_mapping = parse_color_mapping(args.colors)

    print(f"[info] Loading cells from {args.cells}")
    adata_cells = sc.read_h5ad(args.cells)

    print(f"[info] Loading metacells from {args.metacells}")
    adata_meta = sc.read_h5ad(args.metacells)

    # Compute UMAP
    adata_cells = compute_umap_if_needed(adata_cells)

    # Project metacells
    adata_meta = project_metacells(adata_cells, adata_meta)

    # Plot
    output_path = Path(args.output)
    if str(output_path.parent) not in (".", ""):
        output_path.parent.mkdir(parents=True, exist_ok=True)
    plot_clean(
        adata_cells,
        adata_meta,
        str(output_path),
        args.title,
        args.celltype_key,
        color_mapping,
        tuple(args.figsize),
        args.node_size,
        args.node_alpha,
        args.node_linewidth,
        args.outline_darken,
        args.draw_nn,
        args.nn_k,
        args.nn_alpha,
        args.nn_linewidth,
    )


if __name__ == "__main__":
    main()
