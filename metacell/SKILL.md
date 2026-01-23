---
name: metacell
description: Build metacells (metacell aggregation) from single-cell RNA-seq AnnData (.h5ad) using the metacells divide-and-conquer pipeline. Use when constructing metacell objects, choosing metacell size, propagating group annotations, or visualizing metacell embeddings from scRNA-seq data.
---

# Metacell construction

Use the bundled scripts to explore datasets, run the metacell pipeline, and visualize results. Prefer these scripts over reimplementing logic.

## Quick start

1) Explore h5ad inputs to choose parameters.
2) Run the metacell divide-and-conquer pipeline.
3) Visualize metacells on UMAP.

## Scripts

- `scripts/explore.py`
  - Scan a `.h5ad` file or a directory (default: `h5ad/`) and summarize cell/gene counts, likely cell-type column, and suggested target metacell size.
  - Writes a CSV summary (default: `metacell/species_characteristics.csv` next to this skill; configurable via `--output`).
  - Run when you need a quick survey before choosing `--target-metacell-size`.

- `scripts/pipeline.py`
  - Run the metacell divide-and-conquer pipeline on a single `.h5ad`.
  - Writes `<input>.metacells.h5ad` by default.
  - By default also writes `<input>_with_metacells.h5ad` containing `obs["metacell_name"]` (cell → metacell assignment). Control via `--output-cells` / `--no-output-cells`.
  - Optionally propagates a group column from cells to metacells (`--group-key`).
  - Use for the main metacell construction.

- `scripts/visualize.py`
  - Project metacells onto the single-cell UMAP and draw a clean scatter plot.
  - Computes UMAP if missing.
  - Uses a fast edge renderer (LineCollection) and supports filtering/capping edges.
  - Use for publication-style overview plots.

## Typical workflow

1) Explore inputs and choose target size.
   - Run: `python scripts/explore.py` (defaults to `h5ad/`)
   - Or single file: `python scripts/explore.py path/to/sample.h5ad`
   - Optional: `-o/--output path/to/summary.csv` (use `--output -` to write CSV to stdout)
   - Use the suggested `target_size` per dataset.

2) Build metacells for one dataset.
   - Run: `python scripts/pipeline.py path/to/data.h5ad --target-metacell-size 96`
   - This produces:
     - `<input>.metacells.h5ad` (metacell-level AnnData)
     - `<input>_with_metacells.h5ad` (cell-level AnnData with `obs["metacell_name"]`)
   - Add group propagation if needed: `--group-key celltype_anno` (set to empty to skip).
   - Use `--lateral-gene`, `--lateral-gene-pattern`, or `--noisy-gene` when needed.

3) Visualize metacells.
   - Run: `python scripts/visualize.py cells_with_metacells.h5ad cells.metacells.h5ad -o metacells.pdf`
   - Adjust `--celltype-key` to match annotation column in metacell obs.
   - Optional edge controls:
     - `--edge-weight-min 0.1` to drop weak edges
     - `--max-edges 50000` to cap total edges drawn (speeds up large graphs)

## Output documentation requirement

After completing the full workflow, you **must** generate a `.md` documentation file with the same basename as the input file (e.g., input `pmbc.h5ad` → output `pbmc.md`).

### Documentation template

```markdown
## Dataset information
- Input file: `<input_file_path>`
- Number of cells: <n_cells>
- Number of genes: <n_genes>
- Cell types: <n_types> types (<celltype_key>)
- Recommended metacell size: <target_size>

## Run metacell pipeline

```bash
# 1. Explore dataset characteristics
python scripts/explore.py \
    <input_path> \
    --output -

# 2. Build metacells
python scripts/pipeline.py \
    <input_path> \
    --target-metacell-size <target_size> \
    --group-key <celltype_key> \
    --output <output_metacell_path> \
    --output-cells <output_cells_path>

# 3. Visualize metacells
python scripts/visualize.py <cells_with_metacells.h5ad> <metacells.h5ad> \
  -o <output_umap.png> \
  --celltype-key <celltype_key> \
  --colors "<R_style_color_vector>"
```
```

## Notes and checks

- Ensure dependencies are available: `scanpy`, `metacells`, `pandas`, `numpy`, `seaborn`, `matplotlib`.
- The pipeline expects float input; it will cast non-float `adata.X` to float32.
- `visualize.py` requires `obs["metacell_name"]` in the cell-level AnnData to project metacells. Use the `<input>_with_metacells.h5ad` produced by `pipeline.py`.
