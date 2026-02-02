---
name: sc-llm-annotate
description: Annotate single-cell datasets with LLM-guided cell types using user-provided marker genes. Use when working with h5ad or Seurat (RDS) objects and the user requests dotplot-based resolution selection and cell type labeling.
---

# sc-llm-annotate

## Purpose

Guide the workflow to use LLMs and marker genes to annotate cell types in `h5ad` or Seurat (`.rds/.qs/.Rdata`) objects by comparing dotplots across clustering resolutions, selecting a suitable resolution, and adding `cell_type_llm` labels.

## Trigger Scenarios

Use this skill when the user asks to:
- Provide marker genes to label cell types.
- Visualize markers across multiple resolutions via dotplots.
- Choose a resolution for annotation in `h5ad` or Seurat objects.
- Add a new cell type column such as `cell_type_llm`.

## Required Inputs

- A marker gene list, grouped by cell type.
- A dataset path: `h5ad` or Seurat `.rds`.
- Candidate clustering resolutions (if not provided, infer from object metadata).

### Marker format (example)

```
T cell: CD3D CD3E LST1
B cell: MS4A1 CD79A CD79B
Monocyte: LST1 S100A8 S100A9
```

## Workflow

1. **Load dataset**
   - `h5ad`: use `scanpy` / `anndata`.
   - `Seurat`: use Seurat in R.

2. **Identify resolution columns**
   - Find metadata columns that match typical resolution naming (e.g., `leiden_0.5`, `seurat_clusters`, `snn_res.0.8`).
   - If none exist, prompt to run clustering first.

3. **Generate and save two outputs per resolution**
   - **Dotplot images**: For each candidate resolution, plot marker expression by cluster. Keep plots consistent (same marker order and color scale). Save as `dotplot_res0.5.png`, `dotplot_res0.8.png`, etc.
   - **Numeric summary**: Compute per-cluster marker expression matrix (mean/median % expressed and expression level). Save as CSV or TSV (e.g., `marker_matrix_res0.5.csv`).

4. **LLM selects resolution**
   - **If model has vision**: Present the dotplot images to the LLM; it visually inspects and selects the best resolution.
   - **If model lacks vision**: Present the numeric expression matrices; LLM reasons from cluster×marker values.
   - Criteria: distinct marker patterns, minimal mixed signals, balance interpretability and granularity.

5. **LLM assigns cell types**
   - **If model has vision**: Present the dotplot(s) of the selected resolution; LLM reads images and maps clusters to cell types.
   - **If model lacks vision**: Present the numeric marker matrix; LLM reasons from expression values to map clusters to cell types.
   - Produce a mapping: `cluster_id -> cell_type`.

6. **Write annotations**
   - Create a new column `cell_type_llm` in metadata.
   - Preserve any existing labels; do not overwrite unless asked.

7. **Save outputs**
   - Save updated `h5ad` or Seurat object.
   - Export a CSV of cluster-to-cell-type mappings.
   - Keep both dotplot images and numeric marker matrices for reference.

## Output Template

```
Dotplot images: <output_dir>/dotplot_res*.png
Numeric matrices: <output_dir>/marker_matrix_res*.csv
Selected resolution: <resolution_column>
Cluster mapping:
  <cluster_id> -> <cell_type>
  <cluster_id> -> <cell_type>

Added column: cell_type_llm
Saved: <output_path>
```

## Notes

- **Always output both**: dotplot images and per-cluster marker expression matrices (numeric summary). The LLM uses whichever input type it supports.
- **Vision vs numeric fallback**: If the LLM has vision, present dotplot images. If not, use the numeric matrices (cluster × marker expression) for resolution selection and cell type mapping.
- If multiple resolutions look plausible, document why the chosen one is best.
- When markers are ambiguous, label as `Unknown` and flag for follow-up.