---
name: seurat-slim
description: Slim a Seurat object by removing scale data (Seurat v4 scale.data or Seurat v5 scale* layers) to reduce file size. Use when processing a single Seurat .rds or .qs file and you need to strip scaling data across all assays.
---

# Seurat Slim

Use this skill to remove scale data from a single Seurat object stored in a .rds or .qs file.

## Quick start

Run this in R to confirm your Seurat version. If it's not 5.x or newer, configure the environment before continuing:

```r
packageVersion("Seurat")
```

Then run the bundled script:

```bash
scripts/seurat_slim.R <input_file> <output_file>
```

Example:

```bash
scripts/seurat_slim.R /path/to/sample.rds /path/to/sample_slim.rds
```

```bash
scripts/seurat_slim.R /path/to/sample.qs /path/to/sample_slim.qs
```

## Behavior

- Detect Seurat v5 assays by the `layers` slot and remove `scale*` layers.
- Detect Seurat v4 assays by the `scale.data` slot and clear it.
- Skip when the output file already exists.
- Report original size, new size, and reduction percentage.

## Notes

- Require the `Seurat` R package.
- For `.qs` input/output, require the `qs` R package.
- Expect large memory usage for big objects; run on a machine with enough RAM.
