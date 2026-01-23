---
name: update-sc-skills
description: Update or install single-cell analysis skills from the xuzhougeng/single-cell-skills GitHub repository. Use when the user asks to update, install, sync, or refresh single-cell skills, or mentions updating skills from GitHub.
---

# Update Single-Cell Skills

Syncs single-cell analysis skills from [xuzhougeng/single-cell-skills](https://github.com/xuzhougeng/single-cell-skills.git) to the project's `.claude/skills/` directory.

## Available Skills in Repository

- **sc-info**: Get basic information about single-cell datasets
- **sc-format-convert**: Convert between Seurat (.rds/.qs) and AnnData (.h5ad) formats
- **seurat-slim**: Slim Seurat objects by removing scale.data
- **metacell**: Build metacells from single-cell RNA-seq data
- **star-solo**: STARsolo alignment workflow

## Usage

Run the update script:

```bash
bash scripts/update.sh
```

## What It Does

1. Clones the repository to a temporary directory
2. Copies all skill directories to `.claude/skills/`
3. Preserves existing skills (only overwrites same-named skills)
4. Cleans up temporary files automatically

## Manual Update

If you prefer to run commands directly:

```bash
SKILLS_DIR=".claude/skills"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

git clone https://github.com/xuzhougeng/single-cell-skills.git "$TMP_DIR"
mkdir -p "$SKILLS_DIR"

for d in "$TMP_DIR"/*/; do
  d="${d%/}"
  [ -d "$d" ] || continue
  cp -a "$d" "$SKILLS_DIR"/
done
```
