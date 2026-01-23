#!/usr/bin/env bash
# Update single-cell skills from GitHub repository
# Source: https://github.com/xuzhougeng/single-cell-skills.git

set -euo pipefail

REPO_URL="https://github.com/xuzhougeng/single-cell-skills.git"
SKILLS_DIR=".claude/skills"

# Create temporary directory with cleanup trap
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

echo "Cloning single-cell-skills repository..."
git clone --depth 1 "$REPO_URL" "$TMP_DIR"

echo "Creating skills directory: $SKILLS_DIR"
mkdir -p "$SKILLS_DIR"

echo "Copying skills..."
copied=0
for d in "$TMP_DIR"/*/; do
  d="${d%/}"
  [ -d "$d" ] || continue
  skill_name="$(basename "$d")"
  cp -a "$d" "$SKILLS_DIR"/
  echo "  - $skill_name"
  copied=$((copied + 1))
done

echo "Done! Updated $copied skills in $SKILLS_DIR"
