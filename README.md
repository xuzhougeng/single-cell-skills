# 单细胞数据分析相关的 Skills

## 安装方法

建议安装到项目的目录下，而不是全局，避免技能太多

对于claude

```bash
SKILLS_DIR=".claude/skills"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

git clone https://github.com/xuzhougeng/single-cell-skills.git "$TMP_DIR"
mkdir -p "$SKILLS_DIR"

# 将仓库里的各个 skill 目录合并到目标目录（不会删除目标目录已有内容；同名文件会覆盖）
for d in "$TMP_DIR"/*/; do
  d="${d%/}"
  [ -d "$d" ] || continue
  cp -a "$d" "$SKILLS_DIR"/
done
```

对于codex

```bash
SKILLS_DIR=".codex/skills"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

git clone https://github.com/xuzhougeng/single-cell-skills.git "$TMP_DIR"
mkdir -p "$SKILLS_DIR"

# 将仓库里的各个 skill 目录合并到目标目录（不会删除目标目录已有内容；同名文件会覆盖）
for d in "$TMP_DIR"/*/; do
  d="${d%/}"
  [ -d "$d" ] || continue
  cp -a "$d" "$SKILLS_DIR"/
done
```

注：仓库自带一个SKILL, 用于更新技能，因此只需要在cursor/cladue/codex中发送, `update sc skills` 即可更新。

## 使用方法

claude

```bash
将 demoh5ad 转成 Seuat, demo.rds
```

此时会有一个确认，是否需要用skill

```text
 Use skill "sc-format-convert"?                                                                                                                    
 Claude may use instructions, code, or files from this Skill.                                                                                      
                                                                                                                                                   
   Convert single-cell data formats between Seurat (.rds/.qs), AnnData (.h5ad), and 10X (h5/MTX). Use when users ask to convert Seurat <-> h5ad,   
    10X -> Seurat, or 10X -> h5ad. 
```

此外，可以通过用 `/技能名`的方式强制调用

```bash
/sc-format-convert 将 demoh5ad 转成 Seuat, demo.rds
```


## SKILL列表

- sc-info: 了解单细胞数据一些基本信息
- sc-format-convert: 单细胞格式转换, 例如seurat的rds/qs到scanpy的h5ad
- seurat-slim: 对Seurat对象进行瘦身，核心是去掉 scaled.data
