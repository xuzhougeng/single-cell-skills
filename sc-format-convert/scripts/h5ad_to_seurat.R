#!/usr/bin/env Rscript

usage <- function() {
  cat(paste0(
    "Usage: h5ad_to_seurat.R <input.h5ad> -o output.rds [--use-raw] [--spatial] [--no-simplify]\n",
    "\n",
    "Notes:\n",
    "  - This script uses schard (reticulate-free) to load .h5ad into a Seurat object.\n",
    "  - Install schard: devtools::install_github(\"cellgeni/schard\")\n",
    "  - --use-raw maps to schard::h5ad2seurat(..., use.raw=TRUE)\n",
    "  - --spatial uses schard::h5ad2seurat_spatial() (for Visium/spatial h5ad)\n"
  ))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  usage()
  quit(status = 1)
}

input <- NULL
output <- NULL
use_raw <- FALSE
spatial <- FALSE
simplify <- TRUE

i <- 1
while (i <= length(args)) {
  arg <- args[i]
  if (!startsWith(arg, "-") && is.null(input)) {
    input <- arg
  } else if (arg %in% c("-o", "--output")) {
    output <- args[i + 1]
    i <- i + 1
  } else if (arg == "--use-raw") {
    use_raw <- TRUE
  } else if (arg == "--spatial") {
    spatial <- TRUE
  } else if (arg == "--no-simplify") {
    simplify <- FALSE
  } else if (arg == "--keep-intermediate") {
    # Backward-compatible no-op (old SeuratDisk workflow created .h5seurat)
    message("[info] --keep-intermediate is ignored when using schard.")
  } else if (arg == "--help") {
    usage()
    quit(status = 0)
  }
  i <- i + 1
}

if (is.null(input)) {
  usage()
  quit(status = 1)
}

if (!file.exists(input)) {
  stop(paste("Input not found:", input))
}

infer_output <- function(path) {
  base <- sub("\\.h5ad$", "", path)
  return(paste0(base, ".rds"))
}

if (is.null(output)) {
  output <- infer_output(input)
}

if (!requireNamespace("schard", quietly = TRUE)) {
  stop(paste0(
    "R package 'schard' is required.\n",
    "Install it with:\n",
    "  install.packages('devtools')\n",
    "  devtools::install_github('cellgeni/schard')\n"
  ))
}

# Seurat is usually needed to create/use Seurat objects
if (!requireNamespace("Seurat", quietly = TRUE) && !requireNamespace("SeuratObject", quietly = TRUE)) {
  message("[warn] 'Seurat' (or 'SeuratObject') not detected. schard may still work, but Seurat-dependent operations could fail.")
}

seu <- if (isTRUE(spatial)) {
  schard::h5ad2seurat_spatial(input, use.raw = use_raw, simplify = simplify)
} else {
  schard::h5ad2seurat(input, use.raw = use_raw)
}

saveRDS(seu, output)
cat("[ok] Wrote:", output, "\n")
