#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SeuratDisk)
})

usage <- function() {
  cat("Usage: h5ad_to_seurat.R <input.h5ad> -o output.rds [--keep-intermediate]\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  usage()
  quit(status = 1)
}

input <- NULL
output <- NULL
keep_intermediate <- FALSE

i <- 1
while (i <= length(args)) {
  arg <- args[i]
  if (!startsWith(arg, "-") && is.null(input)) {
    input <- arg
  } else if (arg %in% c("-o", "--output")) {
    output <- args[i + 1]
    i <- i + 1
  } else if (arg == "--keep-intermediate") {
    keep_intermediate <- TRUE
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

h5seurat <- sub("\\.h5ad$", ".h5seurat", input)

Convert(input, dest = "h5seurat", overwrite = TRUE)
seu <- LoadH5Seurat(h5seurat)
saveRDS(seu, output)
cat("[ok] Wrote:", output, "\n")

if (!keep_intermediate && file.exists(h5seurat)) {
  file.remove(h5seurat)
}
