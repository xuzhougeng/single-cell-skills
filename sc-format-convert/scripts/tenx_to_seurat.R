#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

usage <- function() {
  cat("Usage: tenx_to_seurat.R <input_10x> -o output.rds [--assay RNA] [--project Project] [--min-cells 3] [--min-features 200]\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  usage()
  quit(status = 1)
}

input <- NULL
output <- NULL
assay <- "RNA"
project <- "Seurat"
min.cells <- 3
min.features <- 200

i <- 1
while (i <= length(args)) {
  arg <- args[i]
  if (!startsWith(arg, "-") && is.null(input)) {
    input <- arg
  } else if (arg %in% c("-o", "--output")) {
    output <- args[i + 1]
    i <- i + 1
  } else if (arg == "--assay") {
    assay <- args[i + 1]
    i <- i + 1
  } else if (arg == "--project") {
    project <- args[i + 1]
    i <- i + 1
  } else if (arg == "--min-cells") {
    min.cells <- as.integer(args[i + 1])
    i <- i + 1
  } else if (arg == "--min-features") {
    min.features <- as.integer(args[i + 1])
    i <- i + 1
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
  if (dir.exists(path)) {
    base <- basename(normalizePath(path))
    return(paste0(base, ".rds"))
  }
  base <- sub("\\.h5$", "", path)
  return(paste0(base, ".rds"))
}

if (is.null(output)) {
  output <- infer_output(input)
}

data <- NULL
if (file.exists(input) && grepl("\\.h5$", input)) {
  data <- Read10X_h5(input)
} else {
  data <- Read10X(input)
}

if (is.list(data)) {
  if ("Gene Expression" %in% names(data)) {
    data <- data[["Gene Expression"]]
  } else if ("RNA" %in% names(data)) {
    data <- data[["RNA"]]
  } else {
    data <- data[[1]]
  }
}

seu <- CreateSeuratObject(
  counts = data,
  assay = assay,
  project = project,
  min.cells = min.cells,
  min.features = min.features
)

saveRDS(seu, output)
cat("[ok] Wrote:", output, "\n")
