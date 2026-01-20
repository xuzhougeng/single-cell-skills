#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: seurat_slim.R <input_file> <output_file>\n")
  quit(status = 1)
}

input_file <- args[[1]]
output_file <- args[[2]]
input_ext <- tolower(tools::file_ext(input_file))
output_ext <- tolower(tools::file_ext(output_file))

if (!(input_ext %in% c("rds", "qs"))) {
  cat("Unsupported input file type:", input_ext, "\n")
  quit(status = 1)
}
if (!(output_ext %in% c("rds", "qs"))) {
  cat("Unsupported output file type:", output_ext, "\n")
  quit(status = 1)
}
if (input_ext == "qs" || output_ext == "qs") {
  if (!requireNamespace("qs", quietly = TRUE)) {
    cat("Package 'qs' is required to handle .qs files\n")
    quit(status = 1)
  }
}

output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

cat("========================================\n")
cat("Processing:", basename(input_file), "\n")
cat("========================================\n")

if (file.exists(output_file)) {
  cat("  Output exists, skipping\n")
  quit(status = 0)
}

cat("  Reading file...\n")
if (input_ext == "qs") {
  obj <- qs::qread(input_file)
} else {
  obj <- readRDS(input_file)
}

original_size <- file.info(input_file)$size / 1024^3
cat("  Original size:", round(original_size, 2), "GB\n")

for (assay_name in Assays(obj)) {
  assay <- obj@assays[[assay_name]]

  if (.hasSlot(assay, "layers")) {
    layers <- Layers(assay)
    scale_layers <- layers[grepl("^scale", layers)]

    if (length(scale_layers) > 0) {
      cat("  Assay [", assay_name, "] (v5) removing layers:", paste(scale_layers, collapse = ", "), "\n")
      for (layer in scale_layers) {
        obj@assays[[assay_name]]@layers[[layer]] <- NULL
      }
    }
  } else if (.hasSlot(assay, "scale.data")) {
    if (length(assay@scale.data) > 0 && nrow(assay@scale.data) > 0) {
      cat("  Assay [", assay_name, "] (v4) removing scale.data\n")
      obj@assays[[assay_name]]@scale.data <- matrix(nrow = 0, ncol = 0)
    }
  }
}

cat("  Saving to:", output_file, "\n")
if (output_ext == "qs") {
  qs::qsave(obj, output_file)
} else {
  saveRDS(obj, output_file)
}

new_size <- file.info(output_file)$size / 1024^3
reduction <- (1 - new_size / original_size) * 100

cat("  New size:", round(new_size, 2), "GB\n")
cat("  Reduced:", round(reduction, 1), "%\n")

rm(obj)
gc(verbose = FALSE)

cat("========================================\n")
cat("All done\n")
cat("Output file:", output_file, "\n")
