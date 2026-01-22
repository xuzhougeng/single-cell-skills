#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
})

usage <- function() {
  cat("Usage: sc_info.R <input_path> [--json-out path] [--md-out path] [--table-out path]\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  usage()
  quit(status = 1)
}

opts <- list(json_out = NULL, md_out = NULL, table_out = NULL)
pos <- c()
i <- 1
while (i <= length(args)) {
  a <- args[[i]]
  if (startsWith(a, "--")) {
    if (i == length(args)) {
      cat("Missing value for", a, "\n")
      usage()
      quit(status = 1)
    }
    val <- args[[i + 1]]
    if (a == "--json-out") {
      opts$json_out <- val
    } else if (a == "--md-out") {
      opts$md_out <- val
    } else if (a == "--table-out") {
      opts$table_out <- val
    } else {
      cat("Unknown option:", a, "\n")
      usage()
      quit(status = 1)
    }
    i <- i + 2
  } else {
    pos <- c(pos, a)
    i <- i + 1
  }
}

if (length(pos) != 1) {
  usage()
  quit(status = 1)
}

input_path <- pos[[1]]
if (!file.exists(input_path)) {
  cat("Input does not exist:", input_path, "\n")
  quit(status = 1)
}

detect_10x_mtx_dir <- function(path) {
  if (!dir.exists(path)) return(FALSE)
  matrix_files <- c("matrix.mtx", "matrix.mtx.gz")
  barcode_files <- c("barcodes.tsv", "barcodes.tsv.gz")
  feature_files <- c("features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz")
  has_matrix <- any(file.exists(file.path(path, matrix_files)))
  has_barcodes <- any(file.exists(file.path(path, barcode_files)))
  has_features <- any(file.exists(file.path(path, feature_files)))
  has_matrix && has_barcodes && has_features
}

summarize_numeric <- function(x) {
  if (length(x) == 0) return(NULL)
  list(
    min = as.numeric(min(x)),
    median = as.numeric(stats::median(x)),
    mean = as.numeric(mean(x)),
    max = as.numeric(max(x))
  )
}

is_sparse <- function(x) {
  inherits(x, "dgCMatrix") || inherits(x, "dgTMatrix") || inherits(x, "dgRMatrix")
}

load_counts <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "h5") {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      cat("Package 'Seurat' is required for 10X .h5\n")
      quit(status = 1)
    }
    counts <- Seurat::Read10X_h5(path)
    if (is.list(counts)) counts <- counts[[1]]
    return(list(type = "10x_h5", counts = counts, obj = NULL))
  }
  if (detect_10x_mtx_dir(path)) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      cat("Package 'Seurat' is required for 10X MTX\n")
      quit(status = 1)
    }
    counts <- Seurat::Read10X(path)
    if (is.list(counts)) counts <- counts[[1]]
    return(list(type = "10x_mtx", counts = counts, obj = NULL))
  }

  if (!(ext %in% c("rds", "qs"))) {
    cat("Unsupported input type for R:", ext, "\n")
    quit(status = 1)
  }

  if (ext == "qs") {
    if (!requireNamespace("qs", quietly = TRUE)) {
      cat("Package 'qs' is required to handle .qs files\n")
      quit(status = 1)
    }
    obj <- qs::qread(path)
  } else {
    obj <- readRDS(path)
  }

  if (inherits(obj, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      cat("Package 'Seurat' is required to read Seurat objects\n")
      quit(status = 1)
    }
    assay_name <- Seurat::DefaultAssay(obj)
    if (is.null(assay_name) || assay_name == "") {
      assay_name <- Seurat::Assays(obj)[[1]]
    }
    counts <- Seurat::GetAssayData(obj, slot = "counts", assay = assay_name)
    return(list(type = "seurat", counts = counts, obj = obj))
  }

  if (inherits(obj, "SummarizedExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      cat("Package 'SummarizedExperiment' is required\n")
      quit(status = 1)
    }
    counts <- SummarizedExperiment::assay(obj)
    return(list(type = "summarized_experiment", counts = counts, obj = obj))
  }

  cat("Unsupported object type in .rds/.qs\n")
  quit(status = 1)
}

loaded <- load_counts(input_path)
counts <- loaded$counts
input_type <- loaded$type
obj <- loaded$obj

n_cells <- ncol(counts)
n_genes <- nrow(counts)
n_nonzero <- if (is_sparse(counts)) Matrix::nnzero(counts) else sum(counts != 0)
sparsity <- 1 - n_nonzero / (n_cells * n_genes)

n_count <- Matrix::colSums(counts)
n_feature <- Matrix::colSums(counts > 0)

qc <- list(
  n_count = summarize_numeric(n_count),
  n_feature = summarize_numeric(n_feature)
)

metadata_fields <- character(0)
assays <- character(0)
reductions <- character(0)
clusters <- character(0)
percent_mt <- NULL

if (!is.null(obj)) {
  if (input_type == "seurat") {
    assays <- Seurat::Assays(obj)
    metadata_fields <- colnames(obj@meta.data)
    reductions <- names(obj@reductions)
    if ("seurat_clusters" %in% metadata_fields) {
      clusters <- "seurat_clusters"
    }
    if ("percent.mt" %in% metadata_fields) {
      percent_mt <- summarize_numeric(obj@meta.data$percent.mt)
    }
  } else if (input_type == "summarized_experiment") {
    if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      metadata_fields <- colnames(SummarizedExperiment::colData(obj))
      assays <- SummarizedExperiment::assayNames(obj)
    }
  }
}

if (!is.null(percent_mt)) {
  qc$percent_mt <- percent_mt
}

output <- list(
  input_path = input_path,
  input_type = input_type,
  engine = "R",
  counts = list(
    n_cells = as.integer(n_cells),
    n_genes = as.integer(n_genes),
    n_nonzero = as.numeric(n_nonzero),
    sparsity = as.numeric(sparsity)
  ),
  qc = qc,
  assays = assays,
  metadata_fields = metadata_fields,
  reductions = reductions,
  clusters = clusters
)

render_table <- function(out) {
  lines <- c(
    "========================================",
    "sc-info (R)",
    "========================================",
    paste0("Input: ", out$input_path),
    paste0("Type: ", out$input_type),
    paste0("Cells: ", out$counts$n_cells),
    paste0("Genes: ", out$counts$n_genes),
    paste0("Nonzero: ", format(out$counts$n_nonzero, scientific = FALSE)),
    paste0("Sparsity: ", sprintf("%.4f", out$counts$sparsity))
  )
  if (length(out$assays) > 0) {
    lines <- c(lines, paste0("Assays: ", paste(out$assays, collapse = ", ")))
  }
  if (length(out$metadata_fields) > 0) {
    lines <- c(lines, paste0("Metadata fields: ", length(out$metadata_fields)))
  }
  if (length(out$reductions) > 0) {
    lines <- c(lines, paste0("Reductions: ", paste(out$reductions, collapse = ", ")))
  }
  if (length(out$clusters) > 0) {
    lines <- c(lines, paste0("Clusters: ", paste(out$clusters, collapse = ", ")))
  }
  lines <- c(lines, "========================================")
  paste(lines, collapse = "\n")
}

table_text <- render_table(output)
cat(table_text, "\n")

if (!is.null(opts$table_out)) {
  writeLines(table_text, opts$table_out)
}

if (!is.null(opts$json_out)) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    cat("Package 'jsonlite' is required for JSON output\n")
    quit(status = 1)
  }
  json <- jsonlite::toJSON(output, auto_unbox = TRUE, pretty = TRUE)
  writeLines(json, opts$json_out)
}

if (!is.null(opts$md_out)) {
  md <- c(
    "# sc-info",
    "",
    paste0("- Input: ", output$input_path),
    paste0("- Type: ", output$input_type),
    paste0("- Engine: ", output$engine),
    "",
    "## Counts",
    "",
    paste0("- Cells: ", output$counts$n_cells),
    paste0("- Genes: ", output$counts$n_genes),
    paste0("- Nonzero: ", format(output$counts$n_nonzero, scientific = FALSE)),
    paste0("- Sparsity: ", sprintf(\"%.4f\", output$counts$sparsity))
  )
  if (length(output$assays) > 0) {
    md <- c(md, "", "## Assays", "", paste0("- ", paste(output$assays, collapse = ", ")))
  }
  if (length(output$metadata_fields) > 0) {
    md <- c(md, "", "## Metadata fields", "", paste0("- ", paste(output$metadata_fields, collapse = ", ")))
  }
  if (length(output$reductions) > 0) {
    md <- c(md, "", "## Reductions", "", paste0("- ", paste(output$reductions, collapse = ", ")))
  }
  if (length(output$clusters) > 0) {
    md <- c(md, "", "## Clusters", "", paste0("- ", paste(output$clusters, collapse = ", ")))
  }
  writeLines(md, opts$md_out)
}
