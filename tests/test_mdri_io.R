source(file.path("R", "mdri_core.R"))

set.seed(20260805)
genes <- unique(unlist(mdri_gene_sets, use.names = FALSE))
mat <- matrix(rnorm(length(genes) * 24), nrow = length(genes),
              dimnames = list(genes, paste0("sample_", seq_len(24))))

expected_z <- t(scale(t(mat)))
expected <- vapply(mdri_gene_sets, function(gs) {
  colMeans(expected_z[intersect(gs, rownames(expected_z)), , drop = FALSE])
}, numeric(ncol(mat)))
observed <- mdri_score(mat)[, names(mdri_gene_sets), drop = FALSE]
stopifnot(isTRUE(all.equal(as.matrix(observed), expected, tolerance = 1e-10, check.attributes = FALSE)))

rds_path <- tempfile(fileext = ".rds")
saveRDS(mat, rds_path)
from_rds <- mdri_read_input(rds_path, "matrix.rds")
stopifnot(isTRUE(all.equal(from_rds, mat, check.attributes = TRUE)))

csv_path <- tempfile(fileext = ".csv")
utils::write.csv(data.frame(gene = rownames(mat), mat, check.names = FALSE), csv_path, row.names = FALSE)
from_csv <- mdri_read_input(csv_path, "matrix.csv")
stopifnot(isTRUE(all.equal(unname(from_csv), unname(mat), tolerance = 1e-10)))

if (requireNamespace("Matrix", quietly = TRUE)) {
  sparse <- Matrix::Matrix(mat, sparse = TRUE)
  sparse_score <- mdri_score(sparse)[, names(mdri_gene_sets), drop = FALSE]
  stopifnot(isTRUE(all.equal(as.matrix(sparse_score), expected, tolerance = 1e-10, check.attributes = FALSE)))
}

if (requireNamespace("SeuratObject", quietly = TRUE) && requireNamespace("Matrix", quietly = TRUE)) {
  counts <- Matrix::Matrix(round(abs(mat) * 10), sparse = TRUE)
  seurat_object <- SeuratObject::CreateSeuratObject(counts = counts)
  seurat_path <- tempfile(fileext = ".rds")
  saveRDS(seurat_object, seurat_path)
  from_seurat <- mdri_read_input(seurat_path, "seurat.rds")
  stopifnot(identical(dim(from_seurat), dim(counts)))
}

if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  sce <- SingleCellExperiment::SingleCellExperiment(list(logcounts = mat))
  sce_path <- tempfile(fileext = ".rds")
  saveRDS(sce, sce_path)
  from_sce <- mdri_read_input(sce_path, "sce.rds")
  stopifnot(isTRUE(all.equal(from_sce, mat, check.attributes = TRUE)))

  if (requireNamespace("zellkonverter", quietly = TRUE)) {
    h5ad_path <- tempfile(fileext = ".h5ad")
    h5ad_written <- tryCatch({
      zellkonverter::writeH5AD(sce, h5ad_path, X_name = "logcounts")
      TRUE
    }, error = function(e) {
      message("H5AD round-trip skipped because the optional Python environment is unavailable: ", conditionMessage(e))
      FALSE
    })
    if (h5ad_written) {
      from_h5ad <- mdri_read_input(h5ad_path, "sce.h5ad")
      stopifnot(identical(dim(from_h5ad), dim(mat)))
    }
  }
}

optional <- c("SeuratObject", "SummarizedExperiment", "SingleCellExperiment", "zellkonverter")
cat("Optional packages:\n")
print(setNames(vapply(optional, requireNamespace, logical(1), quietly = TRUE), optional))
cat("MDRi input and scoring tests passed.\n")
