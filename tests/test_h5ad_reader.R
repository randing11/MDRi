source(file.path("R", "mdri_core.R"))
fixture <- system.file("extdata", "example_anndata.h5ad", package = "zellkonverter")
stopifnot(nzchar(fixture), file.exists(fixture))
mat <- mdri_read_input(fixture, "example_anndata.h5ad")
stopifnot(
  is.matrix(mat) || inherits(mat, "Matrix") || inherits(mat, "DelayedMatrix") || inherits(mat, "DelayedArray"),
  nrow(mat) > 0,
  ncol(mat) > 0
)
cat("Native-R H5AD reader passed:", nrow(mat), "genes x", ncol(mat), "observations\n")
