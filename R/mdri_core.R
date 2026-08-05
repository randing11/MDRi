`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x[1]) || !nzchar(x[1])) y else x

mdri_gene_sets <- list(
  MDR_injury = c("INHBA", "NLRP3", "IL1B", "CXCL8", "S100A8", "S100A9", "FCN1", "SPP1", "VEGFA", "ADM", "LGALS3", "MMP9", "FN1", "PLAUR", "TGFBI"),
  MDR_resolution = c("FOLR2", "LYVE1", "SELENOP", "MRC1", "STAB1", "APOE", "C1QA", "C1QB", "C1QC", "MERTK", "AXL", "GAS6", "IL4I1", "CD36", "HMOX1", "FTH1", "FTL", "SLC40A1", "CD163"),
  MDR_apc_ifn = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA", "CXCL9", "CXCL10", "ISG15", "IFIT1", "STAT1", "IRF1")
)

mdri_prepare_expression <- function(x, gene_sets = mdri_gene_sets) {
  if (is.data.frame(x)) {
    gene_col <- which(tolower(colnames(x)) %in% c("gene", "genes", "symbol", "gene_symbol"))[1]
    if (!is.na(gene_col)) {
      genes <- x[[gene_col]]
      x <- as.matrix(x[, -gene_col, drop = FALSE])
      rownames(x) <- genes
    } else {
      x <- as.matrix(x)
    }
  }
  is_delayed <- inherits(x, "DelayedMatrix") || inherits(x, "DelayedArray")
  if (!is.matrix(x) && !inherits(x, "Matrix") && !is_delayed) {
    stop("Expression input must resolve to a matrix-like object.")
  }
  if (is.null(rownames(x))) stop("Expression input must contain gene names as row names.")

  target_genes <- unique(toupper(unlist(gene_sets, use.names = FALSE)))
  row_hits <- sum(toupper(rownames(x)) %in% target_genes)
  col_hits <- if (is.null(colnames(x))) 0 else sum(toupper(colnames(x)) %in% target_genes)
  if (col_hits > row_hits) x <- t(x)
  if (is.null(colnames(x))) colnames(x) <- paste0("sample_", seq_len(ncol(x)))

  rownames(x) <- toupper(trimws(rownames(x)))
  keep <- !is.na(rownames(x)) & nzchar(rownames(x))
  x <- x[keep, , drop = FALSE]
  if (inherits(x, "sparseMatrix")) {
    x@x <- as.numeric(x@x)
  } else if (!is_delayed) {
    suppressWarnings(storage.mode(x) <- "numeric")
  }

  if (anyDuplicated(rownames(x))) {
    genes <- rownames(x)
    gene_levels <- unique(genes)
    group_id <- match(genes, gene_levels)
    group_n <- tabulate(group_id, nbins = length(gene_levels))
    if (is_delayed) {
      x <- x[!duplicated(rownames(x)), , drop = FALSE]
    } else if (inherits(x, "sparseMatrix") && requireNamespace("Matrix", quietly = TRUE)) {
      aggregator <- Matrix::sparseMatrix(
        i = group_id,
        j = seq_along(group_id),
        x = 1 / group_n[group_id],
        dims = c(length(gene_levels), length(group_id))
      )
      x <- aggregator %*% x
      rownames(x) <- gene_levels
    } else {
      x <- rowsum(as.matrix(x), group = genes, reorder = FALSE)
      x <- x / group_n[match(rownames(x), gene_levels)]
    }
  }
  x
}

mdri_read_expression <- function(path) {
  ext <- tolower(tools::file_ext(path))
  sep <- if (ext == "csv") "," else "\t"
  x <- utils::read.table(
    path,
    header = TRUE,
    sep = sep,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = if (ext == "csv") "\"" else "",
    comment.char = ""
  )
  mdri_prepare_expression(x)
}

mdri_extract_r_object <- function(object, assay = NULL, layer = NULL) {
  if (is.matrix(object) || inherits(object, "Matrix") || is.data.frame(object)) {
    return(mdri_prepare_expression(object))
  }

  if (inherits(object, "Seurat")) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("Reading a Seurat RDS requires the SeuratObject package.")
    }
    assay <- assay %||% SeuratObject::DefaultAssay(object)
    available_layers <- tryCatch(SeuratObject::Layers(object[[assay]]), error = function(e) character())
    if (is.null(layer)) {
      layer <- intersect(c("data", "counts", "scale.data"), available_layers)[1]
      if (is.na(layer) || !nzchar(layer)) layer <- "data"
    }
    mat <- tryCatch(
      SeuratObject::GetAssayData(object, assay = assay, layer = layer),
      error = function(e) SeuratObject::GetAssayData(object, assay = assay, slot = layer)
    )
    return(mdri_prepare_expression(mat))
  }

  if (inherits(object, "SingleCellExperiment") || inherits(object, "SummarizedExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("Reading a SingleCellExperiment RDS requires SummarizedExperiment.")
    }
    available <- SummarizedExperiment::assayNames(object)
    assay_name <- assay %||% intersect(c("logcounts", "data", "counts"), available)[1]
    if (is.na(assay_name) || !nzchar(assay_name)) assay_name <- available[1]
    if (is.na(assay_name) || !nzchar(assay_name)) stop("No assay was found in the uploaded object.")
    return(mdri_prepare_expression(SummarizedExperiment::assay(object, assay_name)))
  }

  stop("Unsupported RDS object. Use a matrix, data.frame, Seurat, SingleCellExperiment, or SummarizedExperiment object.")
}

mdri_read_input <- function(path, filename = basename(path), assay = NULL, layer = NULL) {
  ext <- tolower(tools::file_ext(filename))
  if (ext %in% c("tsv", "txt", "csv")) return(mdri_read_expression(path))
  if (ext == "rds") return(mdri_extract_r_object(readRDS(path), assay = assay, layer = layer))
  if (ext == "h5ad") {
    if (!requireNamespace("zellkonverter", quietly = TRUE)) {
      stop("Reading H5AD requires the optional zellkonverter package.")
    }
    sce <- zellkonverter::readH5AD(path, use_hdf5 = TRUE, reader = "R")
    return(mdri_extract_r_object(sce, assay = assay, layer = layer))
  }
  stop("Unsupported file type: .", ext)
}

mdri_score <- function(expr, gene_sets = mdri_gene_sets, method = "zscore") {
  expr <- mdri_prepare_expression(expr, gene_sets = gene_sets)
  available <- intersect(unique(toupper(unlist(gene_sets, use.names = FALSE))), rownames(expr))
  expr <- expr[available, , drop = FALSE]
  if (method == "log2_zscore") {
    if (inherits(expr, "sparseMatrix")) {
      expr@x <- log2(expr@x + 1)
    } else {
      expr <- log2(expr + 1)
    }
  }

  row_mean <- if (inherits(expr, "Matrix")) Matrix::rowMeans(expr, na.rm = TRUE) else rowMeans(expr, na.rm = TRUE)
  n_obs <- ncol(expr)
  if (n_obs > 1) {
    row_sq <- if (inherits(expr, "Matrix")) Matrix::rowSums(expr ^ 2, na.rm = TRUE) else rowSums(expr ^ 2, na.rm = TRUE)
    row_sd <- sqrt(pmax((row_sq - n_obs * row_mean ^ 2) / (n_obs - 1), 0))
  } else {
    row_sd <- rep(NA_real_, nrow(expr))
  }
  names(row_mean) <- names(row_sd) <- rownames(expr)

  scores <- lapply(names(gene_sets), function(nm) {
    gs <- intersect(toupper(gene_sets[[nm]]), rownames(expr))
    gs <- gs[is.finite(row_sd[gs]) & row_sd[gs] > 0]
    if (length(gs) == 0) return(rep(NA_real_, ncol(expr)))
    expr_gs <- expr[gs, , drop = FALSE]
    if (inherits(expr_gs, "sparseMatrix")) {
      scaled_sum <- Matrix::colSums(Matrix::Diagonal(x = 1 / row_sd[gs]) %*% expr_gs, na.rm = TRUE)
    } else if (inherits(expr_gs, "DelayedMatrix") || inherits(expr_gs, "DelayedArray")) {
      scaled_sum <- colSums(sweep(expr_gs, 1, row_sd[gs], "/"), na.rm = TRUE)
    } else {
      scaled_sum <- colSums(sweep(as.matrix(expr_gs), 1, row_sd[gs], "/"), na.rm = TRUE)
    }
    as.numeric(scaled_sum / length(gs) - mean(row_mean[gs] / row_sd[gs]))
  })
  scores <- as.data.frame(scores, check.names = FALSE)
  colnames(scores) <- names(gene_sets)
  scores$MDR_axis <- scores$MDR_injury - scores$MDR_resolution
  scores$sample <- colnames(expr)
  scores$patient <- substr(scores$sample, 1, 12)
  scores[, c("sample", "patient", names(gene_sets), "MDR_axis")]
}

mdri_gene_coverage <- function(expr, gene_sets = mdri_gene_sets) {
  genes <- toupper(rownames(expr))
  do.call(rbind, lapply(names(gene_sets), function(nm) {
    gs <- toupper(gene_sets[[nm]])
    data.frame(
      program = nm,
      n_total = length(gs),
      n_present = sum(gs %in% genes),
      coverage = sum(gs %in% genes) / length(gs),
      missing = paste(setdiff(gs, genes), collapse = ", ")
    )
  }))
}

mdri_read_clinical <- function(path) {
  clin <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"sample" %in% colnames(clin)) {
    colnames(clin)[1] <- "sample"
  }
  clin$sample <- as.character(clin$sample)
  clin$patient <- if ("patient" %in% colnames(clin)) as.character(clin$patient) else substr(clin$sample, 1, 12)
  clin
}

mdri_merge_clinical <- function(scores, clinical) {
  by_sample <- merge(scores, clinical, by = "sample", all = FALSE)
  if (nrow(by_sample) > 0) return(by_sample)
  clinical <- clinical[!duplicated(clinical$patient), , drop = FALSE]
  merge(scores, clinical, by = "patient", all = FALSE, suffixes = c("", ".clin"))
}

mdri_cox <- function(dat, score, time_col, event_col) {
  stopifnot(requireNamespace("survival", quietly = TRUE))
  use <- dat[, c(score, time_col, event_col), drop = FALSE]
  colnames(use) <- c("score", "time", "event")
  use <- use[complete.cases(use) & use$time > 0, , drop = FALSE]
  if (nrow(use) < 20 || sum(use$event == 1) < 5) return(NULL)
  use$score_z <- as.numeric(scale(use$score))
  fit <- survival::coxph(survival::Surv(time, event) ~ score_z, data = use)
  s <- summary(fit)
  data.frame(
    score = score,
    n = nrow(use),
    events = sum(use$event == 1),
    HR_per_SD = unname(s$coefficients[1, "exp(coef)"]),
    CI_low = unname(s$conf.int[1, "lower .95"]),
    CI_high = unname(s$conf.int[1, "upper .95"]),
    p = unname(s$coefficients[1, "Pr(>|z|)"])
  )
}

mdri_km_data <- function(dat, score, time_col, event_col) {
  stopifnot(requireNamespace("survival", quietly = TRUE))
  use <- dat[, c(score, time_col, event_col), drop = FALSE]
  colnames(use) <- c("score", "time", "event")
  use <- use[complete.cases(use) & use$time > 0, , drop = FALSE]
  if (nrow(use) < 20 || sum(use$event == 1) < 5) return(NULL)
  use$group <- factor(ifelse(use$score >= median(use$score, na.rm = TRUE), "High", "Low"), levels = c("High", "Low"))
  fit <- survival::survfit(survival::Surv(time, event) ~ group, data = use)
  sf <- summary(fit)
  sdat <- data.frame(time = sf$time, survival = sf$surv, group = sub("^group=", "", sf$strata))
  sdif <- survival::survdiff(survival::Surv(time, event) ~ group, data = use)
  p <- 1 - pchisq(sdif$chisq, df = 1)
  list(curve = sdat, p = p, n = nrow(use), events = sum(use$event == 1))
}

mdri_plot_heatmap <- function(scores) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  score_cols <- c("MDR_injury", "MDR_resolution", "MDR_apc_ifn", "MDR_axis")
  x <- scores[, c("sample", score_cols), drop = FALSE]
  long <- reshape(x, varying = score_cols, v.names = "score_value", timevar = "program",
                  times = score_cols, direction = "long")
  long$program <- factor(long$program, levels = score_cols)
  ggplot2::ggplot(long, ggplot2::aes(sample, program, fill = score_value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#b2182b", midpoint = 0) +
    ggplot2::labs(x = NULL, y = NULL, fill = "MDRi score") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
}

mdri_plot_scatter <- function(scores) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  ggplot2::ggplot(scores, ggplot2::aes(MDR_resolution, MDR_injury, color = MDR_axis)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
    ggplot2::geom_point(size = 2.2, alpha = 0.85) +
    ggplot2::scale_color_gradient2(low = "#3674a2", mid = "white", high = "#b95055", midpoint = 0) +
    ggplot2::labs(x = "MDR resolution", y = "MDR injury", color = "Axis") +
    ggplot2::theme_classic(base_size = 12)
}
