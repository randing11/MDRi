mdri_gene_sets <- list(
  MDR_injury = c("INHBA", "NLRP3", "IL1B", "CXCL8", "S100A8", "S100A9", "FCN1", "SPP1", "VEGFA", "ADM", "LGALS3", "MMP9", "FN1", "PLAUR", "TGFBI"),
  MDR_resolution = c("FOLR2", "LYVE1", "SELENOP", "MRC1", "STAB1", "APOE", "C1QA", "C1QB", "C1QC", "MERTK", "AXL", "GAS6", "IL4I1", "CD36", "HMOX1", "FTH1", "FTL", "SLC40A1", "CD163"),
  MDR_apc_ifn = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA", "CXCL9", "CXCL10", "ISG15", "IFIT1", "STAT1", "IRF1")
)

mdri_read_expression <- function(path) {
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  gene_col <- which(tolower(colnames(x)) %in% c("gene", "genes", "symbol", "gene_symbol"))[1]
  if (is.na(gene_col)) gene_col <- 1
  genes <- toupper(trimws(x[[gene_col]]))
  mat <- as.matrix(x[, -gene_col, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")
  rownames(mat) <- genes
  mat <- mat[!is.na(rownames(mat)) & rownames(mat) != "", , drop = FALSE]
  mat <- rowsum(mat, group = rownames(mat), reorder = FALSE) / as.vector(table(rownames(mat))[unique(rownames(mat))])
  mat
}

mdri_score <- function(expr, gene_sets = mdri_gene_sets, method = "zscore") {
  genes <- toupper(rownames(expr))
  rownames(expr) <- genes
  expr[is.na(expr)] <- 0

  if (method == "log2_zscore") {
    expr <- log2(expr + 1)
  }
  z <- t(scale(t(expr)))
  z[is.na(z)] <- 0

  scores <- lapply(names(gene_sets), function(nm) {
    gs <- intersect(toupper(gene_sets[[nm]]), rownames(z))
    if (length(gs) == 0) {
      rep(NA_real_, ncol(z))
    } else {
      colMeans(z[gs, , drop = FALSE], na.rm = TRUE)
    }
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
