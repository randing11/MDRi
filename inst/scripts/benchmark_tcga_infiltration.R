suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(survival)
  library(scales)
})

`%||%` <- function(x, y) if (length(x) == 0 || is.na(x) || !nzchar(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
cmd_args <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", cmd_args, value = TRUE)[1])
script_dir <- if (!is.na(file_arg)) dirname(normalizePath(file_arg, mustWork = FALSE)) else getwd()
root_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)

infil_path <- args[1] %||% file.path(root_dir, "inst", "extdata", "benchmark", "infiltration_estimation_for_tcga.csv.gz")
mdri_path <- args[2] %||% file.path(root_dir, "inst", "extdata", "tcga_mdri_scores_survival.csv")
out_dir <- args[3] %||% file.path(root_dir, "benchmark_output")

if (!file.exists(infil_path)) {
  stop(
    "Infiltration file was not found: ", infil_path, "\n",
    "Usage: Rscript inst/scripts/benchmark_tcga_infiltration.R ",
    "<infiltration_estimation_for_tcga.csv.gz> [tcga_mdri_scores_survival.csv] [output_dir]"
  )
}
pdf_dir <- file.path(out_dir, "pdf")
table_dir <- file.path(out_dir, "tables")
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      strip.background = element_rect(fill = "#252525", color = NA),
      strip.text = element_text(color = "white", face = "bold"),
      legend.key.height = unit(0.35, "cm"),
      legend.key.width = unit(0.35, "cm")
    )
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

mdri <- read.csv(mdri_path, check.names = FALSE)
infil <- read.csv(gzfile(infil_path), check.names = FALSE)
colnames(infil)[1] <- "row_id"

mdri$sample15 <- substr(mdri$sample, 1, 15)
infil$sample15 <- substr(infil$barcode, 1, 15)

mdri_features <- c("MDR_injury", "MDR_resolution", "MDR_apc_ifn", "MDR_axis")
benchmark_features <- c(
  "Macrophage_TIMER",
  "Macrophage/Monocyte_MCPCOUNTER",
  "Monocyte_MCPCOUNTER",
  "Myeloid dendritic cell_MCPCOUNTER",
  "Macrophage_XCELL",
  "Macrophage M1_XCELL",
  "Macrophage M2_XCELL",
  "Monocyte_XCELL",
  "Myeloid dendritic cell_XCELL",
  "immune score_XCELL",
  "stroma score_XCELL",
  "microenvironment score_XCELL",
  "Macrophage_EPIC",
  "Macrophage M1_QUANTISEQ",
  "Macrophage M2_QUANTISEQ",
  "Monocyte_QUANTISEQ",
  "Macrophage M0_CIBERSORT",
  "Macrophage M1_CIBERSORT",
  "Macrophage M2_CIBERSORT"
)
benchmark_features <- intersect(benchmark_features, colnames(infil))

label_map <- c(
  MDR_injury = "MDR injury",
  MDR_resolution = "MDR resolution",
  MDR_apc_ifn = "MDR APC/IFN",
  MDR_axis = "MDR axis",
  "Macrophage_TIMER" = "Macrophage\nTIMER",
  "Macrophage/Monocyte_MCPCOUNTER" = "Mac/Mono\nMCPcounter",
  "Monocyte_MCPCOUNTER" = "Monocyte\nMCPcounter",
  "Myeloid dendritic cell_MCPCOUNTER" = "myDC\nMCPcounter",
  "Macrophage_XCELL" = "Macrophage\nxCell",
  "Macrophage M1_XCELL" = "M1 macrophage\nxCell",
  "Macrophage M2_XCELL" = "M2 macrophage\nxCell",
  "Monocyte_XCELL" = "Monocyte\nxCell",
  "Myeloid dendritic cell_XCELL" = "myDC\nxCell",
  "immune score_XCELL" = "Immune score\nxCell",
  "stroma score_XCELL" = "Stroma score\nxCell",
  "microenvironment score_XCELL" = "Microenv score\nxCell",
  "Macrophage_EPIC" = "Macrophage\nEPIC",
  "Macrophage M1_QUANTISEQ" = "M1 macrophage\nquanTIseq",
  "Macrophage M2_QUANTISEQ" = "M2 macrophage\nquanTIseq",
  "Monocyte_QUANTISEQ" = "Monocyte\nquanTIseq",
  "Macrophage M0_CIBERSORT" = "M0 macrophage\nCIBERSORT",
  "Macrophage M1_CIBERSORT" = "M1 macrophage\nCIBERSORT",
  "Macrophage M2_CIBERSORT" = "M2 macrophage\nCIBERSORT"
)

bench <- mdri %>%
  inner_join(infil %>% select(sample15, barcode, all_of(benchmark_features)), by = "sample15")

write.csv(bench, file.path(table_dir, "MDRi_TCGA_infiltration_merged.csv"), row.names = FALSE)

## 1. Correlation benchmark
cor_df <- expand.grid(mdri = mdri_features, benchmark = benchmark_features, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(
    rho = suppressWarnings(cor(bench[[mdri]], bench[[benchmark]], method = "spearman", use = "pairwise.complete.obs")),
    p = suppressWarnings(cor.test(bench[[mdri]], bench[[benchmark]], method = "spearman", exact = FALSE)$p.value)
  ) %>%
  ungroup() %>%
  mutate(
    mdri_label = factor(label_map[mdri], levels = label_map[mdri_features]),
    benchmark_label = factor(label_map[benchmark], levels = rev(label_map[benchmark_features]))
  )
write.csv(cor_df, file.path(table_dir, "MDRi_infiltration_spearman_correlations.csv"), row.names = FALSE)

p_cor <- ggplot(cor_df, aes(x = mdri_label, y = benchmark_label, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.15) +
  scale_fill_gradient2(low = "#3973AC", mid = "white", high = "#B33A3A", midpoint = 0, limits = c(-1, 1), name = "Spearman rho") +
  labs(x = NULL, y = NULL) +
  theme_pub(8) +
  theme(axis.line = element_blank(), axis.ticks = element_blank())
ggsave(file.path(pdf_dir, "Fig8A_MDRi_vs_infiltration_correlation_heatmap.pdf"), p_cor, width = 6.8, height = 7.6, device = cairo_pdf)

## 2. Univariate Cox benchmark
cox_one <- function(dat, feature, endpoint) {
  time_col <- paste0(endpoint, ".time")
  event_col <- endpoint
  d <- dat %>%
    select(cancer, all_of(c(time_col, event_col, feature))) %>%
    rename(time = all_of(time_col), event = all_of(event_col), score = all_of(feature)) %>%
    filter(is.finite(time), is.finite(event), is.finite(score), time > 0)
  if (nrow(d) < 50 || sum(d$event == 1, na.rm = TRUE) < 10 || sd(d$score, na.rm = TRUE) == 0) return(NULL)
  d$score_z <- as.numeric(scale(d$score))
  fit <- tryCatch(coxph(Surv(time, event) ~ score_z, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  data.frame(
    feature = feature,
    endpoint = endpoint,
    n = nrow(d),
    events = sum(d$event == 1, na.rm = TRUE),
    HR = unname(s$coef[1, "exp(coef)"]),
    log2HR = log2(unname(s$coef[1, "exp(coef)"])),
    lower95 = unname(s$conf.int[1, "lower .95"]),
    upper95 = unname(s$conf.int[1, "upper .95"]),
    p = unname(s$coef[1, "Pr(>|z|)"]),
    c_index = unname(s$concordance[1]),
    stringsAsFactors = FALSE
  )
}

all_features <- c(mdri_features, benchmark_features)
endpoints <- c("OS", "DSS", "PFI", "DFI")
cancers <- sort(unique(bench$cancer))
cox_list <- list()
k <- 1
for (ca in cancers) {
  dat_ca <- bench %>% filter(cancer == ca)
  for (ep in endpoints) {
    for (ft in all_features) {
      res <- cox_one(dat_ca, ft, ep)
      if (!is.null(res)) {
        res$cancer <- ca
        cox_list[[k]] <- res
        k <- k + 1
      }
    }
  }
}
cox_df <- bind_rows(cox_list) %>%
  group_by(endpoint) %>%
  mutate(p_adj_endpoint = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  mutate(
    class = ifelse(feature %in% mdri_features, "MDRi", "Benchmark infiltration"),
    feature_label = label_map[feature],
    direction = ifelse(HR >= 1, "Adverse", "Protective")
  )
write.csv(cox_df, file.path(table_dir, "MDRi_infiltration_panTCGA_univariate_Cox.csv"), row.names = FALSE)

sig_count <- cox_df %>%
  filter(p < 0.05) %>%
  count(feature, feature_label, class, direction, name = "n_significant") %>%
  complete(feature, feature_label, class, direction = c("Adverse", "Protective"), fill = list(n_significant = 0)) %>%
  filter(!is.na(feature)) %>%
  mutate(n_plot = ifelse(direction == "Protective", -n_significant, n_significant))

order_feat <- sig_count %>%
  group_by(feature, feature_label, class) %>%
  summarise(total = sum(abs(n_plot)), .groups = "drop") %>%
  arrange(desc(class == "MDRi"), desc(total)) %>%
  pull(feature_label)

p_count <- sig_count %>%
  mutate(feature_label = factor(feature_label, levels = rev(unique(order_feat)))) %>%
  ggplot(aes(x = feature_label, y = n_plot, fill = direction)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.15) +
  coord_flip() +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = c(Adverse = "#B33A3A", Protective = "#3973AC")) +
  scale_y_continuous(labels = abs) +
  labs(x = NULL, y = "Number of nominally significant cancer-endpoint Cox tests", fill = NULL) +
  theme_pub(8) +
  theme(legend.position = "top")
ggsave(file.path(pdf_dir, "Fig8B_MDRi_benchmark_significant_Cox_counts.pdf"), p_count, width = 7.4, height = 6.8, device = cairo_pdf)

## 3. OS Cox heatmap, compact set
compact_features <- c(
  mdri_features,
  "Macrophage_TIMER",
  "Macrophage/Monocyte_MCPCOUNTER",
  "Macrophage_XCELL",
  "Macrophage_EPIC",
  "Macrophage M2_QUANTISEQ",
  "Macrophage M2_CIBERSORT",
  "immune score_XCELL",
  "microenvironment score_XCELL"
)
compact_features <- intersect(compact_features, all_features)
os_heat <- cox_df %>%
  filter(endpoint == "OS", feature %in% compact_features, p < 0.2) %>%
  mutate(
    cancer = factor(cancer, levels = rev(sort(unique(cancer)))),
    feature_label = factor(label_map[feature], levels = label_map[compact_features]),
    p_cap = pmax(p, 1e-5),
    neglog10p = -log10(p_cap)
  )
p_os <- ggplot(os_heat, aes(x = feature_label, y = cancer)) +
  geom_tile(aes(fill = log2HR), color = "white", linewidth = 0.35, width = 0.95, height = 0.9) +
  geom_point(aes(size = neglog10p), shape = 21, fill = "white", color = "black", stroke = 0.45) +
  scale_fill_gradient2(low = "#3973AC", mid = "white", high = "#B33A3A", midpoint = 0, name = "log2(HR)\nper SD") +
  scale_size_continuous(range = c(0.45, 2.2), name = "-log10(p)") +
  labs(x = NULL, y = NULL) +
  theme_pub(7.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.line = element_blank(), axis.ticks = element_blank())
ggsave(file.path(pdf_dir, "Fig8C_OS_Cox_MDRi_vs_benchmark_heatmap.pdf"), p_os, width = 8.8, height = 7.8, device = cairo_pdf)

## 4. Multivariable benchmark: MDRi adjusted for macrophage/immune estimates
adjusters <- intersect(c("Macrophage/Monocyte_MCPCOUNTER", "Macrophage_TIMER", "immune score_XCELL"), colnames(bench))
multi_one <- function(dat, feature, endpoint, adjusters) {
  time_col <- paste0(endpoint, ".time")
  event_col <- endpoint
  cols <- c(time_col, event_col, feature, adjusters)
  d <- dat %>%
    select(all_of(cols)) %>%
    rename(time = all_of(time_col), event = all_of(event_col), score = all_of(feature)) %>%
    filter(if_all(everything(), ~is.finite(.x)), time > 0)
  if (nrow(d) < 60 || sum(d$event == 1) < 10 || sd(d$score) == 0) return(NULL)
  d$score_z <- as.numeric(scale(d$score))
  for (adj in adjusters) d[[safe_name(adj)]] <- as.numeric(scale(d[[adj]]))
  form <- as.formula(paste("Surv(time, event) ~ score_z +", paste(safe_name(adjusters), collapse = " + ")))
  fit <- tryCatch(coxph(form, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  data.frame(
    feature = feature,
    endpoint = endpoint,
    n = nrow(d),
    events = sum(d$event == 1),
    HR = unname(s$coef["score_z", "exp(coef)"]),
    log2HR = log2(unname(s$coef["score_z", "exp(coef)"])),
    lower95 = unname(s$conf.int["score_z", "lower .95"]),
    upper95 = unname(s$conf.int["score_z", "upper .95"]),
    p = unname(s$coef["score_z", "Pr(>|z|)"]),
    stringsAsFactors = FALSE
  )
}

multi_list <- list()
k <- 1
for (ca in cancers) {
  dat_ca <- bench %>% filter(cancer == ca)
  for (ep in c("OS", "PFI")) {
    for (ft in mdri_features) {
      res <- multi_one(dat_ca, ft, ep, adjusters)
      if (!is.null(res)) {
        res$cancer <- ca
        multi_list[[k]] <- res
        k <- k + 1
      }
    }
  }
}
multi_df <- bind_rows(multi_list) %>%
  mutate(feature_label = label_map[feature], neglog10p = -log10(pmax(p, 1e-5)))
write.csv(multi_df, file.path(table_dir, "MDRi_multivariable_adjusted_for_benchmark_infiltration.csv"), row.names = FALSE)

p_multi <- multi_df %>%
  filter(p < 0.2) %>%
  mutate(
    cancer = factor(cancer, levels = rev(sort(unique(cancer)))),
    feature_label = factor(feature_label, levels = label_map[mdri_features])
  ) %>%
  ggplot(aes(x = feature_label, y = cancer)) +
  geom_tile(aes(fill = log2HR), color = "white", linewidth = 0.35, width = 0.95, height = 0.9) +
  geom_point(aes(size = neglog10p), shape = 21, fill = "white", color = "black", stroke = 0.45) +
  facet_wrap(~endpoint, nrow = 1) +
  scale_fill_gradient2(low = "#3973AC", mid = "white", high = "#B33A3A", midpoint = 0, name = "adjusted\nlog2(HR)") +
  scale_size_continuous(range = c(0.45, 2.2), name = "-log10(p)") +
  labs(x = NULL, y = NULL, caption = paste("Adjusted for:", paste(label_map[adjusters], collapse = ", "))) +
  theme_pub(7.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.line = element_blank(), axis.ticks = element_blank())
ggsave(file.path(pdf_dir, "Fig8D_MDRi_multivariable_adjusted_Cox_heatmap.pdf"), p_multi, width = 8.4, height = 6.8, device = cairo_pdf)

## Summary text
summary_lines <- c(
  paste0("Merged samples: ", nrow(bench)),
  paste0("Cancers: ", paste(sort(unique(bench$cancer)), collapse = ", ")),
  paste0("Benchmark features: ", paste(benchmark_features, collapse = "; ")),
  paste0("Adjusters in multivariable Cox: ", paste(adjusters, collapse = "; "))
)
writeLines(summary_lines, file.path(out_dir, "benchmark_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("DONE MDRi benchmark\n")
cat("Output:", out_dir, "\n")
