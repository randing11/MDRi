suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(ggplot2)
  library(dplyr)
  library(survival)
})

options(shiny.maxRequestSize = 2048 * 1024^2)

`%||%` <- function(x, y) if (is.null(x)) y else x

app_dir <- normalizePath(getwd(), mustWork = FALSE)
if (basename(app_dir) == "shiny") {
  root_dir <- normalizePath(file.path(app_dir, "..", ".."), mustWork = FALSE)
} else {
  root_dir <- normalizePath(file.path(getwd()), mustWork = FALSE)
}
source(file.path(root_dir, "R", "mdri_core.R"))

ext_dir <- file.path(root_dir, "inst", "extdata")
ref_path <- file.path(ext_dir, "tcga_mdri_reference.rds")
reference <- if (file.exists(ref_path)) readRDS(ref_path) else NULL
benchmark_dir <- file.path(ext_dir, "benchmark")
benchmark <- list(
  summary = file.path(benchmark_dir, "benchmark_summary.txt"),
  cor = file.path(benchmark_dir, "MDRi_infiltration_spearman_correlations.csv"),
  cox = file.path(benchmark_dir, "MDRi_infiltration_panTCGA_univariate_Cox.csv"),
  multi = file.path(benchmark_dir, "MDRi_multivariable_adjusted_for_benchmark_infiltration.csv"),
  tam_cor = file.path(benchmark_dir, "MDRi_published_TAM_signature_global_correlations.tsv"),
  tam_cox = file.path(benchmark_dir, "published_TAM_signature_panTCGA_univariate_Cox.tsv"),
  tam_multi = file.path(benchmark_dir, "MDRi_Cox_adjusted_for_TAM_age_stage_purity_CD8.tsv"),
  tam_genes = file.path(benchmark_dir, "published_TAM_signature_gene_sets.tsv")
)
benchmark <- lapply(benchmark, function(path) {
  if (!file.exists(path)) return(NULL)
  if (grepl("\\.txt$", path)) return(readLines(path, warn = FALSE))
  if (grepl("\\.tsv$", path)) return(read.delim(path, check.names = FALSE))
  read.csv(path, check.names = FALSE)
})

score_cols <- c("MDR_injury", "MDR_resolution", "MDR_apc_ifn", "MDR_axis")

ui <- page_navbar(
  title = "MDRi Explorer",
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#2d3462"),

  nav_panel(
    "Overview",
    layout_columns(
      col_widths = c(4, 4, 4),
      card(
        card_header("MDRi framework"),
        tags$p("MDRi quantifies pan-cancer myeloid damage-response states from bulk or single-cell-derived expression matrices."),
        tags$ul(
          tags$li(tags$b("MDR injury:"), " inflammasome, hypoxia, inflammatory remodeling."),
          tags$li(tags$b("MDR resolution:"), " lipid/efferocytosis, resident repair, iron/heme handling."),
          tags$li(tags$b("MDR APC/IFN:"), " antigen presentation and interferon activation."),
          tags$li(tags$b("MDR axis:"), " injury minus resolution.")
        ),
        tags$p("MVP version: upload expression, compute MDRi scores, compare with TCGA reference, and run survival analysis when clinical data are available.")
      ),
      card(
        card_header("Benchmark support"),
        tags$p("The packaged reference includes benchmark analyses against TIMER, MCPcounter, xCell, CIBERSORT, quanTIseq, EPIC and published C1QC+ and SPP1+ TAM signatures."),
        tags$p("Benchmark results are available in the Benchmark tab and can be regenerated from inst/scripts/benchmark_tcga_infiltration.R.")
      ),
      card(
        card_header("Built-in reference"),
        verbatimTextOutput("reference_summary")
      )
    )
  ),

  nav_panel(
    "Score Data",
    layout_sidebar(
      sidebar = sidebar(
        fileInput(
          "expr_file",
          "Expression data",
          accept = c(".tsv", ".txt", ".csv", ".rds", ".h5ad")
        ),
        checkboxInput("use_example", "Use built-in example expression", TRUE),
        textInput("input_assay", "Assay (optional)", value = ""),
        textInput("input_layer", "Layer (optional)", value = ""),
        selectInput("score_method", "Scoring method", choices = c("zscore", "log2_zscore"), selected = "zscore"),
        actionButton("run_score", "Run MDRi scoring", class = "btn-primary"),
        hr(),
        downloadButton("download_scores", "Download scores")
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("Gene coverage"), DTOutput("coverage_table")),
        card(card_header("MDRi score table"), DTOutput("score_table"))
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("Score heatmap"), plotOutput("score_heatmap", height = "360px")),
        card(card_header("Injury-resolution scatter"), plotOutput("score_scatter", height = "360px"))
      )
    )
  ),

  nav_panel(
    "TCGA Reference",
    layout_sidebar(
      sidebar = sidebar(
        selectInput("ref_cancer", "Cancer type", choices = "All"),
        selectInput("ref_score", "MDRi score", choices = score_cols, selected = "MDR_axis")
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("TCGA score distribution"), plotOutput("ref_boxplot", height = "380px")),
        card(card_header("Pan-cancer Cox atlas"), plotOutput("cox_atlas", height = "380px"))
      ),
      card(card_header("TCGA reference score table"), DTOutput("ref_table"))
    )
  ),

  nav_panel(
    "Clinical Survival",
    layout_sidebar(
      sidebar = sidebar(
        fileInput("clin_file", "Clinical TSV/CSV", accept = c(".tsv", ".txt", ".csv")),
        checkboxInput("use_example_clin", "Use built-in example clinical table", TRUE),
        selectInput("surv_score", "Score", choices = score_cols, selected = "MDR_axis"),
        textInput("time_col", "Time column", value = "OS.time"),
        textInput("event_col", "Event column", value = "OS"),
        actionButton("run_surv", "Run survival", class = "btn-primary")
      ),
      layout_columns(
        col_widths = c(5, 7),
        card(card_header("Cox result"), DTOutput("cox_table")),
        card(card_header("Kaplan-Meier"), plotOutput("km_plot", height = "380px"))
      )
    )
  ),

  nav_panel(
    "Benchmark",
    layout_columns(
      col_widths = c(4, 8),
      card(
        card_header("Benchmark summary"),
        verbatimTextOutput("benchmark_summary")
      ),
      card(
        card_header("MDRi versus infiltration estimates"),
        plotOutput("benchmark_cor_plot", height = "520px")
      )
    ),
    layout_columns(
      col_widths = c(6, 6),
      card(card_header("Nominally significant Cox tests"), plotOutput("benchmark_count_plot", height = "460px")),
      card(card_header("MDRi adjusted Cox table"), DTOutput("benchmark_multi_table"))
    ),
    card(card_header("Benchmark Cox result table"), DTOutput("benchmark_cox_table")),
    layout_columns(
      col_widths = c(5, 7),
      card(card_header("Published TAM signature correlations"), plotOutput("tam_cor_plot", height = "330px")),
      card(card_header("Endpoint-FDR significant models"), plotOutput("tam_count_plot", height = "330px"))
    ),
    layout_columns(
      col_widths = c(5, 7),
      card(card_header("Published TAM gene sets"), DTOutput("tam_gene_table")),
      card(card_header("MDRi adjusted for TAM and clinical covariates"), DTOutput("tam_multi_table"))
    )
  ),

  nav_panel(
    "Downloads",
    card(
      card_header("Files included in this MDRi MVP"),
      DTOutput("file_table")
    )
  )
)

server <- function(input, output, session) {
  observe({
    if (!is.null(reference)) {
      choices <- c("All", sort(unique(reference$tcga_scores$cancer)))
      updateSelectInput(session, "ref_cancer", choices = choices, selected = "All")
    }
  })

  output$reference_summary <- renderText({
    if (is.null(reference)) return("No TCGA reference loaded.")
    paste(
      "TCGA samples:", nrow(reference$tcga_scores),
      "\nCancer types:", length(unique(reference$tcga_scores$cancer)),
      "\nCox results:", nrow(reference$tcga_cox),
      "\nKM significant entries:", nrow(reference$tcga_km),
      "\nCreated:", reference$created,
      "\nNote:", reference$note
    )
  })

  output$benchmark_summary <- renderText({
    if (is.null(benchmark$summary)) {
      return("No benchmark summary file was found.")
    }
    paste(benchmark$summary, collapse = "\n")
  })

  output$benchmark_cor_plot <- renderPlot({
    req(benchmark$cor)
    x <- benchmark$cor
    x$mdri_label <- factor(x$mdri_label, levels = unique(x$mdri_label))
    x$benchmark_label <- factor(x$benchmark_label, levels = rev(unique(x$benchmark_label)))
    ggplot(x, aes(mdri_label, benchmark_label, fill = rho)) +
      geom_tile(color = "white", linewidth = 0.25) +
      geom_text(aes(label = sprintf("%.2f", rho)), size = 2.2) +
      scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-1, 1)) +
      labs(x = NULL, y = NULL, fill = "Spearman rho") +
      theme_classic(base_size = 9) +
      theme(axis.line = element_blank(), axis.ticks = element_blank())
  })

  output$benchmark_count_plot <- renderPlot({
    req(benchmark$cox)
    x <- benchmark$cox
    x <- x[is.finite(x$p) & x$p < 0.05, , drop = FALSE]
    if (nrow(x) == 0) return(NULL)
    count_df <- x |>
      group_by(feature_label, class, direction) |>
      summarise(n_significant = n(), .groups = "drop") |>
      mutate(n_plot = ifelse(direction == "Protective", -n_significant, n_significant))
    order_df <- count_df |>
      group_by(feature_label, class) |>
      summarise(total = sum(abs(n_plot)), .groups = "drop") |>
      arrange(desc(class == "MDRi"), desc(total))
    count_df$feature_label <- factor(count_df$feature_label, levels = rev(unique(order_df$feature_label)))
    ggplot(count_df, aes(feature_label, n_plot, fill = direction)) +
      geom_col(width = 0.72, color = "white", linewidth = 0.15) +
      coord_flip() +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.25) +
      scale_fill_manual(values = c(Adverse = "#b2182b", Protective = "#2b6cb0")) +
      scale_y_continuous(labels = abs) +
      labs(x = NULL, y = "Significant cancer-endpoint Cox tests", fill = NULL) +
      theme_classic(base_size = 9) +
      theme(legend.position = "top")
  })

  output$benchmark_multi_table <- renderDT({
    req(benchmark$multi)
    x <- benchmark$multi
    x <- x[order(x$p), c("cancer", "endpoint", "feature_label", "n", "events", "HR", "lower95", "upper95", "p"), drop = FALSE]
    datatable(x, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$benchmark_cox_table <- renderDT({
    req(benchmark$cox)
    x <- benchmark$cox
    x <- x[order(x$p), c("class", "feature_label", "cancer", "endpoint", "n", "events", "HR", "lower95", "upper95", "p", "c_index"), drop = FALSE]
    datatable(x, options = list(pageLength = 15, scrollX = TRUE))
  })

  output$tam_cor_plot <- renderPlot({
    req(benchmark$tam_cor)
    x <- benchmark$tam_cor
    x$mdri <- factor(
      x$mdri,
      levels = score_cols,
      labels = c("MDR injury", "MDR resolution", "MDR APC/IFN", "Injury-resolution axis")
    )
    x$tam_signature <- factor(
      x$tam_signature,
      levels = c("C1QC_TAM", "SPP1_TAM"),
      labels = c("C1QC+ TAM", "SPP1+ TAM")
    )
    ggplot(x, aes(mdri, tam_signature, fill = rho)) +
      geom_tile(color = "white", linewidth = 0.35) +
      geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
      scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-1, 1)) +
      labs(x = NULL, y = NULL, fill = "Spearman rho") +
      theme_classic(base_size = 9) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())
  })

  output$tam_count_plot <- renderPlot({
    req(benchmark$tam_cox, benchmark$tam_multi)
    tam <- benchmark$tam_cox |>
      mutate(class = "Published TAM signatures")
    mdri <- benchmark$tam_multi |>
      mutate(class = "MDRi adjusted for TAM and clinical covariates")
    x <- bind_rows(tam, mdri) |>
      group_by(class, feature, direction) |>
      summarise(n_FDR = sum(q_endpoint < 0.05), .groups = "drop") |>
      mutate(
        n_plot = ifelse(direction == "Protective", -n_FDR, n_FDR),
        feature = recode(
          feature,
          MDR_injury = "MDR injury",
          MDR_resolution = "MDR resolution",
          MDR_apc_ifn = "MDR APC/IFN",
          MDR_axis = "Injury-resolution axis",
          C1QC_TAM = "C1QC+ TAM",
          SPP1_TAM = "SPP1+ TAM"
        )
      )
    ggplot(x, aes(feature, n_plot, fill = direction)) +
      geom_col(width = 0.68) +
      coord_flip() +
      facet_wrap(~class, scales = "free_y", ncol = 1) +
      geom_hline(yintercept = 0, linewidth = 0.25) +
      scale_fill_manual(values = c(Adverse = "#b2182b", Protective = "#2b6cb0")) +
      scale_y_continuous(labels = abs) +
      labs(x = NULL, y = "Endpoint-FDR significant models", fill = NULL) +
      theme_classic(base_size = 8) +
      theme(legend.position = "top", strip.text = element_text(size = 7))
  })

  output$tam_gene_table <- renderDT({
    req(benchmark$tam_genes)
    datatable(benchmark$tam_genes, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$tam_multi_table <- renderDT({
    req(benchmark$tam_multi)
    x <- benchmark$tam_multi
    x <- x[order(x$q_endpoint, x$p), c("cancer", "endpoint", "feature", "n", "events", "HR", "lower95", "upper95", "p", "q_endpoint", "stage_adjusted"), drop = FALSE]
    datatable(x, options = list(pageLength = 10, scrollX = TRUE))
  })

  expression_matrix <- eventReactive(input$run_score, {
    if (isTRUE(input$use_example)) {
      mdri_read_expression(file.path(ext_dir, "example_bulk_expression_MDR_genes.tsv"))
    } else {
      req(input$expr_file)
      mdri_read_input(
        input$expr_file$datapath,
        filename = input$expr_file$name,
        assay = if (nzchar(trimws(input$input_assay))) trimws(input$input_assay) else NULL,
        layer = if (nzchar(trimws(input$input_layer))) trimws(input$input_layer) else NULL
      )
    }
  })

  scores <- eventReactive(input$run_score, {
    expr <- expression_matrix()
    mdri_score(expr, method = input$score_method)
  })

  output$coverage_table <- renderDT({
    req(expression_matrix())
    datatable(mdri_gene_coverage(expression_matrix()), options = list(pageLength = 5, scrollX = TRUE))
  })

  output$score_table <- renderDT({
    req(scores())
    datatable(scores(), options = list(pageLength = 10, scrollX = TRUE))
  })

  output$score_heatmap <- renderPlot({
    req(scores())
    mdri_plot_heatmap(scores())
  })

  output$score_scatter <- renderPlot({
    req(scores())
    mdri_plot_scatter(scores())
  })

  output$download_scores <- downloadHandler(
    filename = function() "MDRi_scores.csv",
    content = function(file) write.csv(scores(), file, row.names = FALSE)
  )

  ref_filtered <- reactive({
    req(reference)
    x <- reference$tcga_scores
    if (!is.null(input$ref_cancer) && input$ref_cancer != "All") x <- x[x$cancer == input$ref_cancer, , drop = FALSE]
    x
  })

  output$ref_boxplot <- renderPlot({
    req(reference, input$ref_score)
    x <- reference$tcga_scores
    top_cancers <- x |>
      dplyr::group_by(cancer) |>
      dplyr::summarise(n = dplyr::n(), med = median(.data[[input$ref_score]], na.rm = TRUE), .groups = "drop") |>
      dplyr::filter(n >= 60) |>
      dplyr::arrange(desc(med)) |>
      dplyr::slice_head(n = 16) |>
      dplyr::pull(cancer)
    x <- x[x$cancer %in% top_cancers, , drop = FALSE]
    x$cancer <- factor(x$cancer, levels = top_cancers)
    ggplot(x, aes(cancer, .data[[input$ref_score]], fill = cancer)) +
      geom_boxplot(outlier.shape = NA, linewidth = 0.25) +
      geom_jitter(width = 0.18, size = 0.25, alpha = 0.25) +
      labs(x = NULL, y = input$ref_score) +
      theme_classic(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  })

  output$cox_atlas <- renderPlot({
    req(reference)
    x <- reference$tcga_cox
    x <- x[x$endpoint %in% c("OS", "DSS", "PFI") & x$score %in% score_cols, , drop = FALSE]
    x$log2HR <- pmax(pmin(log2(x$HR_per_SD), 1.2), -1.2)
    x$neglog10p <- pmin(-log10(x$p_cox), 8)
    top <- x |>
      dplyr::group_by(cancer) |>
      dplyr::summarise(best = max(neglog10p, na.rm = TRUE), .groups = "drop") |>
      dplyr::arrange(desc(best)) |>
      dplyr::slice_head(n = 18) |>
      dplyr::pull(cancer)
    x <- x[x$cancer %in% top, , drop = FALSE]
    x$cancer <- factor(x$cancer, levels = rev(top))
    ggplot(x, aes(score, cancer)) +
      geom_tile(aes(fill = log2HR), colour = "white", linewidth = 0.2) +
      geom_point(aes(size = neglog10p), shape = 21, fill = "white", colour = "black", stroke = 0.2) +
      facet_grid(. ~ endpoint) +
      scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#b2182b", midpoint = 0) +
      scale_size(range = c(0.5, 3.5)) +
      labs(x = NULL, y = NULL, fill = "log2 HR", size = "-log10 p") +
      theme_classic(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())
  })

  output$ref_table <- renderDT({
    req(ref_filtered())
    datatable(ref_filtered()[, c("patient", "cancer", score_cols), drop = FALSE], options = list(pageLength = 10, scrollX = TRUE))
  })

  clinical_data <- eventReactive(input$run_surv, {
    if (isTRUE(input$use_example_clin)) {
      mdri_read_clinical(file.path(ext_dir, "example_clinical.tsv"))
    } else {
      req(input$clin_file)
      path <- input$clin_file$datapath
      if (grepl("\\.csv$", input$clin_file$name, ignore.case = TRUE)) {
        tmp <- read.csv(path, check.names = FALSE)
        tf <- tempfile(fileext = ".tsv")
        write.table(tmp, tf, sep = "\t", quote = FALSE, row.names = FALSE)
        mdri_read_clinical(tf)
      } else {
        mdri_read_clinical(path)
      }
    }
  })

  surv_merged <- eventReactive(input$run_surv, {
    req(scores(), clinical_data())
    mdri_merge_clinical(scores(), clinical_data())
  })

  output$cox_table <- renderDT({
    req(surv_merged())
    if (!all(c(input$time_col, input$event_col, input$surv_score) %in% colnames(surv_merged()))) {
      return(datatable(data.frame(message = "Selected score/time/event columns were not found.")))
    }
    res <- mdri_cox(surv_merged(), input$surv_score, input$time_col, input$event_col)
    if (is.null(res)) return(datatable(data.frame(message = "Not enough samples/events for Cox model.")))
    datatable(res, options = list(dom = "t"))
  })

  output$km_plot <- renderPlot({
    req(surv_merged())
    if (!all(c(input$time_col, input$event_col, input$surv_score) %in% colnames(surv_merged()))) return(NULL)
    kd <- mdri_km_data(surv_merged(), input$surv_score, input$time_col, input$event_col)
    if (is.null(kd)) return(NULL)
    ggplot(kd$curve, aes(time, survival, color = group)) +
      geom_step(linewidth = 0.7) +
      scale_color_manual(values = c(High = "#d25774", Low = "#4490c4")) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = "Time", y = "Survival probability", color = "MDRi group",
           title = paste0(input$surv_score, " median split, log-rank p = ", signif(kd$p, 3))) +
      theme_classic(base_size = 12)
  })

  output$file_table <- renderDT({
    files <- list.files(root_dir, recursive = TRUE, full.names = FALSE)
    datatable(data.frame(file = files), options = list(pageLength = 15))
  })
}

shinyApp(ui, server)
