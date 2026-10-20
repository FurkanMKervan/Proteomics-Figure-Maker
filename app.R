# ==============================================================================
# ProteoPrism - Publication-Grade Figure Maker for Differential Proteomics
# High-End Visualization Suite for Modules 1, 2, and 3
# ==============================================================================

options(shiny.maxRequestSize = 2000 * 1024^2)

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(scales)
  library(grid)
  if (requireNamespace("ggrepel", quietly = TRUE)) library(ggrepel)
  if (requireNamespace("ggsci", quietly = TRUE)) library(ggsci)
  if (requireNamespace("pheatmap", quietly = TRUE)) library(pheatmap)
  if (requireNamespace("RColorBrewer", quietly = TRUE)) library(RColorBrewer)
  if (requireNamespace("viridis", quietly = TRUE)) library(viridis)
  if (requireNamespace("svglite", quietly = TRUE)) library(svglite)
})

# Locate and source core plotting engine
engine_path <- "figure_maker_engine.R"
if (!file.exists(engine_path)) {
  candidates <- c(
    file.path(getwd(), "figure_maker_engine.R"),
    file.path(getwd(), "Figure_Maker", "figure_maker_engine.R"),
    "c:/Users/serife-ayaz/Desktop/Differential_Proteomics_Pipeline/Figure_Maker/figure_maker_engine.R"
  )
  for (cand in candidates) {
    if (file.exists(cand)) { engine_path <- cand; break }
  }
}
if (file.exists(engine_path)) {
  source(engine_path)
}

DEFAULT_DATA_DIR <- "c:/Users/serife-ayaz/Desktop/OMV_Treatment_Sec"

app_theme <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  primary = "#2c3e50",
  secondary = "#18bc9c",
  success = "#18bc9c",
  danger = "#e74c3c",
  info = "#3498db",
  base_font = font_google("Inter"),
  heading_font = font_google("Inter")
)

# ------------------------------------------------------------------------------
# UI DEFINITION
# ------------------------------------------------------------------------------

ui <- page_navbar(
  title = span(style = "font-weight: 800; color: #ffffff; font-size: 1.25rem; letter-spacing: -0.5px;", "ProteoPrism"),
  theme = app_theme,
  bg = "#2c3e50",
  fillable = FALSE,

  # TAB 1: FIGURE MAKER STUDIO
  nav_panel(
    title = "📊 Figure Designer",
    layout_sidebar(
      sidebar = sidebar(
        width = 390,
        open = "always",

        # 1. DATA SOURCE CARD
        accordion(
          open = c("data_source_panel", "plot_type_panel", "styling_panel"),
          accordion_panel(
            "📁 1. Data Source & Import",
            value = "data_source_panel",
            radioButtons(
              "data_input_mode",
              "Select Data Source:",
              choices = c(
                "Pipeline Results Folder" = "folder",
                "Upload Custom File" = "upload",
                "Load Built-in Demo" = "demo"
              ),
              selected = if (dir.exists(DEFAULT_DATA_DIR)) "folder" else "demo"
            ),

            conditionalPanel(
              condition = "input.data_input_mode == 'folder'",
              textInput("folder_path", "Results Folder Path:", value = DEFAULT_DATA_DIR),
              actionButton("scan_folder_btn", "🔄 Scan / Refresh Files", class = "btn-sm btn-outline-primary mb-2"),
              selectInput("folder_file_select", "Select Output File:", choices = NULL),
              fileInput("folder_design_file_upload", "🏷️ Custom Design File (Optional, design.txt / .csv):", accept = c(".txt", ".tsv", ".csv"))
            ),

            conditionalPanel(
              condition = "input.data_input_mode == 'upload'",
              fileInput("uploaded_file", "Upload Expression / DEA File (.csv, .tsv, .xlsx):", accept = c(".csv", ".tsv", ".txt", ".xlsx")),
              fileInput("design_file_upload", "🏷️ Design / Grouping File (Optional, design.txt / .csv):", accept = c(".txt", ".tsv", ".csv"))
            ),

            conditionalPanel(
              condition = "input.data_input_mode == 'demo'",
              selectInput(
                "demo_type",
                "Select Demo Dataset:",
                choices = c(
                  "Enrichment: GO BP Terms" = "go_bp",
                  "GSEA: Hallmark Pathways (NES)" = "gsea",
                  "DEA: Volcano & Protein Abundance" = "dea",
                  "TF Footprinting: DoRothEA Activities" = "tf"
                )
              )
            ),

            uiOutput("data_badge_ui")
          ),

          # 2. PLOT SELECTION & PARAMETERS
          accordion_panel(
            "📈 2. Plot Type & Features",
            value = "plot_type_panel",
            selectInput("selected_plot_type", "Plot Type:", choices = NULL),

            # Sub-controls for Functional Plots
            conditionalPanel(
              condition = "input.selected_plot_type == 'dotplot' || input.selected_plot_type == 'barplot' || input.selected_plot_type == 'lollipop'",
              selectInput("enrich_x_metric", "X-Axis Metric:",
                          choices = c("Gene Ratio" = "GeneRatio", "Gene Count" = "Count", "Fold Enrichment" = "FoldEnrichment", "Rich Factor" = "RichFactor", "-log10(p.adjust)" = "negLog10Padj", "NES" = "NES"),
                          selected = "GeneRatio"),
              selectInput("enrich_color_metric", "Color Metric:",
                          choices = c("Adjusted p-value (FDR)" = "p.adjust", "Raw p-value" = "pvalue", "Fold Enrichment" = "FoldEnrichment", "Count" = "Count"),
                          selected = "p.adjust"),
              sliderInput("enrich_top_n", "Top N Terms to Display:", min = 5, max = 40, value = 15, step = 1),
              sliderInput("enrich_wrap", "Description Wrap Width:", min = 25, max = 80, value = 45, step = 5),

              # EXCLUDE / INCLUDE TERM FILTER
              div(style = "background-color: #f8f9fa; padding: 10px; border-radius: 6px; margin-top: 10px; border: 1px solid #e9ecef;",
                tags$b("🔍 Term Filtering:"),
                radioButtons(
                  "enrich_filter_mode",
                  label = NULL,
                  choices = c(
                    "❌ Exclude Selected Terms" = "exclude",
                    "🎯 Include Selected Terms Only" = "include"
                  ),
                  selected = "exclude"
                ),
                selectizeInput(
                  "enrich_selected_terms",
                  label = NULL,
                  choices = NULL,
                  multiple = TRUE,
                  options = list(placeholder = "Search and select terms to filter...")
                )
              ),

              # P-VALUE STABILITY & ANNOTATIONS
              div(style = "background-color: #f8f9fa; padding: 10px; border-radius: 6px; margin-top: 10px; border: 1px solid #e9ecef;",
                tags$b("⚖️ P-Value & Color Scale Calibration:"),
                selectInput(
                  "color_scale_limits_mode",
                  "Color Scale Limits (Color Stability):",
                  choices = c(
                    "Full Dataset Range (Fixed Limits)" = "full_data",
                    "Dynamic (Displayed Terms Only)" = "dynamic",
                    "Standard Significance (0 to 0.05)" = "standard_0_05"
                  ),
                  selected = "full_data"
                ),
                checkboxInput("enrich_show_p_text", "Display P-Value Values Beside Points / Bars", value = FALSE)
              )
            ),

            # Sub-controls for GSEA Diverging Lollipop
            conditionalPanel(
              condition = "input.selected_plot_type == 'gsea_diverging'",
              sliderInput("gsea_top_n", "Top N Enriched & Suppressed:", min = 5, max = 25, value = 10, step = 1),
              sliderInput("gsea_wrap", "Description Wrap Width:", min = 25, max = 80, value = 45, step = 5),

              div(style = "background-color: #f8f9fa; padding: 10px; border-radius: 6px; margin-top: 10px; border: 1px solid #e9ecef;",
                tags$b("🔍 Term Filtering:"),
                radioButtons(
                  "gsea_filter_mode",
                  label = NULL,
                  choices = c(
                    "❌ Exclude Selected Terms" = "exclude",
                    "🎯 Include Selected Terms Only" = "include"
                  ),
                  selected = "exclude"
                ),
                selectizeInput(
                  "gsea_selected_terms",
                  label = NULL,
                  choices = NULL,
                  multiple = TRUE,
                  options = list(placeholder = "Search and select terms to filter...")
                )
              ),

              div(style = "background-color: #f8f9fa; padding: 10px; border-radius: 6px; margin-top: 10px; border: 1px solid #e9ecef;",
                tags$b("🎨 Colors & Banner Headers:"),
                selectInput(
                  "gsea_color_by",
                  "Point Color (Color By):",
                  choices = c(
                    "FDR / p.adjust (Significance)" = "fdr",
                    "Raw p-value" = "pvalue",
                    "Regulation Status (Activated / Suppressed)" = "status"
                  ),
                  selected = "fdr"
                ),
                conditionalPanel(
                  condition = "input.gsea_color_by == 'fdr' || input.gsea_color_by == 'pvalue'",
                  selectInput(
                    "gsea_palette_choice",
                    "FDR / p-value Color Scale:",
                    choices = c(
                      "Red (Significant) - Light Grey (Non-Sig)" = "red_grey",
                      "Viridis (Purple - Green - Yellow)" = "viridis",
                      "Plasma (Dark Purple - Orange - Yellow)" = "plasma",
                      "Red - Navy Blue Gradient" = "red_blue",
                      "YlOrRd (Yellow - Orange - Red)" = "ylorrd"
                    ),
                    selected = "red_grey"
                  )
                ),
                checkboxInput("gsea_show_top_headers", "Show 'Activated' and 'Suppressed' Top Headers", value = TRUE),
                conditionalPanel(
                  condition = "input.gsea_show_top_headers == true",
                  radioButtons(
                    "gsea_header_style",
                    "Header Style:",
                    choices = c("Banner Box (◄ Suppressed | Activated ►)" = "banner", "Simple Text" = "text"),
                    selected = "banner",
                    inline = TRUE
                  )
                )
              ),

              div(style = "background-color: #f8f9fa; padding: 10px; border-radius: 6px; margin-top: 10px; border: 1px solid #e9ecef;",
                tags$b("🏷️ P-Value Labels (Lollipop):"),
                checkboxInput("gsea_show_p_labels", "Display P-Value Labels Beside Points", value = TRUE),
                conditionalPanel(
                  condition = "input.gsea_show_p_labels == true",
                  selectInput(
                    "gsea_p_label_type",
                    "P-Value Format:",
                    choices = c(
                      "FDR / p.adjust (e.g., FDR=0.012)" = "fdr",
                      "Raw p-value (e.g., p=0.003)" = "pval",
                      "Significance Stars (*, **, ***)" = "stars",
                      "Both (p and FDR)" = "both"
                    ),
                    selected = "fdr"
                  ),
                  sliderInput("gsea_p_label_size", "Label Font Size:", min = 2.5, max = 5.0, value = 3.5, step = 0.2)
                ),
                checkboxInput("gsea_size_by_set", "Scale Point Size by Gene Set Size (setSize)", value = TRUE)
              )
            ),

            # Sub-controls for Volcano Plot
            conditionalPanel(
              condition = "input.selected_plot_type == 'volcano'",
              radioButtons("volcano_p_type", "Significance Metric:", choices = c("Adjusted p-value / FDR" = "adj.P.Val", "Raw p-value" = "P.Value"), inline = TRUE),
              fluidRow(
                column(6, numericInput("volcano_fc_cutoff", "log2FC Cutoff:", value = 1.0, min = 0.1, max = 5.0, step = 0.1)),
                column(6, numericInput("volcano_p_cutoff", "P Cutoff:", value = 0.05, min = 0.0001, max = 0.2, step = 0.01))
              ),
              radioButtons("volcano_label_mode", "Gene Labeling Mode:",
                           choices = c(
                             "Top N by Max Fold Change (FC Max)" = "top_fc",
                             "Top N by Significance (P-value / FDR)" = "top_p",
                             "Top N Combined (|logFC| × -log10 P)" = "top_combined",
                             "Custom Selected Genes" = "custom",
                             "None (Unlabeled)" = "none"
                           ),
                           selected = "top_fc"),
              conditionalPanel(
                condition = "input.volcano_label_mode == 'top_fc' || input.volcano_label_mode == 'top_p' || input.volcano_label_mode == 'top_combined' || input.volcano_label_mode == 'top_n'",
                fluidRow(
                  column(7, sliderInput("volcano_top_labels", "Number of Top Labels (N):", min = 1, max = 30, value = 10, step = 1)),
                  column(5, radioButtons("volcano_balance_mode", "Label Distribution:",
                                         choices = c(
                                           "Both Sides (Up & Down)" = "per_side",
                                           "Overall Highest (|FC|)" = "overall"
                                         ),
                                         selected = "per_side"))
                ),
                checkboxInput("volcano_sig_only", "Label Significant Genes Only (Exceeding Cutoffs)", value = TRUE)
              ),
              conditionalPanel(
                condition = "input.volcano_label_mode == 'custom'",
                selectizeInput("volcano_custom_genes", "Genes to Label:", choices = NULL, multiple = TRUE, options = list(placeholder = "Search and pick genes..."))
              ),
              fluidRow(
                column(6, sliderInput("volcano_label_size", "Label Font Size:", min = 2.0, max = 6.0, value = 3.5, step = 0.2)),
                column(6, sliderInput("volcano_point_size", "Point Size:", min = 1, max = 5, value = 2.2, step = 0.2))
              ),
              sliderInput("volcano_alpha", "Point Opacity (Alpha):", min = 0.2, max = 1.0, value = 0.75, step = 0.05),
              checkboxInput("volcano_sym_x", "Force Symmetric X-Axis", value = TRUE),
              checkboxInput("volcano_show_lines", "Show Cutoff Threshold Lines", value = TRUE)
            ),

            # Sub-controls for Protein Abundance
            conditionalPanel(
              condition = "input.selected_plot_type == 'protein_box'",
              selectizeInput("target_gene_select", "Select Protein / Gene:", choices = NULL, options = list(placeholder = "Search gene name or protein ID...")),
              radioButtons("box_style", "Plot Style:",
                           choices = c(
                             "Profile Plot (Sample-Level Lines)" = "profile_samples",
                             "Profile Plot (Group Mean ± SEM)" = "profile_group",
                             "Boxplot + Jitter" = "boxplot"
                           ),
                           selected = "profile_samples",
                           inline = FALSE),
              conditionalPanel(
                condition = "input.box_style == 'profile_group' || input.box_style == 'profile_samples'",
                fluidRow(
                  column(6, checkboxInput("profile_show_line", "Connect Means with Line", value = TRUE)),
                  column(6, checkboxInput("profile_show_points", "Show Individual Sample Points (Jitter)", value = TRUE))
                )
              )
            ),

            # Sub-controls for Expression Heatmap
            conditionalPanel(
              condition = "input.selected_plot_type == 'heatmap'",
              radioButtons(
                "heatmap_scale",
                "Data Scaling:",
                choices = c("Row Z-Score (Standard)" = "row", "Raw Log2 Abundances" = "none"),
                selected = "row",
                inline = TRUE
              ),
              fluidRow(
                column(6, checkboxInput("heatmap_cluster_rows", "Cluster Rows (Proteins)", value = TRUE)),
                column(6, checkboxInput("heatmap_cluster_cols", "Cluster Columns (Samples)", value = FALSE))
              ),
              fluidRow(
                column(6, checkboxInput("heatmap_show_rownames", "Show Gene / Protein Names", value = TRUE)),
                column(6, checkboxInput("heatmap_show_colnames", "Show Sample Names", value = TRUE))
              ),
              selectInput(
                "heatmap_color_palette",
                "Heatmap Color Palette:",
                choices = c(
                  "Blue - White - Red (RdBu / NPG)" = "RdBu",
                  "Blue - Red (BWR)" = "bwr",
                  "Viridis" = "viridis",
                  "Plasma" = "plasma",
                  "Magma" = "magma"
                ),
                selected = "RdBu"
              ),
              sliderInput("heatmap_top_n", "Maximum Proteins to Display (0 = All):", min = 0, max = 200, value = 50, step = 10)
            ),

            # Sub-controls for Sample Grouping & Manual Naming (Shared by Profile Plot & Heatmap)
            conditionalPanel(
              condition = "input.selected_plot_type == 'protein_box' || input.selected_plot_type == 'heatmap'",
              div(
                style = "background-color: #f8f9fa; padding: 12px; border-radius: 6px; margin-top: 15px; border: 1px solid #dee2e6;",
                div(
                  class = "d-flex align-items-center justify-content-between mb-2",
                  tags$b("👥 Sample Grouping:"),
                  uiOutput("design_status_badge")
                ),
                radioButtons(
                  "group_mode_select",
                  label = NULL,
                  choices = c(
                    "Automatic / Design File" = "auto",
                    "Manual Group Editing" = "manual"
                  ),
                  selected = "auto",
                  inline = TRUE
                ),
                conditionalPanel(
                  condition = "input.group_mode_select == 'manual'",
                  helpText("Assign custom biological names to detected groups:"),
                  uiOutput("manual_group_inputs_ui")
                )
              )
            ),

            # Sub-controls for TF Plots
            conditionalPanel(
              condition = "input.selected_plot_type == 'tf_bar' || input.selected_plot_type == 'tf_lollipop'",
              sliderInput("tf_top_n", "Top N TFs per Regulation Group:", min = 5, max = 30, value = 15, step = 1),
              radioButtons("tf_orientation", "Orientation:", choices = c("Horizontal" = "horizontal", "Vertical" = "vertical"), inline = TRUE),
              div(style = "background-color: #f8f9fa; padding: 10px; border-radius: 6px; margin-top: 10px; border: 1px solid #e9ecef;",
                tags$b("🔍 TF Filtering:"),
                radioButtons(
                  "tf_filter_mode",
                  label = NULL,
                  choices = c(
                    "❌ Exclude Selected TFs" = "exclude",
                    "🎯 Include Selected TFs Only" = "include"
                  ),
                  selected = "exclude"
                ),
                selectizeInput("tf_selected_tfs", label = NULL, choices = NULL, multiple = TRUE, options = list(placeholder = "Search and select TFs to filter..."))
              )
            )
          ),

          # 3. STYLING, PALETTES & TITLES
          accordion_panel(
            "🎨 3. Aesthetics, Palettes & Titles",
            value = "styling_panel",
            selectInput(
              "journal_palette",
              "Color Theme / Journal Palette:",
              choices = c(
                "Nature Publishing Group (NPG)" = "NPG",
                "New England J of Medicine (NEJM)" = "NEJM",
                "The Lancet" = "Lancet",
                "J of Clinical Oncology (JCO)" = "JCO",
                "Science (AAAS)" = "AAAS",
                "ColorBrewer Set1" = "Set1",
                "ColorBrewer Set2" = "Set2",
                "ColorBrewer Dark2" = "Dark2",
                "Viridis (Perceptually Uniform)" = "viridis",
                "Custom Colors" = "custom"
              ),
              selected = "NPG"
            ),

            conditionalPanel(
              condition = "input.journal_palette == 'custom' || input.selected_plot_type == 'volcano' || input.selected_plot_type == 'tf_bar' || input.selected_plot_type == 'gsea_diverging'",
              fluidRow(
                column(4, textInput("col_up", "Up / Act Color:", value = "#E64B35")),
                column(4, textInput("col_down", "Down / Rep Color:", value = "#3C5488")),
                column(4, textInput("col_ns", "NS / Low Color:", value = "#B09C85"))
              )
            ),

            selectInput(
              "plot_theme",
              "Publication Theme:",
              choices = c(
                "Nature Classic (theme_classic)" = "classic",
                "Clean Minimal (theme_minimal)" = "minimal",
                "Black & White (theme_bw)" = "bw",
                "GraphPad Prism Style" = "prism"
              ),
              selected = "classic"
            ),
            sliderInput("base_font_size", "Base Font Size (pt):", min = 9, max = 18, value = 12, step = 1),

            textInput("custom_title", "Plot Title:", value = ""),
            textInput("custom_subtitle", "Subtitle:", value = ""),
            textInput("custom_xlab", "X-Axis Label (Leave blank for default):", value = ""),
            textInput("custom_ylab", "Y-Axis Label (Leave blank for default):", value = "")
          )
        )
      ),

      # MAIN PLOT DISPLAY & EXPORT AREA
      card(
        card_header(
          class = "d-flex justify-content-between align-items-center",
          span("🖼️ High-Resolution Figure Preview", style = "font-weight: 700;"),
          span(uiOutput("plot_info_badge", inline = TRUE))
        ),
        card_body(
          plotOutput("main_plot", height = "580px")
        )
      ),

      # EXPORT CONTROL CARD
      card(
        card_header("💾 Export Publication Figure", style = "font-weight: 700;"),
        card_body(
          fluidRow(
            column(3,
              selectInput(
                "export_format",
                "Export Format:",
                choices = c(
                  "PNG (300 DPI - Standard)" = "png_300",
                  "PNG (600 DPI - High-Res)" = "png_600",
                  "PDF (Vector - Journals)" = "pdf",
                  "SVG (Vector - Illustrator)" = "svg",
                  "TIFF (300 DPI - Lossless)" = "tiff"
                ),
                selected = "pdf"
              )
            ),
            column(3,
              selectInput(
                "export_preset",
                "Journal Dimension Preset:",
                choices = c(
                  "Nature 1-Column (89 mm / 3.5 in)" = "nature_1col",
                  "Nature 1.5-Column (120 mm / 4.7 in)" = "nature_1.5col",
                  "Nature 2-Column (183 mm / 7.2 in)" = "nature_2col",
                  "Standard Publication (8 x 6 in)" = "standard",
                  "Square Figure (7 x 7 in)" = "square",
                  "Custom Size" = "custom"
                ),
                selected = "standard"
              )
            ),
            column(2, numericInput("export_width", "Width (in):", value = 8, min = 2, max = 20, step = 0.5)),
            column(2, numericInput("export_height", "Height (in):", value = 6, min = 2, max = 20, step = 0.5)),
            column(2,
              div(style = "margin-top: 28px;",
                downloadButton("download_fig_btn", "⬇️ Download Figure", class = "btn-success w-100")
              )
            )
          )
        )
      ),

      # DATA TABLE CARD
      card(
        card_header(
          class = "d-flex justify-content-between align-items-center",
          span("📋 Filtered Data Table (Currently Plotted Data)", style = "font-weight: 700;"),
          span(style = "font-size: 0.82rem; color: #6c757d;", "Values in this table reflect exact numeric file data without modification")
        ),
        card_body(
          uiOutput("table_summary_header"),
          DT::dataTableOutput("plot_data_table")
        )
      )
    )
  ),

  # TAB 2: GUIDELINES & JOURNAL SPECS
  nav_panel(
    title = "📖 Journal Specs & Guide",
    card(
      card_header("Scientific Journal Figure Preparation Guidelines"),
      card_body(
        markdown("
### 🏛️ Publication Standards for High-Impact Journals (Nature, Cell, Science)

1. **Resolution & Format:**
   - **Vector format (PDF / SVG):** Always preferred for line drawings, bar plots, dot plots, and volcano plots. Fonts stay perfectly sharp regardless of zoom level.
   - **Bitmaps (PNG / TIFF):** Use **300 DPI** minimum for color figures, **600 DPI** for figures with fine text/details.

2. **Dimensions & Columns:**
   - **1-Column:** 89 mm width (~3.5 inches). Ideal for single volcano plots, TF barplots, or boxplots.
   - **1.5-Column:** 120 mm width (~4.7 inches). Great for medium dot plots.
   - **2-Column:** 183 mm width (~7.2 inches). Recommended for comprehensive pathway dotplots with long descriptions and full figure panels.

3. **Color Palettes:**
   - Use perceptually uniform palettes (Viridis, Plasma) or respected journal palettes (NPG, Lancet, NEJM) that maintain legibility for readers with color vision deficiencies.
        ")
      )
    )
  )
)

# ------------------------------------------------------------------------------
# SERVER LOGIC
# ------------------------------------------------------------------------------

server <- function(input, output, session) {

  raw_data <- reactiveVal(NULL)
  data_type <- reactiveVal("unknown") # "dea", "enrich", "tf", "gsea"

  # Scan folder on startup or button click
  observe({
    folder <- input$folder_path
    if (!is.null(folder) && dir.exists(folder)) {
      csv_files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
      if (length(csv_files) > 0) {
        file_names <- basename(csv_files)
        names(csv_files) <- file_names
        updateSelectInput(session, "folder_file_select", choices = csv_files, selected = csv_files[1])
      }
    }
  }) %>% bindEvent(input$scan_folder_btn, ignoreInit = FALSE)

  # Load data when source changes
  observe({
    mode <- input$data_input_mode
    df <- NULL

    if (mode == "folder") {
      selected_path <- input$folder_file_select
      if (!is.null(selected_path) && file.exists(selected_path)) {
        df <- tryCatch(read.csv(selected_path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
      }
    } else if (mode == "upload") {
      req(input$uploaded_file)
      ext <- tools::file_ext(input$uploaded_file$name)
      path <- input$uploaded_file$datapath
      if (ext %in% c("csv", "txt")) {
        df <- tryCatch(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
      } else if (ext == "tsv") {
        df <- tryCatch(read.delim(path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
      } else if (ext == "xlsx" && requireNamespace("readxl", quietly = TRUE)) {
        df <- tryCatch(as.data.frame(readxl::read_excel(path)), error = function(e) NULL)
      }
    } else if (mode == "demo") {
      dtype <- input$demo_type
      if (dtype == "go_bp") {
        df <- data.frame(
          ID = paste0("GO:000", 10:24),
          Description = c(
            "innate immune response", "defense response", "immune system process",
            "response to external stimulus", "leukocyte homeostasis", "proteolysis",
            "cellular response to cytokine", "antigen processing and presentation",
            "neutrophil activation", "platelet degranulation", "extracellular matrix organization",
            "cell-cell adhesion", "apoptotic process", "inflammatory response", "interferon-gamma signaling"
          ),
          GeneRatio = c("19/77", "25/77", "30/77", "29/77", "6/77", "16/77", "14/77", "11/77", "18/77", "15/77", "22/77", "13/77", "17/77", "20/77", "12/77"),
          Count = c(19, 25, 30, 29, 6, 16, 14, 11, 18, 15, 22, 13, 17, 20, 12),
          FoldEnrichment = c(2.33, 2.00, 1.70, 1.71, 4.28, 2.19, 2.45, 3.12, 2.65, 2.90, 1.85, 2.10, 1.95, 2.40, 3.05),
          pvalue = c(8.6e-5, 9.5e-5, 4.7e-4, 5.5e-4, 8.3e-4, 0.0011, 0.0019, 0.0028, 0.0035, 0.0049, 0.0062, 0.011, 0.015, 0.021, 0.032),
          p.adjust = c(0.0001, 0.00015, 0.0004, 0.0005, 0.0008, 0.0012, 0.002, 0.0035, 0.004, 0.0055, 0.007, 0.012, 0.018, 0.025, 0.038),
          stringsAsFactors = FALSE
        )
      } else if (dtype == "gsea") {
        df <- data.frame(
          ID = paste0("HALLMARK_", 1:12),
          Description = c(
            "Interferon Gamma Response", "IL6 JAK STAT3 Signaling", "TNFA Signaling via NFKB",
            "Inflammatory Response", "Complement", "Coagulation",
            "Protein Secretion", "Apical Surface", "P53 Pathway", "E2F Targets", "G2M Checkpoint", "MYC Targets V1"
          ),
          NES = c(2.14, 1.85, 1.72, 1.58, 1.45, 1.36, -1.53, -1.49, -1.39, -1.35, -1.62, -1.78),
          setSize = c(19, 15, 22, 28, 14, 27, 12, 10, 18, 24, 20, 16),
          pvalue = c(0.0003, 0.0012, 0.0025, 0.0051, 0.012, 0.018, 0.015, 0.022, 0.045, 0.052, 0.008, 0.002),
          p.adjust = c(0.016, 0.025, 0.032, 0.048, 0.065, 0.082, 0.042, 0.051, 0.095, 0.110, 0.035, 0.022),
          stringsAsFactors = FALSE
        )
      } else if (dtype == "dea") {
        genes <- c("HSP90B1", "EEF1A1", "TPM3", "PSMB3", "TPD52L2", "GAPDH", "ACTB", "C3", "IL6", "CASP3", "MMP3", "PTX3", "GSN", "CFH", "CCL2", "HSPD1", paste0("PROT_", 1:40))
        fc <- c(0.72, -1.35, -0.31, 2.67, -3.42, 0.1, -0.05, 1.85, 3.1, -2.1, 2.4, 1.9, 0.8, -1.2, 2.2, 0.9, seq(-1.5, 1.5, length.out = 40))
        pval <- c(0.33, 0.0002, 0.15, 0.001, 0.00001, 0.8, 0.9, 0.0005, 0.0001, 0.002, 0.0008, 0.001, 0.2, 0.01, 0.0015, 0.1, seq(0.01, 0.8, length.out = 40))
        padj <- pmin(pval * 1.5, 1.0)

        s1 <- 20 + fc/2 + 0.1
        s2 <- 20 + fc/2 - 0.1
        c1 <- 20 - fc/2 + 0.1
        c2 <- 20 - fc/2 - 0.1

        df <- data.frame(
          Gene = genes,
          logFC = fc,
          P.Value = pval,
          adj.P.Val = padj,
          Young_OMV_1 = s1,
          Young_OMV_2 = s2,
          Young_Control_1 = c1,
          Young_Control_2 = c2,
          stringsAsFactors = FALSE
        )
      } else if (dtype == "tf") {
        df <- data.frame(
          source = c("RELA", "ETS1", "FOS", "JUN", "NFKB1", "NFYB", "STAT3", "E2F7", "HIF1A", "MAZ", "SP1", "MYC", "TP53", "STAT1"),
          score = c(4.25, 3.74, 3.63, 2.91, 2.41, 2.34, 1.95, -2.61, -2.38, 2.25, -1.85, -2.10, -2.55, 3.12),
          p_value = c(0.00002, 0.0002, 0.0003, 0.0037, 0.016, 0.019, 0.045, 0.009, 0.017, 0.024, 0.065, 0.035, 0.012, 0.001),
          Status = c(rep("Activated", 7), rep("Repressed", 7)),
          stringsAsFactors = FALSE
        )
      }
    }

    if (!is.null(df) && is.data.frame(df) && nrow(df) > 0) {
      raw_data(df)
      cols_lower <- tolower(colnames(df))
      if (any(c("description", "generatio", "enrichment", "richfactor") %in% cols_lower)) {
        if ("nes" %in% cols_lower) {
          data_type("gsea")
        } else {
          data_type("enrich")
        }
      } else if (any(c("source", "tf") %in% cols_lower) && any(c("score", "activity") %in% cols_lower)) {
        data_type("tf")
      } else if ((any(c("logfc", "log2fc", "fold_change") %in% cols_lower) && any(c("p.value", "adj.p.val", "padj", "pval") %in% cols_lower)) ||
                 (any(c("gene.names", "protein.ids", "gene", "protein", "genes") %in% cols_lower) && sum(sapply(df, is.numeric)) >= 2)) {
        data_type("dea")
      } else {
        data_type("generic")
      }
    }
  })

  # --- 1.1 SAMPLE COLUMNS & EXPERIMENTAL DESIGN ---
  current_sample_cols <- reactive({
    df <- raw_data()
    if (is.null(df)) return(character(0))
    anno_pat <- "^(gene|count|row_id|logfc|log2fc|fc|aveexpr|t|b|p[._]?val.*|adj[._]?p.*|padj|fdr|qval.*|sca[._].*|disp.*|prot.*|sig.*|comp.*)$"
    skip_idx <- grepl(anno_pat, colnames(df), ignore.case = TRUE)
    names(df)[sapply(df, is.numeric) & !skip_idx]
  })

  current_design_df <- reactive({
    mode <- input$data_input_mode
    if (mode == "folder") {
      # 1. Custom uploaded file takes priority if provided
      if (!is.null(input$folder_design_file_upload)) {
        return(parse_design_table(input$folder_design_file_upload$datapath))
      }
      # 2. Check if design.txt exists in the active folder or parent folder
      selected_path <- input$folder_file_select
      search_dirs <- c()
      if (!is.null(selected_path) && file.exists(selected_path)) {
        search_dirs <- c(search_dirs, dirname(selected_path))
      }
      if (!is.null(input$folder_path) && dir.exists(input$folder_path)) {
        search_dirs <- c(search_dirs, input$folder_path)
      }
      search_dirs <- unique(search_dirs)
      for (sdir in search_dirs) {
        for (cand in c("design.txt", "design.csv", "Design.txt", "Design.csv")) {
          fpath <- file.path(sdir, cand)
          if (file.exists(fpath)) {
            parsed <- parse_design_table(fpath)
            if (!is.null(parsed)) return(parsed)
          }
        }
      }
    } else if (mode == "upload") {
      if (!is.null(input$design_file_upload)) {
        return(parse_design_table(input$design_file_upload$datapath))
      }
    }
    return(NULL)
  })

  detected_sample_groups <- reactive({
    s_cols <- current_sample_cols()
    if (length(s_cols) == 0) return(character(0))
    des <- current_design_df()
    clean_name <- function(x) gsub("^[X._]+", "", x)

    if (!is.null(des)) {
      matched <- sapply(clean_name(s_cols), function(s) {
        idx <- which(des$Sample == s)
        if (length(idx) > 0) des$Group[idx[1]] else s
      })
      return(setNames(matched, s_cols))
    }

    # Auto prefix detection: remove replicate suffix (_1, _2, .1, etc.)
    clean_grp <- sub("(_r?|\\.)[0-9]+$", "", s_cols)
    clean_grp <- gsub("^X[._]?", "", clean_grp)
    if (all(clean_grp == "") || length(unique(clean_grp)) == 1) {
      return(setNames(s_cols, s_cols))
    }
    if (all(grepl("^[0-9]+$", clean_grp))) {
      clean_grp <- paste("Group", clean_grp)
    }
    return(setNames(clean_grp, s_cols))
  })

  output$design_status_badge <- renderUI({
    s_cols <- current_sample_cols()
    if (length(s_cols) == 0) return(NULL)
    des <- current_design_df()
    grps <- detected_sample_groups()
    n_grps <- length(unique(grps))
    if (!is.null(des)) {
      span(class = "badge bg-success", title = "Automatically mapped from design file", paste0("✓ design.txt (", n_grps, " Groups)"))
    } else {
      span(class = "badge bg-info", title = "Inferred from sample name prefixes", paste0("ℹ Prefix (", n_grps, " Groups)"))
    }
  })

  output$manual_group_inputs_ui <- renderUI({
    s_cols <- current_sample_cols()
    if (length(s_cols) == 0) return(NULL)
    grps <- detected_sample_groups()
    u_grps <- unique(grps)

    inputs <- lapply(seq_along(u_grps), function(i) {
      grp_id <- u_grps[i]
      members <- names(grps)[grps == grp_id]
      members_str <- paste(members, collapse = ", ")
      textInput(
        inputId = paste0("custom_grp_name_", i),
        label = paste0("Group ", i, " (Samples: ", members_str, "):"),
        value = grp_id
      )
    })

    tagList(
      inputs,
      textAreaInput(
        "manual_group_csv",
        "Or Enter All Group Names in Order (Comma-separated):",
        value = paste(grps, collapse = ", "),
        rows = 2
      )
    )
  })

  active_group_mapping <- reactive({
    s_cols <- current_sample_cols()
    if (length(s_cols) == 0) return(NULL)

    if (isTRUE(input$group_mode_select == "manual")) {
      grps <- detected_sample_groups()
      u_grps <- unique(grps)

      # 1. Try individual group text inputs
      manual_map <- character(length(s_cols))
      for (i in seq_along(u_grps)) {
        inp_val <- input[[paste0("custom_grp_name_", i)]]
        val_use <- if (!is.null(inp_val) && nzchar(trimws(inp_val))) trimws(inp_val) else u_grps[i]
        manual_map[grps == u_grps[i]] <- val_use
      }

      if (all(nzchar(manual_map))) {
        return(setNames(manual_map, s_cols))
      }

      # 2. Fallback to CSV textarea if provided
      csv_val <- input$manual_group_csv
      if (!is.null(csv_val) && nzchar(trimws(csv_val))) {
        parts <- trimws(strsplit(csv_val, "[,;\n]+")[[1]])
        if (length(parts) == length(s_cols)) {
          return(setNames(parts, s_cols))
        }
      }
    }

    return(detected_sample_groups())
  })

  # --- 2. UPDATE PLOT CHOICES DYNAMICALLY ---
  observe({
    dtype <- data_type()
    req(dtype)

    choices <- list()
    if (dtype %in% c("enrich", "gsea")) {
      choices <- list(
        "Dot Plot (Enrichment)" = "dotplot",
        "Barplot (Enrichment)" = "barplot",
        "Lollipop Chart (Enrichment)" = "lollipop"
      )
      if (dtype == "gsea" || "NES" %in% colnames(raw_data())) {
        choices[["GSEA Diverging NES Lollipop Chart"]] <- "gsea_diverging"
      }
    } else if (dtype == "dea") {
      choices <- list(
        "Enhanced Volcano Plot" = "volcano",
        "Protein Abundance (Profile Plot / Boxplot)" = "protein_box",
        "Expression Heatmap (Pheatmap / Z-Score)" = "heatmap"
      )
    } else if (dtype == "tf") {
      choices <- list(
        "TF Activity Diverging Barplot" = "tf_bar",
        "TF Activity Lollipop Plot" = "tf_lollipop"
      )
    } else {
      choices <- list("Custom Dot Plot" = "dotplot", "Custom Barplot" = "barplot")
    }

    default_sel <- if (dtype == "dea") {
      if (!is.null(input$data_file) && grepl("absolute", input$data_file$name, ignore.case = TRUE)) "heatmap" else "volcano"
    } else if (dtype == "gsea") "gsea_diverging"
    else if (dtype == "enrich") "dotplot"
    else if (dtype == "tf") "tf_bar"
    else choices[[1]]

    updateSelectInput(session, "selected_plot_type", choices = choices, selected = default_sel)

    df <- raw_data()
    if (!is.null(df)) {
      if (dtype %in% c("enrich", "gsea") && "Description" %in% colnames(df)) {
        terms <- unique(df$Description)
        updateSelectizeInput(session, "enrich_selected_terms", choices = terms, server = TRUE)
        updateSelectizeInput(session, "gsea_selected_terms", choices = terms, server = TRUE)
      } else if (dtype == "dea") {
        std_dea <- standardize_dea_df(df)
        cols_to_keep <- intersect(c("Gene", "Protein_ID", "Display_Name"), colnames(std_dea))
        unique_entries <- unique(std_dea[, cols_to_keep, drop = FALSE])

        # Build choices: clean Display_Name as label and value
        gene_choices <- setNames(as.character(unique_entries$Display_Name), as.character(unique_entries$Display_Name))
        gene_choices <- gene_choices[!is.na(gene_choices) & gene_choices != ""]

        updateSelectizeInput(session, "volcano_custom_genes", choices = gene_choices, server = TRUE)
        updateSelectizeInput(session, "target_gene_select", choices = gene_choices, selected = gene_choices[1], server = TRUE)
      } else if (dtype == "tf") {
        tf_col <- grep("^(source|tf)$", colnames(df), ignore.case = TRUE, value = TRUE)
        if (length(tf_col) > 0) {
          tfs <- unique(na.omit(df[[tf_col[1]]]))
          updateSelectizeInput(session, "tf_selected_tfs", choices = tfs, server = TRUE)
        }
      }
    }
  })

  observe({
    preset <- input$export_preset
    req(preset)
    if (preset == "nature_1col") {
      updateNumericInput(session, "export_width", value = 3.5)
      updateNumericInput(session, "export_height", value = 4.0)
    } else if (preset == "nature_1.5col") {
      updateNumericInput(session, "export_width", value = 4.7)
      updateNumericInput(session, "export_height", value = 5.0)
    } else if (preset == "nature_2col") {
      updateNumericInput(session, "export_width", value = 7.2)
      updateNumericInput(session, "export_height", value = 6.0)
    } else if (preset == "standard") {
      updateNumericInput(session, "export_width", value = 8.0)
      updateNumericInput(session, "export_height", value = 6.0)
    } else if (preset == "square") {
      updateNumericInput(session, "export_width", value = 7.0)
      updateNumericInput(session, "export_height", value = 7.0)
    }
  })

  # --- 3. UI BADGES ---
  output$data_badge_ui <- renderUI({
    dtype <- data_type()
    badge_label <- switch(dtype,
      "enrich" = "Module 2: Functional / Pathway Enrichment",
      "gsea" = "Module 2: GSEA Results (NES)",
      "dea" = "Module 1: Differential Expression / Quant",
      "tf" = "Module 3: DoRothEA TF Footprinting",
      "Generic Data"
    )
    badge_color <- switch(dtype,
      "enrich" = "bg-primary",
      "gsea" = "bg-info",
      "dea" = "bg-danger",
      "tf" = "bg-success",
      "bg-secondary"
    )

    div(
      class = "mt-2 p-2 rounded text-white text-center",
      class = badge_color,
      style = "font-size: 0.82rem; font-weight: 600;",
      paste("Detected:", badge_label)
    )
  })

  output$plot_info_badge <- renderUI({
    p_type <- input$selected_plot_type
    req(p_type)
    span(
      class = "badge bg-light text-dark border",
      style = "font-size: 0.85rem;",
      paste("Active:", toupper(p_type))
    )
  })

  # --- 4. REACTIVE GGPLOT OBJECT GENERATOR ---
  current_plot <- reactive({
    df <- raw_data()
    req(df)
    dtype <- data_type()
    req(dtype)
    ptype <- input$selected_plot_type
    req(ptype)

    # Validation: Prevent race conditions when switching datasets
    valid_plot_map <- list(
      "dea" = c("volcano", "protein_box", "heatmap"),
      "enrich" = c("dotplot", "barplot", "lollipop", "gsea_diverging"),
      "gsea" = c("dotplot", "barplot", "lollipop", "gsea_diverging"),
      "tf" = c("tf_bar", "tf_lollipop"),
      "generic" = c("dotplot", "barplot", "volcano")
    )
    expected_plots <- valid_plot_map[[dtype]]
    if (!is.null(expected_plots) && !ptype %in% expected_plots) {
      return(NULL)
    }

    user_title <- if (nzchar(input$custom_title)) input$custom_title else NULL
    user_sub <- if (nzchar(input$custom_subtitle)) input$custom_subtitle else NULL
    user_xlab <- if (nzchar(input$custom_xlab)) input$custom_xlab else NULL
    user_ylab <- if (nzchar(input$custom_ylab)) input$custom_ylab else NULL

    palette_choice <- input$journal_palette
    up_c <- input$col_up
    down_c <- input$col_down
    ns_c <- input$col_ns
    th_name <- input$plot_theme
    b_size <- input$base_font_size

    # Calculate locked color limits if requested
    calc_p_limits <- function(col_name) {
      mode <- input$color_scale_limits_mode
      if (is.null(mode) || mode == "dynamic") return(NULL)
      if (mode == "standard_0_05") return(c(0, 0.05))
      if (mode == "full_data" && col_name %in% colnames(df)) {
        vals <- suppressWarnings(as.numeric(df[[col_name]]))
        vals <- vals[!is.na(vals) & vals >= 0]
        if (length(vals) > 0) return(c(min(vals), max(vals)))
      }
      return(NULL)
    }

    p <- NULL

    if (ptype == "dotplot") {
      p_lims <- calc_p_limits(input$enrich_color_metric)
      p <- plot_enrichment_dot(
        df = df,
        x_metric = input$enrich_x_metric,
        color_metric = input$enrich_color_metric,
        top_n = input$enrich_top_n,
        selected_terms = input$enrich_selected_terms,
        filter_mode = input$enrich_filter_mode,
        palette_type = if (palette_choice == "viridis") "viridis" else "gradient",
        gradient_low = up_c,
        gradient_high = down_c,
        wrap_length = input$enrich_wrap,
        p_limits = p_lims,
        show_p_labels = input$enrich_show_p_text,
        title = if (!is.null(user_title)) user_title else "Functional Pathway Enrichment",
        subtitle = user_sub,
        xlab = user_xlab,
        ylab = if (!is.null(user_ylab)) user_ylab else "",
        theme_name = th_name,
        base_size = b_size
      )
    } else if (ptype == "barplot") {
      p_lims <- calc_p_limits(input$enrich_color_metric)
      p <- plot_enrichment_bar(
        df = df,
        x_metric = input$enrich_x_metric,
        color_metric = input$enrich_color_metric,
        top_n = input$enrich_top_n,
        selected_terms = input$enrich_selected_terms,
        filter_mode = input$enrich_filter_mode,
        palette_type = if (palette_choice == "viridis") "viridis" else "gradient",
        gradient_low = up_c,
        gradient_high = down_c,
        wrap_length = input$enrich_wrap,
        p_limits = p_lims,
        show_p_labels = input$enrich_show_p_text,
        title = if (!is.null(user_title)) user_title else "Pathway Enrichment Barplot",
        subtitle = user_sub,
        xlab = user_xlab,
        ylab = if (!is.null(user_ylab)) user_ylab else "",
        theme_name = th_name,
        base_size = b_size
      )
    } else if (ptype == "lollipop") {
      p_lims <- calc_p_limits(input$enrich_color_metric)
      p <- plot_enrichment_lollipop(
        df = df,
        x_metric = input$enrich_x_metric,
        color_metric = input$enrich_color_metric,
        top_n = input$enrich_top_n,
        selected_terms = input$enrich_selected_terms,
        filter_mode = input$enrich_filter_mode,
        palette_type = if (palette_choice == "viridis") "viridis" else "gradient",
        gradient_low = up_c,
        gradient_high = down_c,
        wrap_length = input$enrich_wrap,
        p_limits = p_lims,
        show_p_labels = input$enrich_show_p_text,
        title = if (!is.null(user_title)) user_title else "Pathway Enrichment Lollipop",
        subtitle = user_sub,
        xlab = user_xlab,
        ylab = if (!is.null(user_ylab)) user_ylab else "",
        theme_name = th_name,
        base_size = b_size
      )
    } else if (ptype == "gsea_diverging") {
      gsea_c_metric <- if (!is.null(input$gsea_color_by) && input$gsea_color_by == "pvalue") "pvalue" else "p.adjust"
      p_lims <- calc_p_limits(gsea_c_metric)

      # Determine GSEA color palette
      gsea_pal_choice <- if (!is.null(input$gsea_palette_choice)) input$gsea_palette_choice else "red_grey"
      gsea_pal_type <- "gradient"
      gsea_vir_opt <- "D"
      gsea_low <- "#DC0000"
      gsea_high <- "#C0C0C0"

      if (gsea_pal_choice == "red_grey") {
        gsea_pal_type <- "gradient"
        gsea_low <- "#DC0000"  # Vibrant Red for significant
        gsea_high <- "#C0C0C0" # Neutral Light Grey for non-significant
      } else if (gsea_pal_choice == "viridis") {
        gsea_pal_type <- "viridis"
        gsea_vir_opt <- "D"
      } else if (gsea_pal_choice == "plasma") {
        gsea_pal_type <- "viridis"
        gsea_vir_opt <- "C"
      } else if (gsea_pal_choice == "red_blue") {
        gsea_pal_type <- "gradient"
        gsea_low <- up_c
        gsea_high <- down_c
      } else if (gsea_pal_choice == "ylorrd") {
        gsea_pal_type <- "brewer"
      }

      p <- plot_gsea_diverging(
        df = df,
        top_n = input$gsea_top_n,
        selected_terms = input$gsea_selected_terms,
        filter_mode = input$gsea_filter_mode,
        color_by = if (!is.null(input$gsea_color_by)) input$gsea_color_by else "fdr",
        palette_type = gsea_pal_type,
        gradient_low = gsea_low,
        gradient_high = gsea_high,
        viridis_opt = gsea_vir_opt,
        p_limits = p_lims,
        up_color = up_c,
        down_color = down_c,
        size_by_set = input$gsea_size_by_set,
        show_p_labels = input$gsea_show_p_labels,
        p_label_type = input$gsea_p_label_type,
        p_label_size = input$gsea_p_label_size,
        show_top_headers = if (!is.null(input$gsea_show_top_headers)) input$gsea_show_top_headers else TRUE,
        header_style = if (!is.null(input$gsea_header_style)) input$gsea_header_style else "banner",
        wrap_length = input$gsea_wrap,
        title = if (!is.null(user_title)) user_title else "GSEA Normalized Enrichment Scores (NES)",
        subtitle = user_sub,
        theme_name = th_name,
        base_size = b_size
      )
    } else if (ptype == "volcano") {
      p <- plot_volcano(
        df = df,
        p_val_type = input$volcano_p_type,
        fc_cutoff = input$volcano_fc_cutoff,
        p_cutoff = input$volcano_p_cutoff,
        up_color = up_c,
        down_color = down_c,
        ns_color = ns_c,
        point_size = input$volcano_point_size,
        point_alpha = input$volcano_alpha,
        top_n_labels = if (!is.null(input$volcano_top_labels)) input$volcano_top_labels else 10,
        selected_genes = input$volcano_custom_genes,
        label_type = if (!is.null(input$volcano_label_mode)) input$volcano_label_mode else "top_fc",
        label_distribution = if (!is.null(input$volcano_balance_mode)) input$volcano_balance_mode else "per_side",
        only_sig_labels = if (!is.null(input$volcano_sig_only)) input$volcano_sig_only else TRUE,
        label_size = if (!is.null(input$volcano_label_size)) input$volcano_label_size else 3.5,
        symmetric_x = input$volcano_sym_x,
        show_cutoffs = input$volcano_show_lines,
        title = if (!is.null(user_title)) user_title else "Differential Expression Volcano Plot",
        subtitle = user_sub,
        xlab = if (!is.null(user_xlab)) user_xlab else "log2(Fold Change)",
        ylab = user_ylab,
        theme_name = th_name,
        base_size = b_size
      )
    } else if (ptype == "protein_box") {
      req(input$target_gene_select)
      show_l <- if (!is.null(input$profile_show_line)) input$profile_show_line else TRUE
      show_pt <- if (!is.null(input$profile_show_points)) input$profile_show_points else TRUE
      p <- plot_protein_abundance(
        df = df,
        gene_name = input$target_gene_select,
        plot_type = input$box_style,
        group_mapping = active_group_mapping(),
        palette = if (palette_choice %in% names(JOURNAL_PALETTES)) palette_choice else "NPG",
        show_line = show_l,
        show_points = show_pt,
        title = user_title,
        theme_name = th_name,
        base_size = b_size
      )
    } else if (ptype == "heatmap") {
      top_n_val <- if (!is.null(input$heatmap_top_n) && input$heatmap_top_n > 0) input$heatmap_top_n else NULL
      p <- plot_expression_heatmap(
        df = df,
        group_mapping = active_group_mapping(),
        scale = if (!is.null(input$heatmap_scale)) input$heatmap_scale else "row",
        cluster_rows = if (!is.null(input$heatmap_cluster_rows)) input$heatmap_cluster_rows else TRUE,
        cluster_cols = if (!is.null(input$heatmap_cluster_cols)) input$heatmap_cluster_cols else FALSE,
        show_rownames = if (!is.null(input$heatmap_show_rownames)) input$heatmap_show_rownames else TRUE,
        show_colnames = if (!is.null(input$heatmap_show_colnames)) input$heatmap_show_colnames else TRUE,
        top_n = top_n_val,
        palette = if (!is.null(input$heatmap_color_palette)) input$heatmap_color_palette else "RdBu",
        title = user_title
      )
    } else if (ptype == "tf_bar") {
      p <- plot_tf_barplot(
        df = df,
        top_n = input$tf_top_n,
        selected_tfs = input$tf_selected_tfs,
        filter_mode = input$tf_filter_mode,
        act_color = up_c,
        rep_color = down_c,
        orientation = input$tf_orientation,
        title = if (!is.null(user_title)) user_title else "Transcription Factor Footprinting Activities",
        subtitle = user_sub,
        xlab = if (!is.null(user_xlab)) user_xlab else "TF Activity Score",
        ylab = if (!is.null(user_ylab)) user_ylab else "Transcription Factor",
        theme_name = th_name,
        base_size = b_size
      )
    } else if (ptype == "tf_lollipop") {
      p <- plot_tf_lollipop(
        df = df,
        top_n = input$tf_top_n,
        selected_tfs = input$tf_selected_tfs,
        filter_mode = input$tf_filter_mode,
        act_color = up_c,
        rep_color = down_c,
        title = if (!is.null(user_title)) user_title else "Transcription Factor Activity Lollipop",
        subtitle = user_sub,
        theme_name = th_name,
        base_size = b_size
      )
    }

    return(p)
  })

  # --- 5. RENDER PLOT ---
  output$main_plot <- renderPlot({
    p <- current_plot()
    req(p)
    if (inherits(p, "gtable") || inherits(p, "grob")) {
      grid::grid.newpage()
      grid::grid.draw(p)
    } else {
      print(p)
    }
  }, res = 100)

  # --- 6. RENDER SUMMARY HEADER & DATA TABLE ---
  output$table_summary_header <- renderUI({
    df <- raw_data()
    req(df)
    dtype <- data_type()
    req(dtype)
    ptype <- input$selected_plot_type
    req(ptype)

    valid_plot_map <- list(
      "dea" = c("volcano", "protein_box", "heatmap"),
      "enrich" = c("dotplot", "barplot", "lollipop", "gsea_diverging"),
      "gsea" = c("dotplot", "barplot", "lollipop", "gsea_diverging"),
      "tf" = c("tf_bar", "tf_lollipop")
    )
    if (!is.null(valid_plot_map[[dtype]]) && !ptype %in% valid_plot_map[[dtype]]) {
      return(NULL)
    }

    if (ptype == "gsea_diverging" && "NES" %in% colnames(df)) {
      n_act <- sum(df$NES > 0, na.rm = TRUE)
      n_supp <- sum(df$NES < 0, na.rm = TRUE)
      div(
        class = "mb-3 p-2 bg-light rounded border d-flex gap-3 align-items-center",
        span(class = "badge bg-danger p-2 fs-6", style = "font-weight: 700;", paste("🔴 Activated Pathways (NES > 0):", n_act)),
        span(class = "badge bg-primary p-2 fs-6", style = "font-weight: 700;", paste("🔵 Suppressed Pathways (NES < 0):", n_supp)),
        span(class = "text-muted ms-auto", style = "font-size: 0.85rem;", "The first column indicates regulation Status.")
      )
    } else if (ptype %in% c("tf_bar", "tf_lollipop") && "score" %in% colnames(df)) {
      n_act <- sum(df$score > 0, na.rm = TRUE)
      n_rep <- sum(df$score < 0, na.rm = TRUE)
      div(
        class = "mb-3 p-2 bg-light rounded border d-flex gap-3 align-items-center",
        span(class = "badge bg-danger p-2 fs-6", style = "font-weight: 700;", paste("🔴 Activated TFs:", n_act)),
        span(class = "badge bg-primary p-2 fs-6", style = "font-weight: 700;", paste("🔵 Repressed TFs:", n_rep))
      )
    } else if (ptype == "volcano") {
      dea_df <- standardize_dea_df(df)
      p_col <- if (input$volcano_p_type == "P.Value") "P.Value" else "adj.P.Val"
      p_use <- pmax(dea_df[[p_col]], 1e-300)
      fc_cut <- input$volcano_fc_cutoff
      p_cut <- input$volcano_p_cutoff
      n_up <- sum(dea_df$logFC >= fc_cut & p_use <= p_cut, na.rm = TRUE)
      n_down <- sum(dea_df$logFC <= -fc_cut & p_use <= p_cut, na.rm = TRUE)
      n_ns <- nrow(dea_df) - n_up - n_down
      div(
        class = "mb-3 p-2 bg-light rounded border d-flex gap-3 align-items-center",
        span(class = "badge bg-danger p-2 fs-6", style = "font-weight: 700;", paste("🔴 Enriched / Up (logFC ≥", fc_cut, "):", n_up)),
        span(class = "badge bg-primary p-2 fs-6", style = "font-weight: 700;", paste("🔵 Reduced / Down (logFC ≤ -", fc_cut, "):", n_down)),
        span(class = "badge bg-secondary p-2 fs-6", style = "font-weight: 700;", paste("⚪ Not Significant:", n_ns)),
        span(class = "text-muted ms-auto", style = "font-size: 0.85rem;", paste("Filter: |logFC| ≥", fc_cut, "& P ≤", p_cut))
      )
    } else if (ptype == "heatmap") {
      dea_df <- standardize_dea_df(df)
      s_cols <- current_sample_cols()
      grps <- active_group_mapping()
      n_grps <- if (!is.null(grps)) length(unique(grps)) else 0
      div(
        class = "mb-3 p-2 bg-light rounded border d-flex gap-3 align-items-center",
        span(class = "badge bg-danger p-2 fs-6", style = "font-weight: 700;", paste("Total Proteins:", nrow(dea_df))),
        span(class = "badge bg-primary p-2 fs-6", style = "font-weight: 700;", paste("Number of Samples:", length(s_cols))),
        span(class = "badge bg-success p-2 fs-6", style = "font-weight: 700;", paste("Number of Groups:", n_grps)),
        if ("Comparison" %in% colnames(dea_df) && !is.na(dea_df$Comparison[1])) span(class = "badge bg-secondary p-2 fs-6", style = "font-weight: 700;", paste("Comparison:", dea_df$Comparison[1])) else NULL
      )
    } else {
      NULL
    }
  })

  output$plot_data_table <- DT::renderDataTable({
    df <- raw_data()
    req(df)
    dtype <- data_type()
    req(dtype)
    ptype <- input$selected_plot_type
    req(ptype)

    valid_plot_map <- list(
      "dea" = c("volcano", "protein_box", "heatmap"),
      "enrich" = c("dotplot", "barplot", "lollipop", "gsea_diverging"),
      "gsea" = c("dotplot", "barplot", "lollipop", "gsea_diverging"),
      "tf" = c("tf_bar", "tf_lollipop")
    )
    if (!is.null(valid_plot_map[[dtype]]) && !ptype %in% valid_plot_map[[dtype]]) {
      return(NULL)
    }

    # Filter table to match current plot filtering if applicable
    table_df <- df
    if (ptype %in% c("dotplot", "barplot", "lollipop") && "Description" %in% colnames(df)) {
      sel <- input$enrich_selected_terms
      fmode <- input$enrich_filter_mode
      if (!is.null(sel) && length(sel) > 0) {
        if (fmode == "exclude") table_df <- table_df %>% filter(!Description %in% sel)
        else table_df <- table_df %>% filter(Description %in% sel)
      }
    } else if (ptype == "gsea_diverging" && "Description" %in% colnames(df)) {
      sel <- input$gsea_selected_terms
      fmode <- input$gsea_filter_mode
      if (!is.null(sel) && length(sel) > 0) {
        if (fmode == "exclude") table_df <- table_df %>% filter(!Description %in% sel)
        else table_df <- table_df %>% filter(Description %in% sel)
      }
      if ("NES" %in% colnames(table_df)) {
        table_df$Status <- ifelse(table_df$NES > 0, "Activated", "Suppressed")
        first_cols <- intersect(c("Status", "Description", "NES", "p.adjust", "pvalue", "setSize", "Count"), colnames(table_df))
        rest_cols <- setdiff(colnames(table_df), first_cols)
        table_df <- table_df[, c(first_cols, rest_cols)]
      }
    } else if (ptype == "volcano") {
      dea_df <- standardize_dea_df(table_df)
      p_col <- if (input$volcano_p_type == "P.Value") "P.Value" else "adj.P.Val"
      p_use <- pmax(dea_df[[p_col]], 1e-300)
      fc_cut <- input$volcano_fc_cutoff
      p_cut <- input$volcano_p_cutoff

      reg <- rep("Not Significant", nrow(dea_df))
      reg[dea_df$logFC >= fc_cut & p_use <= p_cut] <- "Enriched / Up"
      reg[dea_df$logFC <= -fc_cut & p_use <= p_cut] <- "Reduced / Down"

      table_df$Regulation <- reg
      first_cols <- intersect(c("Regulation", "Gene", "logFC", "adj.P.Val", "P.Value", "AveExpr"), colnames(table_df))
      rest_cols <- setdiff(colnames(table_df), first_cols)
      table_df <- table_df[, c(first_cols, rest_cols)]
    } else if (ptype == "heatmap") {
      dea_df <- standardize_dea_df(table_df)
      first_cols <- intersect(c("Gene", "Display_Name", "logFC", "adj.P.Val", "P.Value", "Comparison"), colnames(dea_df))
      rest_cols <- setdiff(colnames(dea_df), first_cols)
      table_df <- dea_df[, c(first_cols, rest_cols)]
    } else if (ptype %in% c("tf_bar", "tf_lollipop")) {
      tf_col <- grep("^(source|tf)$", colnames(df), ignore.case = TRUE, value = TRUE)
      if (length(tf_col) > 0) {
        sel <- input$tf_selected_tfs
        fmode <- input$tf_filter_mode
        if (!is.null(sel) && length(sel) > 0) {
          if (fmode == "exclude") table_df <- table_df[!table_df[[tf_col[1]]] %in% sel, ]
          else table_df <- table_df[table_df[[tf_col[1]]] %in% sel, ]
        }
      }
    }

    dt <- DT::datatable(
      table_df,
      options = list(pageLength = 10, scrollX = TRUE, autoWidth = TRUE),
      rownames = FALSE,
      class = "compact stripe hover"
    )

    if ("Status" %in% colnames(table_df)) {
      dt <- dt %>% DT::formatStyle(
        "Status",
        backgroundColor = DT::styleEqual(c("Activated", "Suppressed"), c("#FDF0ED", "#EEF2F7")),
        color = DT::styleEqual(c("Activated", "Suppressed"), c("#C0392B", "#2C3E50")),
        fontWeight = "bold"
      )
    }
    if ("Regulation" %in% colnames(table_df)) {
      dt <- dt %>% DT::formatStyle(
        "Regulation",
        backgroundColor = DT::styleEqual(c("Enriched / Up", "Reduced / Down", "Not Significant"), c("#FDF0ED", "#EEF2F7", "#F8F9FA")),
        color = DT::styleEqual(c("Enriched / Up", "Reduced / Down", "Not Significant"), c("#C0392B", "#2C3E50", "#7F8C8D")),
        fontWeight = "bold"
      )
    }

    dt
  })

  # --- 7. DOWNLOAD HANDLER ---
  output$download_fig_btn <- downloadHandler(
    filename = function() {
      fmt <- input$export_format
      ext <- switch(fmt,
        "png_300" = "png",
        "png_600" = "png",
        "pdf" = "pdf",
        "svg" = "svg",
        "tiff" = "tiff",
        "png"
      )
      paste0("Figure_", input$selected_plot_type, "_", Sys.Date(), ".", ext)
    },
    content = function(file) {
      p <- current_plot()
      req(p)

      fmt <- input$export_format
      w <- input$export_width
      h <- input$export_height

      dpi <- 300
      out_ext <- "png"

      if (fmt == "png_300") {
        dpi <- 300; out_ext <- "png"
      } else if (fmt == "png_600") {
        dpi <- 600; out_ext <- "png"
      } else if (fmt == "pdf") {
        out_ext <- "pdf"
      } else if (fmt == "svg") {
        out_ext <- "svg"
      } else if (fmt == "tiff") {
        dpi <- 300; out_ext <- "tiff"
      }

      export_figure(
        plot_obj = p,
        file_path = file,
        format = out_ext,
        width = w,
        height = h,
        units = "in",
        dpi = dpi
      )
    }
  )
}

shinyApp(ui = ui, server = server)
