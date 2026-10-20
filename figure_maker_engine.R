# ==============================================================================
# ProteoPrism - Core Plotting Engine
# High-Quality, Publication-Ready Visualization Functions for Differential Proteomics
# Compatible with outputs from:
# 1. Proteomics_Pipeline.Rmd (DEA, Quant, PCA, Heatmaps, Boxplots)
# 2. Ultimate_Functional_Analysis.Rmd (GO, KEGG, Reactome, Hallmark, GSEA)
# 3. Proteomic_Footprinting.Rmd (DoRothEA TF Activities)
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(scales)
  if (requireNamespace("ggrepel", quietly = TRUE)) library(ggrepel)
  if (requireNamespace("ggsci", quietly = TRUE)) library(ggsci)
  if (requireNamespace("pheatmap", quietly = TRUE)) library(pheatmap)
  if (requireNamespace("RColorBrewer", quietly = TRUE)) library(RColorBrewer)
  if (requireNamespace("viridis", quietly = TRUE)) library(viridis)
})

# ------------------------------------------------------------------------------
# 1. PALETTE DEFINITIONS & THEME UTILITIES
# ------------------------------------------------------------------------------

JOURNAL_PALETTES <- list(
  "NPG" = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"),
  "NEJM" = c("#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1", "#6F99AD", "#FFDC91", "#EE4C97"),
  "Lancet" = c("#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F", "#FDAF91", "#AD002A", "#ADB6B6", "#1B1919"),
  "JCO" = c("#0073C2", "#EFC000", "#868686", "#CD534C", "#7AA6DC", "#003C67", "#8F7700", "#3B3B3B", "#A73030"),
  "AAAS" = c("#3B4992", "#EE0000", "#008B45", "#631879", "#008280", "#BB0021", "#5F559B", "#A20056"),
  "Set1" = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF"),
  "Set2" = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"),
  "Dark2" = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
)

get_theme <- function(theme_name = "classic", base_size = 12, font_family = "") {
  t <- switch(theme_name,
    "classic" = theme_classic(base_size = base_size, base_family = font_family),
    "minimal" = theme_minimal(base_size = base_size, base_family = font_family),
    "bw" = theme_bw(base_size = base_size, base_family = font_family),
    "prism" = {
      theme_classic(base_size = base_size, base_family = font_family) +
        theme(
          axis.line = element_line(linewidth = 0.8, color = "black"),
          axis.ticks = element_line(linewidth = 0.8, color = "black"),
          panel.grid = element_blank(),
          strip.background = element_rect(fill = "grey92", color = "black", linewidth = 0.8)
        )
    },
    theme_classic(base_size = base_size, base_family = font_family)
  )
  
  t + theme(
    plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 8)),
    plot.subtitle = element_text(hjust = 0.5, margin = margin(b = 8)),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold")
  )
}

parse_generatio <- function(vec) {
  if (is.numeric(vec)) return(vec)
  sapply(as.character(vec), function(x) {
    if (is.na(x) || x == "") return(NA_real_)
    if (grepl("/", x)) {
      parts <- as.numeric(strsplit(x, "/")[[1]])
      if (length(parts) == 2 && parts[2] != 0) return(parts[1] / parts[2])
    }
    val <- suppressWarnings(as.numeric(x))
    return(val)
  }, USE.NAMES = FALSE)
}

format_pval_label <- function(pval, padj = NULL, mode = "fdr") {
  if (mode == "fdr") {
    val <- if (!is.null(padj) && !is.na(padj)) padj else pval
    if (is.na(val)) return("")
    if (val < 0.001) sprintf("FDR=%.1e", val) else sprintf("FDR=%.3f", val)
  } else if (mode == "pval") {
    if (is.na(pval)) return("")
    if (pval < 0.001) sprintf("p=%.1e", pval) else sprintf("p=%.3f", pval)
  } else if (mode == "stars") {
    val <- if (!is.null(padj) && !is.na(padj)) padj else pval
    if (is.na(val)) return("")
    if (val < 0.001) "***" else if (val < 0.01) "**" else if (val < 0.05) "*" else "ns"
  } else if (mode == "both") {
    p_str <- if (!is.na(pval)) (if (pval < 0.001) sprintf("p=%.1e", pval) else sprintf("p=%.3f", pval)) else ""
    fdr_str <- if (!is.null(padj) && !is.na(padj)) (if (padj < 0.001) sprintf("FDR=%.1e", padj) else sprintf("FDR=%.3f", padj)) else ""
    paste(p_str, fdr_str, sep = "\n")
  } else {
    ""
  }
}

# ------------------------------------------------------------------------------
# 2. MODULE 2: FUNCTIONAL & PATHWAY ENRICHMENT PLOTS
# ------------------------------------------------------------------------------

standardize_enrich_df <- function(df) {
  df <- as.data.frame(df)
  if (!"Description" %in% colnames(df)) {
    desc_cols <- grep("^(term|name|pathway|description|id)$", colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(desc_cols) > 0) df$Description <- df[[desc_cols[1]]]
  }
  
  if (!"p.adjust" %in% colnames(df)) {
    padj_cols <- grep("^(padj|p_adjust|p.adj|fdr|qvalue|q_value|adj\\.p\\.val|adj_pval|sca\\.adj\\.pval)$", colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(padj_cols) > 0) df$p.adjust <- as.numeric(df[[padj_cols[1]]])
    else if ("pvalue" %in% colnames(df)) df$p.adjust <- as.numeric(df$pvalue)
    else if ("P.Value" %in% colnames(df)) df$p.adjust <- as.numeric(df$P.Value)
    else df$p.adjust <- rep(0.05, nrow(df))
  } else {
    df$p.adjust <- as.numeric(df$p.adjust)
  }
  
  if (is.null(df$p.adjust) || length(df$p.adjust) == 0) {
    df$p.adjust <- rep(0.05, nrow(df))
  }
  
  if (!"pvalue" %in% colnames(df)) {
    pval_cols <- grep("^(pval|p_val|p.val|pvalue|p\\.value|sca\\.p\\.value)$", colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(pval_cols) > 0) df$pvalue <- as.numeric(df[[pval_cols[1]]])
    else if ("p.adjust" %in% colnames(df)) df$pvalue <- df$p.adjust
    else df$pvalue <- rep(0.05, nrow(df))
  } else {
    df$pvalue <- as.numeric(df$pvalue)
  }
  
  if (is.null(df$pvalue) || length(df$pvalue) == 0) {
    df$pvalue <- df$p.adjust
  }
  
  if (!"Count" %in% colnames(df)) {
    cnt_cols <- grep("^(count|setsize|size|genesize)$", colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(cnt_cols) > 0) df$Count <- as.numeric(df[[cnt_cols[1]]])
    else df$Count <- 10
  } else {
    df$Count <- as.numeric(df$Count)
  }
  
  if ("GeneRatio" %in% colnames(df)) {
    df$GeneRatio_num <- parse_generatio(df$GeneRatio)
  } else {
    df$GeneRatio_num <- df$Count / max(df$Count, na.rm = TRUE)
  }
  
  if (!"FoldEnrichment" %in% colnames(df)) {
    fe_cols <- grep("^(foldenrichment|fold_enrichment|fold.enrichment|enrichment)$", colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(fe_cols) > 0) df$FoldEnrichment <- as.numeric(df[[fe_cols[1]]])
    else df$FoldEnrichment <- NA_real_
  } else {
    df$FoldEnrichment <- as.numeric(df$FoldEnrichment)
  }
  
  if (!"RichFactor" %in% colnames(df)) {
    rf_cols <- grep("^(richfactor|rich_factor|rich.factor)$", colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(rf_cols) > 0) df$RichFactor <- as.numeric(df[[rf_cols[1]]])
    else df$RichFactor <- NA_real_
  } else {
    df$RichFactor <- as.numeric(df$RichFactor)
  }
  
  df$negLog10Padj <- -log10(pmax(df$p.adjust, 1e-50))
  
  if ("NES" %in% colnames(df)) {
    df$NES <- as.numeric(df$NES)
  }
  
  return(df)
}

#' Plot Functional Enrichment Dot Plot
plot_enrichment_dot <- function(df,
                                x_metric = "GeneRatio",
                                color_metric = "p.adjust",
                                size_metric = "Count",
                                top_n = 15,
                                selected_terms = NULL,
                                filter_mode = "exclude",
                                sort_by = "p.adjust",
                                palette_type = "gradient",
                                gradient_low = "#E64B35",
                                gradient_high = "#3C5488",
                                viridis_opt = "D",
                                wrap_length = 50,
                                point_size_range = c(3, 9),
                                p_limits = NULL,
                                show_p_labels = FALSE,
                                title = "Functional Enrichment Analysis",
                                subtitle = NULL,
                                xlab = NULL,
                                ylab = "",
                                theme_name = "classic",
                                base_size = 12) {
  data <- standardize_enrich_df(df)
  
  if (!is.null(selected_terms) && length(selected_terms) > 0) {
    if (filter_mode == "exclude") {
      data <- data %>% filter(!Description %in% selected_terms)
    } else {
      data <- data %>% filter(Description %in% selected_terms)
    }
  }
  
  if (sort_by == "p.adjust") {
    data <- data %>% arrange(p.adjust)
  } else if (sort_by == "Count") {
    data <- data %>% arrange(desc(Count))
  } else if (sort_by == "FoldEnrichment" && !all(is.na(data$FoldEnrichment))) {
    data <- data %>% arrange(desc(FoldEnrichment))
  } else if (sort_by == "NES" && "NES" %in% colnames(data)) {
    data <- data %>% arrange(desc(abs(NES)))
  } else {
    data <- data %>% arrange(p.adjust)
  }
  data <- head(data, top_n)
  
  if (nrow(data) == 0) stop("No terms match the selected criteria.")
  
  x_col <- switch(x_metric,
    "GeneRatio" = "GeneRatio_num",
    "Count" = "Count",
    "FoldEnrichment" = if ("FoldEnrichment" %in% colnames(data) && !all(is.na(data$FoldEnrichment))) "FoldEnrichment" else "GeneRatio_num",
    "RichFactor" = if ("RichFactor" %in% colnames(data) && !all(is.na(data$RichFactor))) "RichFactor" else "GeneRatio_num",
    "negLog10Padj" = "negLog10Padj",
    "NES" = if ("NES" %in% colnames(data)) "NES" else "GeneRatio_num",
    "GeneRatio_num"
  )
  if (is.null(xlab)) {
    xlab <- switch(x_metric,
      "GeneRatio" = "Gene Ratio",
      "Count" = "Gene Count",
      "FoldEnrichment" = "Fold Enrichment",
      "RichFactor" = "Rich Factor",
      "negLog10Padj" = "-log10(p.adjust)",
      "NES" = "Normalized Enrichment Score (NES)",
      x_metric
    )
  }
  
  data$TermWrapped <- stringr::str_wrap(data$Description, width = wrap_length)
  data$TermWrapped <- factor(data$TermWrapped, levels = rev(data$TermWrapped))
  
  size_col <- if (size_metric %in% colnames(data)) size_metric else "Count"
  
  p <- ggplot(data, aes(x = .data[[x_col]], y = .data[["TermWrapped"]], size = .data[[size_col]], color = .data[[color_metric]])) +
    geom_point(alpha = 0.9) +
    scale_size_continuous(range = point_size_range, name = size_metric) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    get_theme(theme_name, base_size = base_size)
  
  if (palette_type == "gradient") {
    if (color_metric %in% c("p.adjust", "pvalue")) {
      p <- p + scale_color_gradient(low = gradient_low, high = gradient_high, name = color_metric, limits = p_limits)
    } else {
      p <- p + scale_color_gradient(low = gradient_high, high = gradient_low, name = color_metric, limits = p_limits)
    }
  } else if (palette_type == "viridis") {
    p <- p + scale_color_viridis_c(option = viridis_opt, direction = if (color_metric %in% c("p.adjust", "pvalue")) -1 else 1, name = color_metric, limits = p_limits)
  } else if (palette_type == "brewer") {
    p <- p + scale_color_distiller(palette = "RdYlBu", direction = if (color_metric %in% c("p.adjust", "pvalue")) 1 else -1, name = color_metric, limits = p_limits)
  }
  
  if (show_p_labels) {
    data$P_Text <- mapply(format_pval_label, data$pvalue, data$p.adjust, MoreArgs = list(mode = "fdr"))
    p <- p + geom_text(data = data, aes(x = .data[[x_col]], y = .data[["TermWrapped"]], label = P_Text),
                       hjust = -0.3, vjust = 0.5, size = 3.2, color = "black", inherit.aes = FALSE) +
      scale_x_continuous(expand = expansion(mult = c(0.05, 0.2)))
  }
  
  return(p)
}

#' Plot Functional Enrichment Bar Plot
plot_enrichment_bar <- function(df,
                                x_metric = "Count",
                                color_metric = "p.adjust",
                                top_n = 15,
                                selected_terms = NULL,
                                filter_mode = "exclude",
                                sort_by = "p.adjust",
                                orientation = "horizontal",
                                palette_type = "gradient",
                                gradient_low = "#E64B35",
                                gradient_high = "#3C5488",
                                viridis_opt = "D",
                                wrap_length = 50,
                                p_limits = NULL,
                                show_p_labels = FALSE,
                                title = "Functional Enrichment Barplot",
                                subtitle = NULL,
                                xlab = NULL,
                                ylab = "",
                                theme_name = "classic",
                                base_size = 12) {
  data <- standardize_enrich_df(df)
  
  if (!is.null(selected_terms) && length(selected_terms) > 0) {
    if (filter_mode == "exclude") {
      data <- data %>% filter(!Description %in% selected_terms)
    } else {
      data <- data %>% filter(Description %in% selected_terms)
    }
  }
  
  if (sort_by == "p.adjust") data <- data %>% arrange(p.adjust)
  else if (sort_by == "Count") data <- data %>% arrange(desc(Count))
  else if (sort_by == "FoldEnrichment") data <- data %>% arrange(desc(FoldEnrichment))
  else data <- data %>% arrange(p.adjust)
  data <- head(data, top_n)
  
  if (nrow(data) == 0) stop("No terms match the selected criteria.")
  
  val_col <- switch(x_metric,
    "Count" = "Count",
    "GeneRatio" = "GeneRatio_num",
    "FoldEnrichment" = "FoldEnrichment",
    "RichFactor" = "RichFactor",
    "negLog10Padj" = "negLog10Padj",
    "Count"
  )
  if (is.null(xlab)) {
    xlab <- switch(x_metric,
      "Count" = "Gene Count",
      "GeneRatio" = "Gene Ratio",
      "FoldEnrichment" = "Fold Enrichment",
      "RichFactor" = "Rich Factor",
      "negLog10Padj" = "-log10(p.adjust)",
      x_metric
    )
  }
  
  data$TermWrapped <- stringr::str_wrap(data$Description, width = wrap_length)
  data$TermWrapped <- factor(data$TermWrapped, levels = rev(data$TermWrapped))
  
  p <- ggplot(data, aes(x = .data[[val_col]], y = .data[["TermWrapped"]], fill = .data[[color_metric]])) +
    geom_col(width = 0.7, color = "black", linewidth = 0.3) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    get_theme(theme_name, base_size = base_size)
  
  if (palette_type == "gradient") {
    p <- p + scale_fill_gradient(low = gradient_low, high = gradient_high, name = color_metric, limits = p_limits)
  } else if (palette_type == "viridis") {
    p <- p + scale_fill_viridis_c(option = viridis_opt, direction = -1, name = color_metric, limits = p_limits)
  } else if (palette_type == "brewer") {
    p <- p + scale_fill_distiller(palette = "YlOrRd", direction = 1, name = color_metric, limits = p_limits)
  }
  
  if (show_p_labels) {
    data$P_Text <- mapply(format_pval_label, data$pvalue, data$p.adjust, MoreArgs = list(mode = "fdr"))
    p <- p + geom_text(data = data, aes(x = .data[[val_col]], y = .data[["TermWrapped"]], label = P_Text),
                       hjust = -0.2, vjust = 0.5, size = 3.2, color = "black", inherit.aes = FALSE) +
      scale_x_continuous(expand = expansion(mult = c(0.05, 0.2)))
  }
  
  if (orientation == "vertical") {
    p <- p + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  return(p)
}

#' Plot Functional Enrichment Lollipop Chart
plot_enrichment_lollipop <- function(df,
                                     x_metric = "Count",
                                     color_metric = "p.adjust",
                                     top_n = 15,
                                     selected_terms = NULL,
                                     filter_mode = "exclude",
                                     sort_by = "p.adjust",
                                     point_size = 5,
                                     palette_type = "gradient",
                                     gradient_low = "#E64B35",
                                     gradient_high = "#3C5488",
                                     wrap_length = 50,
                                     p_limits = NULL,
                                     show_p_labels = FALSE,
                                     title = "Enrichment Lollipop Chart",
                                     subtitle = NULL,
                                     xlab = NULL,
                                     ylab = "",
                                     theme_name = "classic",
                                     base_size = 12) {
  data <- standardize_enrich_df(df)
  
  if (!is.null(selected_terms) && length(selected_terms) > 0) {
    if (filter_mode == "exclude") {
      data <- data %>% filter(!Description %in% selected_terms)
    } else {
      data <- data %>% filter(Description %in% selected_terms)
    }
  }
  
  if (sort_by == "p.adjust") data <- data %>% arrange(p.adjust)
  else data <- data %>% arrange(desc(Count))
  data <- head(data, top_n)
  
  if (nrow(data) == 0) stop("No terms match the selected criteria.")
  
  val_col <- switch(x_metric,
    "Count" = "Count",
    "GeneRatio" = "GeneRatio_num",
    "FoldEnrichment" = "FoldEnrichment",
    "RichFactor" = "RichFactor",
    "negLog10Padj" = "negLog10Padj",
    "Count"
  )
  if (is.null(xlab)) xlab <- val_col
  
  data$TermWrapped <- stringr::str_wrap(data$Description, width = wrap_length)
  data$TermWrapped <- factor(data$TermWrapped, levels = rev(data$TermWrapped))
  
  p <- ggplot(data, aes(x = .data[[val_col]], y = .data[["TermWrapped"]])) +
    geom_segment(aes(x = 0, xend = .data[[val_col]], y = .data[["TermWrapped"]], yend = .data[["TermWrapped"]]), color = "grey60", linewidth = 0.8) +
    geom_point(aes(color = .data[[color_metric]]), size = point_size) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    get_theme(theme_name, base_size = base_size)
  
  if (palette_type == "gradient") {
    p <- p + scale_color_gradient(low = gradient_low, high = gradient_high, name = color_metric, limits = p_limits)
  } else {
    p <- p + scale_color_viridis_c(direction = -1, name = color_metric, limits = p_limits)
  }
  
  if (show_p_labels) {
    data$P_Text <- mapply(format_pval_label, data$pvalue, data$p.adjust, MoreArgs = list(mode = "fdr"))
    p <- p + geom_text(data = data, aes(x = .data[[val_col]], y = .data[["TermWrapped"]], label = P_Text),
                       hjust = -0.3, vjust = 0.5, size = 3.2, color = "black", inherit.aes = FALSE) +
      scale_x_continuous(expand = expansion(mult = c(0.05, 0.2)))
  }
  
  return(p)
}

#' Plot GSEA Diverging NES Lollipop Chart (Colored by FDR, with Top Headers)
plot_gsea_diverging <- function(df,
                                top_n = 10,
                                selected_terms = NULL,
                                filter_mode = "exclude",
                                color_by = "fdr",          # "fdr", "pvalue", "status"
                                palette_type = "gradient", # "gradient", "viridis", "brewer"
                                gradient_low = "#E64B35",  # Most significant FDR (Red)
                                gradient_high = "#3C5488", # Higher FDR (Navy/Blue)
                                viridis_opt = "D",
                                p_limits = NULL,
                                up_color = "#E64B35",      # If color_by == "status"
                                down_color = "#3C5488",    # If color_by == "status"
                                point_size = 5.0,
                                size_by_set = TRUE,
                                show_p_labels = TRUE,
                                p_label_type = "fdr",
                                p_label_size = 3.5,
                                show_top_headers = TRUE,   # "Activated" and "Suppressed" at top of plot
                                header_style = "banner",   # "banner" or "text"
                                wrap_length = 45,
                                title = "GSEA Normalized Enrichment Scores (NES)",
                                subtitle = NULL,
                                theme_name = "classic",
                                base_size = 12) {
  data <- standardize_enrich_df(df)
  if (!"NES" %in% colnames(data)) stop("GSEA data must contain 'NES' column.")
  
  if (!is.null(selected_terms) && length(selected_terms) > 0) {
    if (filter_mode == "exclude") {
      data <- data %>% filter(!Description %in% selected_terms)
    } else {
      data <- data %>% filter(Description %in% selected_terms)
    }
  }
  
  pos <- data %>% filter(NES > 0) %>% arrange(desc(NES)) %>% head(top_n)
  neg <- data %>% filter(NES < 0) %>% arrange(NES) %>% head(top_n)
  plot_df <- bind_rows(pos, neg)
  
  if (nrow(plot_df) == 0) stop("No GSEA terms match criteria.")
  
  plot_df <- plot_df %>% arrange(NES)
  plot_df$Status <- ifelse(plot_df$NES > 0, "Activated", "Suppressed")
  plot_df$TermWrapped <- stringr::str_wrap(plot_df$Description, width = wrap_length)
  plot_df$TermWrapped <- factor(plot_df$TermWrapped, levels = plot_df$TermWrapped)
  
  # Format p-value labels
  plot_df$P_Label <- mapply(format_pval_label, plot_df$pvalue, plot_df$p.adjust, MoreArgs = list(mode = p_label_type))
  
  p <- ggplot(plot_df, aes(x = NES, y = TermWrapped)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6)
  
  # Stems
  if (color_by == "status") {
    p <- p + geom_segment(aes(x = 0, xend = NES, y = TermWrapped, yend = TermWrapped, color = Status), linewidth = 0.9)
  } else {
    p <- p + geom_segment(aes(x = 0, xend = NES, y = TermWrapped, yend = TermWrapped), color = "grey65", linewidth = 0.9)
  }
  
  # Points & Color Mapping
  if (color_by == "fdr") {
    if (size_by_set && "Count" %in% colnames(plot_df)) {
      p <- p + geom_point(aes(color = p.adjust, size = Count), alpha = 0.95) +
        scale_size_continuous(range = c(3.5, 8.5), name = "Gene Count")
    } else {
      p <- p + geom_point(aes(color = p.adjust), size = point_size, alpha = 0.95)
    }
    
    if (palette_type == "gradient") {
      p <- p + scale_color_gradient(low = gradient_low, high = gradient_high, name = "FDR", limits = p_limits)
    } else if (palette_type == "viridis") {
      p <- p + scale_color_viridis_c(option = viridis_opt, direction = -1, name = "FDR", limits = p_limits)
    } else if (palette_type == "brewer") {
      p <- p + scale_color_distiller(palette = "RdYlBu", direction = 1, name = "FDR", limits = p_limits)
    }
  } else if (color_by == "pvalue") {
    if (size_by_set && "Count" %in% colnames(plot_df)) {
      p <- p + geom_point(aes(color = pvalue, size = Count), alpha = 0.95) +
        scale_size_continuous(range = c(3.5, 8.5), name = "Gene Count")
    } else {
      p <- p + geom_point(aes(color = pvalue), size = point_size, alpha = 0.95)
    }
    
    if (palette_type == "gradient") {
      p <- p + scale_color_gradient(low = gradient_low, high = gradient_high, name = "p-value", limits = p_limits)
    } else if (palette_type == "viridis") {
      p <- p + scale_color_viridis_c(option = viridis_opt, direction = -1, name = "p-value", limits = p_limits)
    } else {
      p <- p + scale_color_distiller(palette = "RdYlBu", direction = 1, name = "p-value", limits = p_limits)
    }
  } else {
    # Status coloring
    if (size_by_set && "Count" %in% colnames(plot_df)) {
      p <- p + geom_point(aes(color = Status, size = Count), alpha = 0.95) +
        scale_size_continuous(range = c(3.5, 8.5), name = "Gene Count")
    } else {
      p <- p + geom_point(aes(color = Status), size = point_size, alpha = 0.95)
    }
    p <- p + scale_color_manual(values = c("Activated" = up_color, "Suppressed" = down_color), name = "Status")
  }
  
  p <- p +
    labs(title = title, subtitle = subtitle, x = "Normalized Enrichment Score (NES)", y = "") +
    get_theme(theme_name, base_size = base_size)
  
  # Add p-value text labels next to the dots
  if (show_p_labels) {
    pos_sub <- plot_df %>% filter(NES > 0)
    neg_sub <- plot_df %>% filter(NES < 0)
    
    if (nrow(pos_sub) > 0) {
      p <- p + geom_text(data = pos_sub, aes(x = NES, y = TermWrapped, label = P_Label),
                         hjust = -0.3, vjust = 0.5, size = p_label_size, fontface = "bold", color = "black")
    }
    if (nrow(neg_sub) > 0) {
      # Position slightly above stem towards 0 to avoid overlapping with left y-axis text
      p <- p + geom_text(data = neg_sub, aes(x = NES, y = TermWrapped, label = P_Label),
                         hjust = -0.2, vjust = -0.7, size = p_label_size, fontface = "bold", color = "black")
    }
    
    p <- p + scale_x_continuous(expand = expansion(mult = c(0.12, 0.22)))
  }
  
  # Top headers: "Activated" and "Suppressed" written at the top of the chart!
  if (show_top_headers) {
    max_y <- length(levels(plot_df$TermWrapped))
    neg_x <- if (any(plot_df$NES < 0)) mean(c(min(plot_df$NES), 0)) else -0.8
    pos_x <- if (any(plot_df$NES > 0)) mean(c(max(plot_df$NES), 0)) else 0.8
    y_pos <- max_y + 0.95
    
    if (header_style == "banner") {
      p <- p +
        annotate("label", x = neg_x, y = y_pos, label = "◄ Suppressed",
                 fill = "#EEF2F7", color = "#2C3E50", fontface = "bold", size = 4.2, linewidth = 0.5) +
        annotate("label", x = pos_x, y = y_pos, label = "Activated ►",
                 fill = "#FDF0ED", color = "#C0392B", fontface = "bold", size = 4.2, linewidth = 0.5)
    } else {
      p <- p +
        annotate("text", x = neg_x, y = y_pos, label = "◄ Suppressed",
                 color = "#2C3E50", fontface = "bold", size = 4.4) +
        annotate("text", x = pos_x, y = y_pos, label = "Activated ►",
                 color = "#C0392B", fontface = "bold", size = 4.4)
    }
    
    p <- p +
      coord_cartesian(clip = "off") +
      scale_y_discrete(expand = expansion(mult = c(0.04, 0.16))) +
      theme(plot.margin = margin(t = 28, r = 20, b = 10, l = 10))
  }
  
  return(p)
}

# ------------------------------------------------------------------------------
# 3. MODULE 1: DIFFERENTIAL EXPRESSION & PROTEOMICS PLOTS
# ------------------------------------------------------------------------------

standardize_dea_df <- function(df) {
  df <- as.data.frame(df)
  cols <- colnames(df)
  
  # 1. Detect Gene Symbol column (prioritize Gene.names, Gene_names, Symbol, then Gene)
  gene_symbol_col <- NULL
  pref_gene_patterns <- c(
    "^gene[._]?names?$",
    "^genenames?$",
    "^symbol$",
    "^gene[._]?symbol$",
    "^gene$"
  )
  
  for (pat in pref_gene_patterns) {
    matches <- grep(pat, cols, ignore.case = TRUE, value = TRUE)
    for (m in matches) {
      vals <- as.character(df[[m]])
      non_empty <- vals[!is.na(vals) & vals != ""]
      # Exclude if it's strictly row indices (e.g. 1, 2, 3...)
      if (length(non_empty) > 0 && !all(grepl("^[0-9]+$", head(non_empty, 50)))) {
        gene_symbol_col <- m
        break
      }
    }
    if (!is.null(gene_symbol_col)) break
  }
  
  # 2. Detect Protein ID column
  protein_id_col <- NULL
  pref_prot_patterns <- c(
    "^protein[._]?ids?$",
    "^proteinids?$",
    "^original[._]?ids?$",
    "^uniprot[._]?ids?$",
    "^protein$"
  )
  for (pat in pref_prot_patterns) {
    matches <- grep(pat, cols, ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) {
      protein_id_col <- matches[1]
      break
    }
  }
  
  # Extract Gene and Protein_ID
  has_real_prot <- !is.null(protein_id_col)
  raw_prot <- if (has_real_prot) as.character(df[[protein_id_col]]) else ""
  raw_gene <- if (!is.null(gene_symbol_col)) as.character(df[[gene_symbol_col]]) else raw_prot
  
  # Clean gene names: if NA or empty or "NA", fallback to Protein ID if available, else Gene_X
  bad_gene <- is.na(raw_gene) | raw_gene == "" | raw_gene == "NA"
  if (any(bad_gene)) {
    if (has_real_prot) {
      raw_gene[bad_gene] <- raw_prot[bad_gene]
    } else {
      raw_gene[bad_gene] <- paste0("Gene_", which(bad_gene))
    }
  }
  
  if (!has_real_prot) {
    raw_prot <- raw_gene
  }
  
  df$Gene <- raw_gene
  df$Protein_ID <- raw_prot
  
  # Check if Protein_ID is fake / placeholder (like Protein_1, Protein 1, etc.)
  is_fake_prot <- is.na(raw_prot) | raw_prot == "" | grepl("^Protein[ _]?[0-9]+$", raw_prot, ignore.case = TRUE)
  
  # Only add Protein_ID in parentheses if:
  # 1. It is a real protein ID (not fake/placeholder)
  # 2. It differs from Gene
  # 3. Gene symbol is duplicated across rows (multiple isoforms/entries)
  dup_genes <- duplicated(raw_gene) | duplicated(raw_gene, fromLast = TRUE)
  use_paren <- !is_fake_prot & (raw_gene != raw_prot) & dup_genes
  
  df$Display_Name <- ifelse(use_paren,
                            paste0(df$Gene, " (", df$Protein_ID, ")"),
                            df$Gene)
  
  if (!"logFC" %in% colnames(df)) {
    fc_cols <- grep("^(logfc|log2fc|fc|log2_fold_change)$", colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(fc_cols) > 0) df$logFC <- as.numeric(df[[fc_cols[1]]])
    else stop("DEA table must contain 'logFC' or 'log2FC' column.")
  } else {
    df$logFC <- as.numeric(df$logFC)
  }
  
  if (!"adj.P.Val" %in% colnames(df)) {
    padj_cols <- grep("^(adj\\.p\\.val|adj_pval|padj|fdr|qval|sca\\.adj\\.pval)$", colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(padj_cols) > 0) df$adj.P.Val <- as.numeric(df[[padj_cols[1]]])
    else if ("P.Value" %in% colnames(df)) df$adj.P.Val <- as.numeric(df$P.Value)
    else df$adj.P.Val <- 0.05
  } else {
    df$adj.P.Val <- as.numeric(df$adj.P.Val)
  }
  
  if (!"P.Value" %in% colnames(df)) {
    pval_cols <- grep("^(p\\.value|pval|p_val|pvalue|sca\\.p\\.value)$", colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(pval_cols) > 0) df$P.Value <- as.numeric(df[[pval_cols[1]]])
    else df$P.Value <- df$adj.P.Val
  } else {
    df$P.Value <- as.numeric(df$P.Value)
  }
  
  return(df)
}

#' Plot Publication Volcano Plot
plot_volcano <- function(df,
                         p_val_type = "adj.P.Val",
                         fc_cutoff = 1.0,
                         p_cutoff = 0.05,
                         up_color = "#E64B35",
                         down_color = "#3C5488",
                         ns_color = "#B09C85",
                         point_size = 2,
                         point_alpha = 0.75,
                         top_n_labels = 10,
                         selected_genes = NULL,
                         label_type = "top_fc",
                         label_criterion = "fc",
                         label_distribution = "per_side",
                         only_sig_labels = TRUE,
                         label_size = 3.5,
                         symmetric_x = TRUE,
                         show_cutoffs = TRUE,
                         title = "Volcano Plot",
                         subtitle = NULL,
                         xlab = "log2(Fold Change)",
                         ylab = NULL,
                         theme_name = "classic",
                         base_size = 12) {
  data <- standardize_dea_df(df)
  
  p_col <- if (p_val_type == "P.Value") "P.Value" else "adj.P.Val"
  data$P_Use <- pmax(data[[p_col]], 1e-300)
  data$NegLog10P <- -log10(data$P_Use)
  
  if (is.null(ylab)) {
    ylab <- if (p_val_type == "P.Value") "-log10(p-value)" else "-log10(Adjusted p-value / FDR)"
  }
  
  data$Regulation <- "Not Significant"
  data$Regulation[data$logFC >= fc_cutoff & data$P_Use <= p_cutoff] <- "Enriched / Up"
  data$Regulation[data$logFC <= -fc_cutoff & data$P_Use <= p_cutoff] <- "Reduced / Down"
  data$Regulation <- factor(data$Regulation, levels = c("Enriched / Up", "Reduced / Down", "Not Significant"))
  data$Row_ID <- seq_len(nrow(data))
  data$Label <- ""
  
  if (label_type == "custom" && !is.null(selected_genes) && length(selected_genes) > 0) {
    match_rows <- which(data$Gene %in% selected_genes | data$Protein_ID %in% selected_genes | data$Display_Name %in% selected_genes)
    if (only_sig_labels) {
      match_rows <- match_rows[data$Regulation[match_rows] != "Not Significant"]
    }
    if (length(match_rows) > 0) {
      data$Label[match_rows] <- data$Gene[match_rows]
    }
  } else if (label_type %in% c("top_fc", "top_p", "top_combined", "top_n") && top_n_labels > 0) {
    # Determine ranking criterion
    criterion <- if (label_type == "top_fc") {
      "fc"
    } else if (label_type == "top_p") {
      "p_val"
    } else if (label_type == "top_combined") {
      "combined"
    } else {
      # label_type == "top_n"
      if (!is.null(label_criterion) && label_criterion %in% c("fc", "p_val", "combined")) {
        label_criterion
      } else {
        "fc"
      }
    }
    
    # Candidate pool
    candidate_data <- if (only_sig_labels) {
      sig_df <- data %>% filter(Regulation != "Not Significant")
      if (nrow(sig_df) > 0) sig_df else data
    } else {
      data
    }
    
    candidate_data$CombinedScore <- abs(candidate_data$logFC) * candidate_data$NegLog10P
    
    if (label_distribution == "overall") {
      if (criterion == "fc") {
        chosen_rows <- candidate_data %>% arrange(desc(abs(logFC))) %>% filter(!duplicated(Gene)) %>% head(top_n_labels)
      } else if (criterion == "p_val") {
        chosen_rows <- candidate_data %>% arrange(P_Use) %>% filter(!duplicated(Gene)) %>% head(top_n_labels)
      } else {
        chosen_rows <- candidate_data %>% arrange(desc(CombinedScore)) %>% filter(!duplicated(Gene)) %>% head(top_n_labels)
      }
    } else {
      # per_side
      up_pool <- if (only_sig_labels && any(candidate_data$Regulation == "Enriched / Up")) {
        candidate_data %>% filter(Regulation == "Enriched / Up")
      } else {
        candidate_data %>% filter(logFC > 0)
      }
      
      down_pool <- if (only_sig_labels && any(candidate_data$Regulation == "Reduced / Down")) {
        candidate_data %>% filter(Regulation == "Reduced / Down")
      } else {
        candidate_data %>% filter(logFC < 0)
      }
      
      if (criterion == "fc") {
        up_chosen <- up_pool %>% arrange(desc(logFC)) %>% filter(!duplicated(Gene)) %>% head(top_n_labels)
        down_chosen <- down_pool %>% arrange(logFC) %>% filter(!duplicated(Gene)) %>% head(top_n_labels)
      } else if (criterion == "p_val") {
        up_chosen <- up_pool %>% arrange(P_Use) %>% filter(!duplicated(Gene)) %>% head(top_n_labels)
        down_chosen <- down_pool %>% arrange(P_Use) %>% filter(!duplicated(Gene)) %>% head(top_n_labels)
      } else {
        up_chosen <- up_pool %>% arrange(desc(CombinedScore)) %>% filter(!duplicated(Gene)) %>% head(top_n_labels)
        down_chosen <- down_pool %>% arrange(desc(CombinedScore)) %>% filter(!duplicated(Gene)) %>% head(top_n_labels)
      }
      chosen_rows <- bind_rows(up_chosen, down_chosen)
    }
    
    if (nrow(chosen_rows) > 0) {
      data$Label[chosen_rows$Row_ID] <- chosen_rows$Gene
    }
  }
  
  color_map <- c(
    "Enriched / Up" = up_color,
    "Reduced / Down" = down_color,
    "Not Significant" = ns_color
  )
  
  p <- ggplot(data, aes(x = logFC, y = NegLog10P, color = Regulation)) +
    geom_point(alpha = point_alpha, size = point_size) +
    scale_color_manual(values = color_map) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab, color = "Regulation") +
    get_theme(theme_name, base_size = base_size)
  
  if (show_cutoffs) {
    p <- p +
      geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "grey40", linewidth = 0.5) +
      geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "grey40", linewidth = 0.5)
  }
  
  if (symmetric_x) {
    max_x <- max(abs(data$logFC), na.rm = TRUE) * 1.05
    p <- p + xlim(-max_x, max_x)
  }
  
  if (any(data$Label != "") && requireNamespace("ggrepel", quietly = TRUE)) {
    labeled_sub <- data %>% filter(Label != "")
    p <- p + ggrepel::geom_text_repel(
      data = labeled_sub,
      aes(label = Label),
      color = "black",
      size = label_size,
      fontface = "bold",
      box.padding = 0.4,
      point.padding = 0.3,
      segment.color = "grey50",
      segment.size = 0.3,
      max.overlaps = 50
    )
  }
  
  return(p)
}

#' Parse experimental design table (e.g. design.txt or design.csv)
parse_design_table <- function(path_or_df) {
  if (is.null(path_or_df)) return(NULL)
  des <- if (is.character(path_or_df) && file.exists(path_or_df)) {
    tryCatch({
      d <- read.table(path_or_df, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      if (ncol(d) < 2) d <- read.csv(path_or_df, header = TRUE, stringsAsFactors = FALSE)
      d
    }, error = function(e) NULL)
  } else if (is.data.frame(path_or_df)) {
    path_or_df
  } else {
    NULL
  }
  if (is.null(des) || ncol(des) < 2) return(NULL)
  
  cols_lower <- tolower(colnames(des))
  samp_col <- grep("^(sample|samples|samplename|id)$", cols_lower, value = TRUE)
  grp_col <- grep("^(group|groups|condition|conditions|treatment)$", cols_lower, value = TRUE)
  
  if (length(samp_col) == 0) samp_col <- colnames(des)[2]
  else samp_col <- colnames(des)[which(cols_lower == samp_col[1])]
  
  if (length(grp_col) == 0) grp_col <- colnames(des)[1]
  else grp_col <- colnames(des)[which(cols_lower == grp_col[1])]
  
  clean_name <- function(x) gsub("^[X._]+", "", x)
  data.frame(
    Sample = clean_name(as.character(des[[samp_col]])),
    Group = as.character(des[[grp_col]]),
    stringsAsFactors = FALSE
  )
}

#' Plot Protein Abundance Across Samples / Conditions
plot_protein_abundance <- function(df,
                                   gene_name,
                                   sample_cols = NULL,
                                   group_mapping = NULL,
                                   plot_type = "profile_samples",
                                   palette = "NPG",
                                   show_line = TRUE,
                                   show_points = TRUE,
                                   title = NULL,
                                   theme_name = "classic",
                                   base_size = 12) {
  data <- standardize_dea_df(df)
  
  # Search target protein by Display_Name, Gene, or Protein_ID
  target_row <- data %>% filter(
    tolower(Display_Name) == tolower(gene_name) |
    tolower(Gene) == tolower(gene_name) |
    tolower(Protein_ID) == tolower(gene_name)
  )
  
  if (nrow(target_row) == 0) {
    # Fallback to regex word boundary search
    target_row <- data %>% filter(
      grepl(paste0("^", gene_name, "($|[^a-zA-Z0-9])"), Gene, ignore.case = TRUE) |
      grepl(paste0("^", gene_name, "($|[^a-zA-Z0-9])"), Display_Name, ignore.case = TRUE) |
      grepl(paste0("^", gene_name, "($|[^a-zA-Z0-9])"), Protein_ID, ignore.case = TRUE)
    )
  }
  
  if (nrow(target_row) == 0) stop(paste("Gene / Protein ID '", gene_name, "' not found in data."))
  
  target_row <- target_row[1, ]
  disp_gene <- target_row$Gene
  disp_name <- target_row$Display_Name
  disp_prot <- target_row$Protein_ID
  
  if (is.null(sample_cols)) {
    anno_pat <- "^(gene|count|row_id|logfc|log2fc|fc|aveexpr|t|b|p[._]?val.*|adj[._]?p.*|padj|fdr|qval.*|sca[._].*|disp.*|prot.*|sig.*|comp.*)$"
    skip_idx <- grepl(anno_pat, colnames(df), ignore.case = TRUE)
    sample_cols <- names(df)[sapply(df, is.numeric) & !skip_idx]
  }
  
  if (length(sample_cols) == 0) stop("No numeric sample intensity columns found.")
  
  intensities <- as.numeric(target_row[sample_cols])
  plot_df <- data.frame(
    Sample = sample_cols,
    Intensity = intensities,
    stringsAsFactors = FALSE
  )
  
  clean_name <- function(x) gsub("^[X._]+", "", x)
  
  if (!is.null(group_mapping)) {
    if (is.data.frame(group_mapping)) {
      des_df <- group_mapping
      cols_l <- tolower(colnames(des_df))
      s_col <- grep("^(sample|samples|id)$", cols_l, value = TRUE)
      g_col <- grep("^(group|groups|condition|conditions)$", cols_l, value = TRUE)
      s_name <- if (length(s_col) > 0) colnames(des_df)[which(cols_l == s_col[1])] else colnames(des_df)[1]
      g_name <- if (length(g_col) > 0) colnames(des_df)[which(cols_l == g_col[1])] else colnames(des_df)[2]
      nm_clean <- clean_name(as.character(des_df[[s_name]]))
      grp_vals <- as.character(des_df[[g_name]])
      plot_df$Group <- sapply(clean_name(plot_df$Sample), function(s) {
        idx <- which(nm_clean == s)
        if (length(idx) > 0) grp_vals[idx[1]] else s
      })
    } else if (is.list(group_mapping)) {
      plot_df$Group <- "Unknown"
      for (grp in names(group_mapping)) {
        grp_members <- clean_name(group_mapping[[grp]])
        match_idx <- clean_name(plot_df$Sample) %in% grp_members
        plot_df$Group[match_idx] <- grp
      }
    } else if (is.character(group_mapping) && !is.null(names(group_mapping))) {
      nm_clean <- clean_name(names(group_mapping))
      plot_df$Group <- sapply(clean_name(plot_df$Sample), function(s) {
        idx <- which(nm_clean == s)
        if (length(idx) > 0) as.character(group_mapping[idx[1]]) else s
      })
    } else if (is.character(group_mapping) && length(group_mapping) == length(sample_cols)) {
      plot_df$Group <- group_mapping
    } else {
      plot_df$Group <- plot_df$Sample
    }
  } else {
    clean_grp <- sub("(_r?|\\.)[0-9]+$", "", plot_df$Sample)
    clean_grp <- gsub("^X[._]?", "", clean_grp)
    if (all(clean_grp == "") || length(unique(clean_grp)) == 1) {
      if (length(sample_cols) == 2 && "Comparison" %in% colnames(target_row) && !is.na(target_row$Comparison) && target_row$Comparison != "") {
        comp_parts <- strsplit(as.character(target_row$Comparison), "_vs_")[[1]]
        if (length(comp_parts) == 2) {
          plot_df$Group <- comp_parts
        } else {
          plot_df$Group <- plot_df$Sample
        }
      } else {
        plot_df$Group <- plot_df$Sample
      }
    } else {
      if (all(grepl("^[0-9]+$", clean_grp))) {
        clean_grp <- paste("Group", clean_grp)
      }
      plot_df$Group <- clean_grp
    }
  }
  
  plot_df$Group <- factor(plot_df$Group, levels = unique(plot_df$Group))
  
  # Ensure clean title without fake (Protein 1)
  if (is.null(title) || title == "") {
    title <- paste("Abundance Profile of", disp_name)
  }
  
  cols <- JOURNAL_PALETTES[[palette]]
  if (is.null(cols)) cols <- JOURNAL_PALETTES$NPG
  
  # Normalize plot_type: remove barplot and violin, map them to profile
  if (plot_type %in% c("bar_sem", "bar_sd", "barplot", "bar")) {
    plot_type <- "profile_group"
  } else if (plot_type == "violin") {
    plot_type <- "boxplot"
  }
  
  if (plot_type == "profile_samples") {
    # Order samples by group
    plot_df <- plot_df %>% arrange(Group, Sample)
    plot_df$Sample <- factor(plot_df$Sample, levels = unique(plot_df$Sample))
    
    p <- ggplot(plot_df, aes(x = Sample, y = Intensity))
    if (show_line) {
      p <- p + geom_line(aes(group = 1), color = "grey55", linewidth = 0.8)
    }
    p <- p +
      geom_point(aes(color = Group, fill = Group), size = 4, shape = 21, stroke = 1) +
      scale_color_manual(values = rep(cols, length.out = length(unique(plot_df$Group)))) +
      scale_fill_manual(values = rep(cols, length.out = length(unique(plot_df$Group)))) +
      labs(title = title, x = "Sample", y = "Log2 Intensity / Abundance", color = "Group", fill = "Group") +
      get_theme(theme_name, base_size = base_size) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "right"
      )
  } else if (plot_type == "profile_group") {
    summary_df <- plot_df %>%
      group_by(Group) %>%
      summarise(
        Mean = mean(Intensity, na.rm = TRUE),
        SEM = {
          s <- sd(Intensity, na.rm = TRUE)
          n <- sum(!is.na(Intensity))
          if (n > 1) s / sqrt(n) else 0
        },
        .groups = "drop"
      ) %>%
      mutate(
        ymin = Mean - SEM,
        ymax = Mean + SEM
      )
    summary_df$Group <- factor(summary_df$Group, levels = levels(plot_df$Group))
    summary_df <- summary_df %>% arrange(Group)
    
    p <- ggplot(summary_df, aes(x = Group, y = Mean))
    if (show_line) {
      p <- p + geom_line(aes(group = 1), color = "grey40", linewidth = 1.0)
    }
    p <- p + geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.18, linewidth = 0.8, color = "black")
    if (show_points) {
      p <- p + geom_jitter(data = plot_df, aes(x = Group, y = Intensity, fill = Group), width = 0.15, size = 2.5, shape = 21, color = "black", alpha = 0.65, inherit.aes = FALSE)
    }
    p <- p +
      geom_point(aes(fill = Group), size = 5, shape = 21, color = "black", stroke = 1.1) +
      scale_fill_manual(values = rep(cols, length.out = length(unique(plot_df$Group)))) +
      labs(title = title, x = "Condition / Group", y = "Log2 Abundance (Mean ± SEM)") +
      get_theme(theme_name, base_size = base_size) +
      theme(legend.position = "none")
  } else {
    # Default: Boxplot + Jitter
    p <- ggplot(plot_df, aes(x = Group, y = Intensity, fill = Group)) +
      geom_boxplot(alpha = 0.75, outlier.shape = NA, width = 0.45, color = "black", linewidth = 0.6) +
      geom_jitter(width = 0.15, size = 2.5, shape = 21, fill = "white", color = "black", alpha = 0.85) +
      scale_fill_manual(values = rep(cols, length.out = length(unique(plot_df$Group)))) +
      labs(title = title, x = "Condition / Group", y = "Log2 Abundance / Intensity") +
      get_theme(theme_name, base_size = base_size) +
      theme(legend.position = "none")
  }
  
  return(p)
}

#' Plot PCA (Principal Component Analysis)
plot_pca <- function(df,
                     sample_cols = NULL,
                     group_mapping = NULL,
                     palette = "NPG",
                     show_labels = TRUE,
                     show_ellipse = TRUE,
                     title = "PCA of Proteomics Samples",
                     theme_name = "classic",
                     base_size = 12) {
  if (is.null(sample_cols)) {
    anno_pat <- "^(gene|count|row_id|logfc|log2fc|fc|aveexpr|t|b|p[._]?val.*|adj[._]?p.*|padj|fdr|qval.*|sca[._].*|disp.*|prot.*|sig.*|comp.*)$"
    skip_idx <- grepl(anno_pat, colnames(df), ignore.case = TRUE)
    sample_cols <- names(df)[sapply(df, is.numeric) & !skip_idx]
  }
  if (length(sample_cols) < 3) stop("At least 3 samples are needed for PCA.")
  
  mat <- as.matrix(df[, sample_cols, drop = FALSE])
  mat <- mat[complete.cases(mat), , drop = FALSE]
  vars <- apply(mat, 1, var, na.rm = TRUE)
  mat <- mat[!is.na(vars) & vars > 0, , drop = FALSE]
  if (nrow(mat) < 2) stop("Not enough non-zero variance features for PCA.")
  
  pca_res <- prcomp(t(mat), scale. = TRUE)
  var_exp <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)
  
  pca_df <- data.frame(
    Sample = sample_cols,
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    stringsAsFactors = FALSE
  )
  
  if (!is.null(group_mapping) && is.list(group_mapping)) {
    pca_df$Group <- "Unknown"
    for (grp in names(group_mapping)) {
      pca_df$Group[pca_df$Sample %in% group_mapping[[grp]]] <- grp
    }
  } else {
    clean_grp <- gsub("[0-9r._]+$", "", pca_df$Sample)
    clean_grp <- gsub("^X[._]?", "", clean_grp)
    if (all(clean_grp == "") || length(unique(clean_grp)) == 1) {
      if ("Comparison" %in% colnames(df) && !is.na(df$Comparison[1]) && df$Comparison[1] != "") {
        comp_parts <- strsplit(as.character(df$Comparison[1]), "_vs_")[[1]]
        if (length(comp_parts) == 2) {
          half <- length(sample_cols) / 2
          pca_df$Group <- c(rep(comp_parts[1], half), rep(comp_parts[2], length(sample_cols) - half))
        } else {
          pca_df$Group <- pca_df$Sample
        }
      } else {
        pca_df$Group <- pca_df$Sample
      }
    } else {
      pca_df$Group <- clean_grp
    }
  }
  
  cols <- JOURNAL_PALETTES[[palette]]
  if (is.null(cols)) cols <- JOURNAL_PALETTES$NPG
  
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
    geom_point(size = 4, alpha = 0.9) +
    scale_color_manual(values = rep(cols, length.out = length(unique(pca_df$Group)))) +
    scale_fill_manual(values = rep(cols, length.out = length(unique(pca_df$Group)))) +
    labs(
      title = title,
      x = paste0("PC1 (", var_exp[1], "% variance)"),
      y = paste0("PC2 (", var_exp[2], "% variance)")
    ) +
    get_theme(theme_name, base_size = base_size)
  
  # Only draw ellipses if there are multiple groups AND each group has at least 3 samples
  grp_counts <- table(pca_df$Group)
  can_ellipse <- length(grp_counts) > 1 && all(grp_counts >= 3)
  if (show_ellipse && can_ellipse) {
    p <- p + stat_ellipse(geom = "polygon", alpha = 0.15, linetype = "dashed", level = 0.95)
  }
  
  if (show_labels && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      aes(label = Sample),
      color = "black",
      size = 3.5,
      box.padding = 0.3
    )
  }
  
  return(p)
}

#' Plot Expression Heatmap with Hierarchical Clustering
plot_expression_heatmap <- function(df,
                                    sample_cols = NULL,
                                    group_mapping = NULL,
                                    scale = "row",
                                    cluster_rows = TRUE,
                                    cluster_cols = FALSE,
                                    show_rownames = TRUE,
                                    show_colnames = TRUE,
                                    top_n = NULL,
                                    palette = "RdBu",
                                    title = NULL,
                                    fontsize = 10) {
  std_df <- standardize_dea_df(df)
  
  if (is.null(sample_cols)) {
    anno_pat <- "^(gene|count|row_id|logfc|log2fc|fc|aveexpr|t|b|p[._]?val.*|adj[._]?p.*|padj|fdr|qval.*|sca[._].*|disp.*|prot.*|sig.*|comp.*)$"
    skip_idx <- grepl(anno_pat, colnames(df), ignore.case = TRUE)
    sample_cols <- names(df)[sapply(df, is.numeric) & !skip_idx]
  }
  if (length(sample_cols) < 2) stop("At least 2 sample columns required for heatmap.")
  
  # Filter top_n if requested
  if (!is.null(top_n) && is.numeric(top_n) && top_n > 0 && top_n < nrow(std_df)) {
    if ("adj.P.Val" %in% colnames(std_df)) {
      std_df <- std_df[order(std_df$adj.P.Val, -abs(std_df$logFC)), ]
    } else {
      std_df <- std_df[order(-abs(std_df$logFC)), ]
    }
    std_df <- head(std_df, top_n)
  }
  
  mat <- as.matrix(std_df[, sample_cols, drop = FALSE])
  mat <- mat[complete.cases(mat), , drop = FALSE]
  
  row_labels <- std_df$Display_Name[seq_len(nrow(mat))]
  rownames(mat) <- make.unique(as.character(row_labels))
  
  clean_name <- function(x) gsub("^[X._]+", "", x)
  
  if (!is.null(group_mapping)) {
    if (is.data.frame(group_mapping)) {
      des_df <- group_mapping
      cols_l <- tolower(colnames(des_df))
      s_col <- grep("^(sample|samples|id)$", cols_l, value = TRUE)
      g_col <- grep("^(group|groups|condition|conditions)$", cols_l, value = TRUE)
      s_name <- if (length(s_col) > 0) colnames(des_df)[which(cols_l == s_col[1])] else colnames(des_df)[1]
      g_name <- if (length(g_col) > 0) colnames(des_df)[which(cols_l == g_col[1])] else colnames(des_df)[2]
      nm_clean <- clean_name(as.character(des_df[[s_name]]))
      grp_vals <- as.character(des_df[[g_name]])
      col_groups <- sapply(clean_name(sample_cols), function(s) {
        idx <- which(nm_clean == s)
        if (length(idx) > 0) grp_vals[idx[1]] else s
      })
    } else if (is.list(group_mapping)) {
      col_groups <- rep("Unknown", length(sample_cols))
      for (grp in names(group_mapping)) {
        grp_members <- clean_name(group_mapping[[grp]])
        match_idx <- clean_name(sample_cols) %in% grp_members
        col_groups[match_idx] <- grp
      }
    } else if (is.character(group_mapping) && !is.null(names(group_mapping))) {
      nm_clean <- clean_name(names(group_mapping))
      col_groups <- sapply(clean_name(sample_cols), function(s) {
        idx <- which(nm_clean == s)
        if (length(idx) > 0) as.character(group_mapping[idx[1]]) else s
      })
    } else if (is.character(group_mapping) && length(group_mapping) == length(sample_cols)) {
      col_groups <- group_mapping
    } else {
      col_groups <- sample_cols
    }
  } else {
    clean_grp <- sub("(_r?|\\.)[0-9]+$", "", sample_cols)
    clean_grp <- gsub("^X[._]?", "", clean_grp)
    if (all(clean_grp == "") || length(unique(clean_grp)) == 1) {
      if (length(sample_cols) == 2 && "Comparison" %in% colnames(std_df) && !is.na(std_df$Comparison[1]) && std_df$Comparison[1] != "") {
        comp_parts <- strsplit(as.character(std_df$Comparison[1]), "_vs_")[[1]]
        if (length(comp_parts) == 2) {
          col_groups <- comp_parts
        } else {
          col_groups <- sample_cols
        }
      } else {
        col_groups <- sample_cols
      }
    } else {
      if (all(grepl("^[0-9]+$", clean_grp))) {
        clean_grp <- paste("Group", clean_grp)
      }
      col_groups <- clean_grp
    }
  }
  
  unique_grps <- unique(col_groups)
  anno_col <- data.frame(
    Group = factor(col_groups, levels = unique_grps),
    row.names = sample_cols
  )
  
  # Coordinated publication colors for annotation groups
  npg_cols <- JOURNAL_PALETTES$NPG
  group_colors <- setNames(rep(npg_cols, length.out = length(unique_grps)), unique_grps)
  anno_colors <- list(Group = group_colors)
  
  # Colors
  if (palette == "viridis") {
    color_pal <- viridis::viridis(100)
  } else if (palette == "plasma") {
    color_pal <- viridis::plasma(100)
  } else if (palette == "magma") {
    color_pal <- viridis::magma(100)
  } else if (palette == "bwr") {
    color_pal <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
  } else {
    # Default: RdBu / NPG style blue-white-red
    color_pal <- colorRampPalette(c("#3C5488", "white", "#E64B35"))(100)
  }
  
  main_title <- if (!is.null(title) && title != "") {
    title
  } else if ("Comparison" %in% colnames(std_df) && !is.na(std_df$Comparison[1])) {
    paste("Expression Heatmap -", std_df$Comparison[1])
  } else {
    "Expression Heatmap"
  }
  
  # Remove 0 variance rows before scaling
  vars <- apply(mat, 1, var, na.rm = TRUE)
  mat <- mat[!is.na(vars) & vars > 0, , drop = FALSE]
  if (nrow(mat) == 0) stop("No rows with non-zero variance for heatmap.")
  
  ph <- pheatmap::pheatmap(
    mat,
    scale = scale,
    annotation_col = anno_col,
    annotation_colors = anno_colors,
    cluster_rows = cluster_rows && nrow(mat) > 1,
    cluster_cols = cluster_cols && ncol(mat) > 1,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    color = color_pal,
    main = main_title,
    fontsize = fontsize,
    border_color = NA,
    silent = TRUE
  )
  
  return(ph$gtable)
}

# ------------------------------------------------------------------------------
# 4. MODULE 3: TRANSCRIPTION FACTOR FOOTPRINTING PLOTS
# ------------------------------------------------------------------------------

standardize_tf_df <- function(df) {
  df <- as.data.frame(df)
  
  if (!"source" %in% colnames(df)) {
    tf_cols <- grep("^(source|tf|transcription_factor|gene)$", colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(tf_cols) > 0) df$source <- as.character(df[[tf_cols[1]]])
    else stop("TF table must contain a 'source' or 'TF' column.")
  } else {
    df$source <- as.character(df$source)
  }
  
  if (!"score" %in% colnames(df)) {
    score_cols <- grep("^(score|activity|nes|t_stat)$", colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(score_cols) > 0) df$score <- as.numeric(df[[score_cols[1]]])
    else stop("TF table must contain a 'score' column.")
  } else {
    df$score <- as.numeric(df$score)
  }
  
  if (!"p_value" %in% colnames(df)) {
    pval_cols <- grep("^(p_value|pvalue|p.value|padj)$", colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(pval_cols) > 0) df$p_value <- as.numeric(df[[pval_cols[1]]])
    else df$p_value <- 0.05
  } else {
    df$p_value <- as.numeric(df$p_value)
  }
  
  if (!"Status" %in% colnames(df)) {
    df$Status <- ifelse(df$score > 0, "Activated", "Repressed")
  }
  
  return(df)
}

#' Plot TF Activity Diverging Barplot
plot_tf_barplot <- function(df,
                            top_n = 15,
                            selected_tfs = NULL,
                            filter_mode = "exclude",
                            act_color = "#E64B35",
                            rep_color = "#3C5488",
                            orientation = "horizontal",
                            title = "Inferred Transcription Factor Activities",
                            subtitle = NULL,
                            xlab = "TF Activity Score",
                            ylab = "Transcription Factor",
                            theme_name = "classic",
                            base_size = 12) {
  data <- standardize_tf_df(df)
  
  if (!is.null(selected_tfs) && length(selected_tfs) > 0) {
    if (filter_mode == "exclude") {
      data <- data %>% filter(!source %in% selected_tfs)
    } else {
      data <- data %>% filter(source %in% selected_tfs)
    }
  }
  
  act <- data %>% filter(score > 0) %>% arrange(desc(score)) %>% head(top_n)
  rep <- data %>% filter(score < 0) %>% arrange(score) %>% head(top_n)
  plot_df <- bind_rows(act, rep) %>% arrange(score)
  
  plot_df$source <- factor(plot_df$source, levels = plot_df$source)
  
  p <- ggplot(plot_df, aes(x = score, y = source, fill = Status)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.6) +
    scale_fill_manual(values = c("Activated" = act_color, "Repressed" = rep_color)) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab, fill = "TF Status") +
    get_theme(theme_name, base_size = base_size)
  
  if (orientation == "vertical") {
    p <- p + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  return(p)
}

#' Plot TF Activity Lollipop
plot_tf_lollipop <- function(df,
                             top_n = 15,
                             selected_tfs = NULL,
                             filter_mode = "exclude",
                             act_color = "#E64B35",
                             rep_color = "#3C5488",
                             point_size = 5,
                             title = "TF Activity Lollipop Plot",
                             subtitle = NULL,
                             theme_name = "classic",
                             base_size = 12) {
  data <- standardize_tf_df(df)
  
  if (!is.null(selected_tfs) && length(selected_tfs) > 0) {
    if (filter_mode == "exclude") {
      data <- data %>% filter(!source %in% selected_tfs)
    } else {
      data <- data %>% filter(source %in% selected_tfs)
    }
  }
  
  act <- data %>% filter(score > 0) %>% arrange(desc(score)) %>% head(top_n)
  rep <- data %>% filter(score < 0) %>% arrange(score) %>% head(top_n)
  plot_df <- bind_rows(act, rep) %>% arrange(score)
  plot_df$source <- factor(plot_df$source, levels = plot_df$source)
  
  p <- ggplot(plot_df, aes(x = score, y = source, color = Status)) +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey40") +
    geom_segment(aes(x = 0, xend = score, y = source, yend = source), linewidth = 0.8) +
    geom_point(size = point_size) +
    scale_color_manual(values = c("Activated" = act_color, "Repressed" = rep_color)) +
    labs(title = title, subtitle = subtitle, x = "TF Activity Score", y = "Transcription Factor", color = "Status") +
    get_theme(theme_name, base_size = base_size)
  
  return(p)
}

# ------------------------------------------------------------------------------
# 5. UNIVERSAL HIGH-RES EXPORT HELPER
# ------------------------------------------------------------------------------

export_figure <- function(plot_obj,
                          file_path,
                          format = "png",
                          width = 8,
                          height = 6,
                          units = "in",
                          dpi = 300) {
  if (inherits(plot_obj, "pheatmap")) {
    plot_obj <- plot_obj$gtable
  }
  
  ext <- tolower(format)
  if (!grepl(paste0("\\.", ext, "$"), file_path)) {
    file_path <- paste0(file_path, ".", ext)
  }
  
  if (ext %in% c("png", "tiff", "tif", "jpeg", "jpg")) {
    ggsave(filename = file_path, plot = plot_obj, width = width, height = height, units = units, dpi = dpi)
  } else if (ext == "pdf") {
    ggsave(filename = file_path, plot = plot_obj, width = width, height = height, units = units, device = cairo_pdf)
  } else if (ext == "svg") {
    if (requireNamespace("svglite", quietly = TRUE)) {
      ggsave(filename = file_path, plot = plot_obj, width = width, height = height, units = units, device = "svg")
    } else {
      ggsave(filename = file_path, plot = plot_obj, width = width, height = height, units = units)
    }
  } else {
    ggsave(filename = file_path, plot = plot_obj, width = width, height = height, units = units, dpi = dpi)
  }
  
  return(file_path)
}
