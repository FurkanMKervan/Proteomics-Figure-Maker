# ProteoFigure Studio - Publication-Grade Figure Maker

**ProteoFigure Studio** is a visualization and publication-quality figure generation tool designed for differential proteomics data analysis.

It seamlessly integrates outputs from proteomics analysis workflows (`Proteomics_Pipeline.Rmd`, `Ultimate_Functional_Analysis.Rmd`, and `Proteomic_Footprinting.Rmd`), providing diverse plotting options, academic journal color palettes, and full customization capabilities.

---

## Key Features

- **Functional Enrichment Plots:** Dot Plot, Horizontal/Vertical Barplot, Lollipop Chart, and GSEA Waterfall plot.
- **Differential Expression Visualizations:** Enhanced Volcano Plot (custom gene labeling with `ggrepel`), Individual Protein Boxplots / Violin Plots, and PCA with confidence ellipses.
- **Transcription Factor (TF) Dynamics:** Diverging Barplots, TF Lollipops, and activity score distributions.
- **Journal-Grade Palettes:** Nature Publishing Group (NPG), NEJM, The Lancet, JAMA, JCO, Science (AAAS), Viridis, ColorBrewer, and custom HEX color pickers.
- **Publication-Ready Exports:** Vector formats (PDF, SVG for Adobe Illustrator), high-resolution raster images (300 & 600 DPI PNG, TIFF), with Nature single-column (89 mm) and double-column (183 mm) dimension presets.

---

## Usage

### 1. Interactive Interface (Shiny App) - Recommended

Launch the interactive application using RStudio or the R console:

```r
source("run_figure_maker.R")
```

Or run `app.R` directly:
```r
shiny::runApp()
```

Workflow:
1. Select your results directory or upload analysis tables.
2. Choose your plot type, metrics (GeneRatio, Count, FoldEnrichment), and color palette.
3. Inspect the live preview and adjust dimensions.
4. Export publication-ready figures.

### 2. Programmatic Usage in RMarkdown / Scripts

Call functions from `figure_maker_engine.R` directly in your R workflows:

```r
source("figure_maker_engine.R")

# 1. Enrichment Dot Plot
df_go <- read.csv("Ultimate_GO_BP_Enriched_Young_OMV_vs_Young_Control.csv")
p_dot <- plot_enrichment_dot(df_go, x_metric = "GeneRatio", color_metric = "p.adjust", top_n = 15)

# 2. Volcano Plot with Custom Gene Labels
df_dea <- read.csv("all_proteins_Young_OMV_vs_Young_Control.csv")
p_volc <- plot_volcano(df_dea, selected_genes = c("IL6", "C3", "CASP3"), label_type = "custom")

# 3. Export High-Resolution Figure
export_figure(p_dot, file_path = "Figure1.pdf", format = "pdf", width = 7.2, height = 6.0)
```

For detailed examples and tutorials, see [Figure_Maker_Tutorial.Rmd](Figure_Maker_Tutorial.Rmd).

---

## 📁 Repository Structure

```text
Proteomics Figure Maker/
├── app.R                     # Interactive Shiny Application
├── figure_maker_engine.R     # Core modular plotting functions engine
├── run_figure_maker.R        # Application launcher script
├── Figure_Maker_Tutorial.Rmd # Step-by-step tutorial notebook
├── LICENSE                   # MIT License
└── README.md                 # Project documentation
```

---

## Prerequisites

Required packages:
`shiny`, `bslib`, `DT`, `ggplot2`, `dplyr`, `tidyr`, `stringr`, `scales`, `ggrepel`, `ggsci`, `viridis`, `RColorBrewer`, `svglite`

Missing packages will be installed automatically when executing `run_figure_maker.R`.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
