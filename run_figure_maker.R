# ==============================================================================
# ProteoPrism Launcher
# Double click or run in RStudio / Console to launch the application
# ==============================================================================

# Ensure required packages are present
required_pkgs <- c("shiny", "bslib", "DT", "ggplot2", "dplyr", "tidyr", "stringr", "scales")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

# Optional high-res / palette packages
optional_pkgs <- c("ggrepel", "ggsci", "viridis", "RColorBrewer", "svglite")
opt_missing <- optional_pkgs[!sapply(optional_pkgs, requireNamespace, quietly = TRUE)]
if (length(opt_missing) > 0) {
  message("Installing optional visualization packages: ", paste(opt_missing, collapse = ", "))
  install.packages(opt_missing, repos = "https://cloud.r-project.org")
}

# Determine app directory
app_dir <- tryCatch({
  # If executed via RStudio or source()
  dirname(normalizePath(sys.frame(1)$ofile))
}, error = function(e) {
  # Fallback locations
  locs <- c(
    "c:/Users/serife-ayaz/Desktop/Differential_Proteomics_Pipeline/Figure_Maker",
    file.path(getwd(), "Figure_Maker"),
    getwd()
  )
  for (l in locs) {
    if (file.exists(file.path(l, "app.R"))) return(l)
  }
  return(getwd())
})

message("Starting ProteoPrism from: ", app_dir)
shiny::runApp(app_dir, launch.browser = TRUE)
