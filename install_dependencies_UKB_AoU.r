# install_dependencies_UKB_AoU.r
# Run this ONCE to install all required R packages for the RB1 analysis in the UK Biobank and All of Us cohorts and generating associated plots

required_packages <- c(
  # Core data wrangling / plotting
  "tidyverse",
  "ggplot2",   
  "readxl",
  
  # Survival analysis
  "survival",
  "survminer",
  
  # Plot styling / extensions
  "ggthemes",
  "ggsci",
  "ggpubr",
  "ggbreak",
  "cowplot",
  
  # Development tools
  "devtools"
)

installed <- rownames(installed.packages())
to_install <- setdiff(required_packages, installed)

if (length(to_install) > 0) {
  install.packages(to_install, dependencies = TRUE)
}

message("All required packages are installed.")