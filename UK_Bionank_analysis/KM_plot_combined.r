
library(ggpubr)
library(survival)
library(survminer)
library(tidyverse)
library(ggpubr)
library(survival)
library(survminer)
library(dplyr)
library(ggthemes)
library(ggsci)
library(stringr)
library(ggtext)

# ------------------------------------------------------------
# 1) Read data
# ------------------------------------------------------------

df_age_diag <- read.delim('all_data_with_diag_age_stringent_permissive_wgs_only.tsv', stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 2) Select the specific columns needed for plotting
# ------------------------------------------------------------

df_plot_1 <- df_age_diag %>%
  select(all_of(c('eid', 'recruit_age', 'dob', 'stringent_carrier', 
                  'stringent_final', 'stringent_agediag',
                  'permissive_final', 'permissive_agediag')))

# ------------------------------------------------------------
# 3) Define censoring/reference time
# ------------------------------------------------------------

# Reference year is the date of last HES/cancer registry update.
# Anyone without an event by this date is effectively censored here.
reference_year <- ymd('2023-07-01')
df_age_diag$dob <- as.Date(df_age_diag$dob)

# Ensure DOB is a Date (important for lubridate interval math)
df_plot_1 <- df_plot_1 %>%
  mutate(reference_age = floor(interval(dob, reference_year)/years(1)))

# ------------------------------------------------------------
# 4) Reorder columns (purely for readability / inspection)
# ------------------------------------------------------------

df_plot_1 <- df_plot_1 %>%
  relocate('dob', .after = 'eid') %>%
  relocate('stringent_agediag', .after = 'recruit_age') %>%
  relocate('stringent_final', .after = 'stringent_agediag') %>%
  relocate('permissive_agediag', .after = 'stringent_final') %>%
  relocate('permissive_final', .after = 'permissive_agediag') %>%
  relocate('reference_age', .after = 'permissive_final')


# ------------------------------------------------------------
# 5) Replace missing diagnosis ages (NA) with censoring age
# ------------------------------------------------------------

# Replace NA values in Age_at_Diagnosis with the last observed age before censoring
# Interpretation: if Age_at_Diagnosis is NA, person did NOT have the event
# observed; set their "time" to their last observed follow-up age (censoring).
df_plot_1$stringent_agediag[is.na(df_plot_1$stringent_agediag)] <- df_plot_1$reference_age[is.na(df_plot_1$stringent_agediag)]
df_plot_1$permissive_agediag[is.na(df_plot_1$permissive_agediag)] <- df_plot_1$reference_age[is.na(df_plot_1$permissive_agediag)]

# ------------------------------------------------------------
# 6) Build a long-format dataset (one row per person per phenotype)
# ------------------------------------------------------------

df_long <- bind_rows(
  # Outcome 1: Rb (retinoblastoma)
  df_plot_1 %>%
    transmute(ID = eid,
              time = stringent_agediag,
              event = stringent_final,
              phenotype = "Rb",
              carrier_status   = if_else(stringent_carrier == 1, "RB1+", "RB1-")),
  # Outcome 2: broader cancer phenotype
  df_plot_1 %>%
    transmute(ID = eid,
              time = permissive_agediag,
              event = permissive_final,
              phenotype = "Ocular & Extraocular cancers",
              carrier_status   = if_else(stringent_carrier == 1, "RB1+", "RB1-")))

# Set factor ordering (controls legend order and strata labeling order)
df_long <- df_long %>%
  mutate(
    phenotype = factor(phenotype, levels = c("Rb", "Ocular & Extraocular cancers")),
    carrier_status = factor(carrier_status, levels = c("RB1+", "RB1-"))) #0 = RB1- and 1 = RB1+

# ------------------------------------------------------------
# 7) Fit Kaplan–Meier curves
# ------------------------------------------------------------

fit <- survfit(Surv(time, event) ~ carrier_status + phenotype, data = df_long)

# ------------------------------------------------------------
# 8) Plot survival curves with survminer
# ------------------------------------------------------------

plot_1 <- ggsurvplot(
  fit, data=df_long,
  color = "phenotype", # map phenotype to color
  linetype = "carrier_status", # map carrier status to line type
  palette = c("blue", "red", "blue", "red"),
  
  xlim = c(0, 60),
  xlab = "Age (years)",
  ylab = 'Cancer free survival',
  censor = FALSE, # don't show censor tick marks
  conf.int = FALSE,
  break.time.by = 10,
  
  ggtheme = theme_minimal() +
    theme(legend.direction = "vertical",
          legend.title = element_text(size = 26, face = "bold",margin = margin(b = 18, r = 25, unit = "pt")),
          legend.text = element_text(size = 26, face = "italic", margin = margin(l = 25, unit = "pt")),
          legend.key.width = unit(2, 'cm'),
          legend.spacing.x = unit(25, "pt"),
          legend.box.spacing = unit(10, "pt"),
          legend.box.background = element_rect(fill = "white", colour = "gray20"),
          legend.box.margin = margin(2, 42, 2, 20, unit = "pt"),
          axis.title.y = element_text(size = 26, face = "bold", margin = margin(r = 25)),
          axis.title.x = element_text(size = 26, face = "bold" , margin = margin(t = 20)),
          axis.text.y = element_text(size = 26, color = "black"),
          axis.text.x = element_text(size = 26, color = "black")),
  
  risk.table = TRUE,
  risk.table.height = 0.18,
  risk.table.fontsize = 7,
  risk.table.width = 0.55,
  legend.title = "phenotype"  # legend title for the color aesthetic
  )


# ------------------------------------------------------------
# 9) Customize the risk table separately
# ------------------------------------------------------------

risk_table_narrow <- plot_1$table +
  scale_y_discrete(labels =  c(
    "RB1- \nOcular & Extraocular cancers",
    "RB1- Rb",
    "RB1+ \nOcular & Extraocular cancers",
    "RB1+ Rb")) +
  
  theme(plot.margin = margin(35, 2, 2, 2),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(color = "grey60", linewidth = 0.3),
        panel.background = element_rect(fill = "grey95", color = NA),
        
        plot.title = element_text(size = 25, face = "bold", hjust = 0),
        
        axis.text.y = element_markdown(size = 20, lineheight = 1, margin = margin(r = 30, unit = "pt")),
        
        axis.title.y = element_markdown(size = 25, margin = margin(r = 70)),
        
        axix.text.x = element_markdown(size = 20,  margin = margin(t = 20)),
        
        axis.title.x = element_markdown(color = "black", size = 25, margin = margin(t = 20)))


# ------------------------------------------------------------
# 10) Combine plots: survival curve alone, and curve + centered risk table
# ------------------------------------------------------------

combined_plot <- ggarrange(plot_1$plot)

#combined_plot_2 <- ggarrange(plot_1$plot, plot_1$table, ncol = 1, heights = c(2, 1))

combined_plot_2 <- ggarrange(
  plot_1$plot,
  
  cowplot::ggdraw() +
    cowplot::draw_plot(risk_table_narrow, x = 0.2, width = 0.8, height = 1),
  
  ncol = 1,
  heights = c(3, 1))

# ------------------------------------------------------------
# 11) Save to disk
# ------------------------------------------------------------

ggsave("KM_plot_combined_no_rt.png", plot = combined_plot, width = 28, height = 16, dpi = 600)
ggsave("KM_plot_combined_with_rt.png", plot = combined_plot_2, width = 22, height = 20, dpi = 600)
