
library(tidyverse)
library(ggplot2)
library(cowplot)

# -------------------------------------------------------------------
# 1) Load cohort
# -------------------------------------------------------------------

final_cohort <- read.delim('RB1_analysis_total_cohort.tsv', stringsAsFactors = FALSE)

# -------------------------------------------------------------------
# 2) Restrict to carriers and build derived phenotype-by-age variables
# -------------------------------------------------------------------

# Filter to carriers and keep only the columns used downstream
carriers <- final_cohort %>%
  filter(qcd_carrier_strict == 1) %>%
  select('qcd_carrier_strict', 'Rb_phenotype', 'age_diag_max')

# Create a binary variable: diagnosed by age 60 (1) vs after 60 / not by 60 (0)
carriers <- carriers %>%
  mutate(Rb_phenotype_by_60 = ifelse(age_diag_max <= 60, 1, 0))

carriers$Rb_phenotype_by_60[is.na(carriers$Rb_phenotype_by_60)]=0

# -------------------------------------------------------------------
# 3) Penetrance to age 60 as a pie chart
# -------------------------------------------------------------------

pie_chart <- carriers %>% 
  group_by(Rb_phenotype_by_60) %>% 
  summarise(carriers_sum = sum(qcd_carrier_strict), .groups = "drop") %>% 
  mutate(penetrance = carriers_sum/sum(carriers_sum) * 100, label_color = ifelse(Rb_phenotype_by_60 == 1, "white", "black")) %>% 
  
  
  ggplot(aes(x = "", y = penetrance, fill = factor(Rb_phenotype_by_60))) + 
  geom_col(width = 1, color = "black", linewidth = 0.6) + 
  coord_polar(theta ="y") + 
  geom_text(aes( label = paste0(carriers_sum, " (", round(penetrance, 1)," %)"), color = label_color),
            position = position_stack(vjust = 0.5), size = 17) + 
  scale_color_identity() +
  labs(title = "", fill = NULL, size = 12) + 
  scale_fill_manual(values = c("0" ="#9DE0A6", "1" = "#A21702"), 
                    labels = c( "0" = "Do not have Rb/Eye cancer", "1" = "Have Rb/Eye cancer")) + 
  
  theme_void() + # removes axes, ticks, background grid
  
  theme(plot.title = element_text(hjust = 0, vjust = 2.6, size = 14),
        legend.text = element_text(size = 20, margin = margin(5,5,5,5, "pt")),
        
        # Put legend inside plot area
        legend.position = "inside",
        legend.position.inside = c(0.24, 0.18),
        legend.justification = c("right", "top"),
        
        # Style legend box
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
        legend.key = element_rect(fill = "black", color = "black"),
        
        # Ensure white background for saved output
        panel.background = element_rect(fill = "white", colour = NA, linewidth = 0),
        plot.background  = element_rect(fill = "white", colour = NA, linewidth = 0)) # semi-transparent legend box


pie_chart

pie_chart_2 <- pie_chart +
  theme(legend.position = "none")


ggsave("penetrance_plot.png", plot = pie_chart, width = 16, height = 11, dpi = 600)
ggsave("penetrance_plot_no_legend.png", plot = pie_chart_2, width = 15, height = 10, dpi = 600)

