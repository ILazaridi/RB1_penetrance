
library(tidyverse)
library(ggplot2)
library(cowplot)


# ---------------------------
# Data: prevalence counts
# ---------------------------

# Creates a wide-format data frame:
# - columns are cohorts (UKB, AoU, Combined)
# - rows are groups (observed, expected upper bound, expected lower bound)

df_2 <- data.frame("UKB" = c(21, 44, 5),
                   "AoU" = c(28, 29, 3),
                   "Combined" = c(49, 73, 8),
                   row.names = c("observed", "expected_upper", "expected_lower"), 
                   stringsAsFactors = FALSE)

df_2  # prints to console

# ---------------------------
# Reshape to long format for ggplot
# ---------------------------

# Converts row names (observed/expected_*) into a real column called "group"
# Then pivots so you get one row per (group, cohort) with a "count"

df_long_2 <- df_2 %>%
  rownames_to_column("group") %>%
  pivot_longer(
    cols = -group,
    names_to = "cohort",
    values_to = "count")

# Sets the x-axis order (otherwise ggplot might sort alphabetically)

df_long_2$cohort <- factor(df_long_2$cohort, levels = c("UKB", "AoU", "Combined"))

df_long_2 # prints long data


# Split into two datasets:
# - df_obs: the bars
# - df_exp: the horizontal expected lines (upper/lower)

df_obs <- df_long_2 %>% filter(group == "observed")
df_exp <- df_long_2 %>% filter(group != "observed")

df_obs
df_exp

# ---------------------------
# Plot parameters
# ---------------------------

# Controls bar width and the "half-width" used to draw the expected line segments
# so the segment spans the width of the bar

bar_width <- 0.3
half_w <- bar_width / 2

# ---------------------------
# Build plot
# ---------------------------

prevalence_graph_2 <- ggplot() +
  
  # ---- Observed bars ----
  geom_col(
    data = df_obs,
    aes(x = cohort, y = count, fill = group),
    width = bar_width,
    fill = 'lightgray',
    show.legend = FALSE) +
  
  # Label observed bars with sample size text
  # vjust > 0 moves text downward (relative to y position)
  geom_text(
    data = df_obs,
    aes(x = cohort, y = count, label = paste0("(n=", count, ")")),
    vjust = 1.5,   # moves labels just below the bars
    size = 8,
    fontface = "bold", 
    color = "black") +
  
  # ---- Expected prevalence lines ----
  # Uses numeric x positions so you can draw a horizontal segment centered on each bar.
  # x = cohort index - half_w, xend = cohort index + half_w
  geom_segment(
    data = df_exp,
    aes(x = as.numeric(factor(cohort)) - half_w,
        xend = as.numeric(factor(cohort)) + half_w,
        y = count,
        yend = count,
        color = group),
    linewidth = 1.5,
    inherit.aes = FALSE) +
  
  # Labels for expected lines (nudged right of the segment)
  geom_text(
    data = df_exp,
    aes(x = as.numeric(factor(cohort)) + 0.20,   # nudges label to the right
        y = count,
        label = ifelse(group == "expected_upper", "Expected- upper", "Expected- lower")),
    inherit.aes = FALSE,
    hjust = 0,
    size = 8,
    color = "black") +
  
  # Axis labels (x is blank on purpose)
  labs(x = NULL, y = "Prevalence",
       color = "Expected",
       size = 15) +
  
  # Both expected_* groups are set to black, so they’re visually identical.
  # (You’re distinguishing them via their text labels instead)
  
  scale_color_manual(
    values = c("expected_upper" = "black",
               "expected_lower" = "black")) +
  
  # Adds padding so right-side labels don’t get clipped
  scale_x_discrete(expand = expansion(mult = c(0.18, 0.3))) +
  
  # Sets y ticks every 5 units from 0 to 80
  scale_y_continuous(breaks = seq(0, 80, by = 5)) +
  
  scale_fill_viridis_d(option = "viridis") +
  
  # Plot title
  labs(title = "Prevalence of retinoblastoma in the UK Biobank (UKB) and All of Us (AoU)", fill = NULL, size = 15) +
  
  # Set theme
  theme_cowplot() +
  
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white"),
    
    # Typography sizing and spacing
    axis.title.y = element_text(margin = margin(r = 15), size = 25),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(margin = margin(t = 10), size = 25),
    plot.title = element_text(hjust = 0.5, size=25))

prevalence_graph_2

# ---------------------------
# Save output
# ---------------------------

# Saves to your current working directory

ggsave("prevalence_plot_UKB_AOU_combined.png", plot = prevalence_graph_2, width = 20, height = 12, dpi = 600)

