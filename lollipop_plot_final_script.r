
# Required packages
library(ggplot2)
library(tidyverse)
library(ggbreak)
library(readxl)



# Import datasets

## ---- ClinVar and UKB/AoU variant data ----
# Excel file in the format (potentially change to tsv file input)?? 
# - Tag: ClinVar or UKB/AoU (distinguish between different variant sources)
# - genomic_position: genomic coordinates of variant (build 38),
# - clinical_significance
# - y-axis: (+) values for clinvar assertion level and (-) values for UKB/AoU allele count


lollipop_input <- read_excel("lollipop_plot_data.xlsx", sheet="Sheet6") %>%
  mutate(y_axis = as.numeric(y_axis), 
         genomic_position = as.numeric(genomic_position))


## ---- RB1 exon coordinate boundaries ----
# Columns: "type", "start", "end" and "ID"

#Specifying coding regions of exon 1 only (not including 5' UTR region)
# Will replace exon1 from GFF file which contains UTR
row_exon_1 <- data.frame(type="exon",
                         start=48303913,
                         end=48304049,
                         ID="exon-NM_000321.3-1",
                         exon_num=1)

#Specifying coding regions of exon 27 only (not including 3' UTR region)
# Will replace exon1 from GFF file which contains UTR
row_exon_27 <- data.frame(type="exon",
                          start=48479998,
                          end=48480071,
                          ID="exon-NM_000321.3-27",
                          exon_num=27)


exon_coordinates <- read_tsv("RB1_exon_coordinates.tsv") %>%
  select(type, start, end, ID) %>%
  slice(-c(1, 2, 28)) %>% # remove rows 1, 2 and 28 that correspond to full mRNA coordinates, exon1 coordinates and exon 27 coordinates respectively
  mutate(exon_num = sub(".*exon-NM_000321\\.3-", "", ID),
         exon_num = as.numeric(exon_num)) %>%
  bind_rows(row_exon_1, row_exon_27) %>% # replace values with previously specified coordinates for exon1/2
  arrange(exon_num)


# Plot asthetics

# Specify clinical classification colour as in the input file
classification_colors <- c("Pathogenic/Likely Pathogenic" = "red",  "Unreported" = "gray80")

lollipop_input$clinical_significance <- factor(lollipop_input$clinical_significance,
                                               levels = names(classification_colors))

classification_fill <- classification_colors

# Shape of point for each variant based on whether it was from the ClinVar dataset or the variants from UKB/AoU
classification_shapes <- c("ClinVar" = 24, "UKB/AoU" = 21)


## ---- Axis limits used in the plot ----

# Specify start and end of the gene (this region includes the UTRs)
gene_start <- 48303751
gene_end <- 48481890

# NMD escape region genomic coordinates (5' and 3') 
# +150bp from start of coding region (exon1)
# -55bp of penultimate exon (exon26) and all of final exon (exon27)
nmd_5prime <- c(start = 48303913, end = 48304063)
nmd_3prime <- c(start = 48477349, end = 48480071)


# Will be used to plot the exon boxes on the image
# Exon start/end coordinates to plot on x axis (specifies the length of the exon) and ymin/max to specify height of boxes
df_rect <- data.frame(xmin=exon_coordinates$start,
                      xmax=exon_coordinates$end,
                      ymin=-0.04,
                      ymax=0.04)




# Base plot
p1 <- ggplot(lollipop_input, 
             aes(y = y_axis, 
                 x = genomic_position, 
                 label ="", 
                 fill= clinical_significance, 
                 shape=Tag, 
                 colour = clinical_significance)) +
  
  # Lolipop bar size (vertical lines starting at baseline and ending at the corresponding y position)
  geom_segment(aes(x = genomic_position, 
                   y = 0, 
                   xend = genomic_position, 
                   yend = y_axis), 
               color = "black", 
               linewidth = 0.5)  +
  
  # Lollipop point size
  geom_point(size=2) + 
  
  # Y axis: baseline set at 0
  geom_hline(yintercept = 0, 
             size = 0.5, 
             color = "black") +
  
  theme_minimal() + 
  
  # Custom y axis specifying ClinVar submission level and UKB/AoU allele counts
  scale_y_continuous(limits = c(-1.5, 1.5), 
                     breaks = seq(-1.5, 1.5, 0.5),
                     labels = c("", "2", "1", "0", "2*", "", ""), # Specifies each label for the x axis
                     expand = expansion(mult = c(0.0, 0.01))) +
  
  # Specify colour/fill of points based on clinical classification (Pathogenic/Likely Pathogenic or Unreported)
  scale_fill_manual(values = classification_fill, 
                    drop = FALSE) +
  
  scale_colour_manual(values = classification_colors, 
                      drop = FALSE) +
  
  # Specify point shape based on variant source (ClinVar or UKB/AoU)
  scale_shape_manual(values = classification_shapes) +
  
  # Modify axes and titles
  labs(x = "Genomic position (GRCh38/hg38)",
       y = expression(paste("    UKB/AoU (allele count)                   ClinVar (assertion level)")),
       title = expression(plain("GENE:")*italic(" RB1")), 
       subtitle = "MANE select transcript ID: ENST00000267163.6/NM_000321.3") + 
  
  # Updating plot theme
  theme(legend.background = element_rect(),
        
        legend.title = element_blank(),
        
        legend.text = element_text(colour = "black", size = 10),
        
        legend.key = element_rect(color = "grey"),
        
        legend.key.size = unit(0.2, "cm"),
        
        legend.key.width = unit(0.5,"cm"),
        
        legend.position =  "top",
        
        legend.justification = "left",
        
        legend.box.just = "left",
        
        axis.line.y = element_line(color = "black", linewidth = 0.5), 
        
        axis.title.y=element_text(colour="black",size=10),
        
        axis.text.y=element_text(color="black", size=8, hjust = 1, vjust = 1),
        
        axis.ticks.y=element_line(color="black", linewidth = 0.5),
        
        axis.title.x = element_text(colour="black",size=10),
        
        axis.text.x=element_text(color="black", size=8, angle = 90, hjust = 1, vjust = 1),
        
        axis.ticks.x=element_line(color="black", linewidth = 0.5),
        
        title = element_text(colour = "black", size = 12),
        
        subtitle = element_text(colour = "black", size = 10),
        
        panel.grid.minor.x = element_blank(),
        
        panel.grid.major.x = element_blank(),
        
        panel.grid.minor.y = element_blank(),
        
        panel.grid.major.y = element_blank(),
        
        plot.margin = margin(t = 20, r = 20, b = 30, l = 60))


# Adds the NMD escape regions onto the plot
# Specify the 5' and 3' NMD escape regions via xmin and xmax
p2<-p1 +
  annotate("rect", xmin=48303913, xmax = 48304063, ymin = -1.5, ymax = 1.5, alpha = .2,fill = "dodgerblue1") + #5' NMD escape region
  annotate("rect", xmin =48477349, xmax = 48480071, ymin = -1.5, ymax = 1.5, alpha = .2,fill = "dodgerblue1") + #3' NMD escape region
  theme(plot.margin = unit(c(0, 0, 0, 0.4), "cm")) 


# Adds the exons onto main plot
p3 <- p2 +
  geom_rect(data = df_rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill='hotpink',
            alpha=1,
            colour='hotpink',
            size=1.5,
            inherit.aes = FALSE)


# Add in breaks between exons
n_exon <- nrow(exon_coordinates) #Number of exons in the gene


break_lines <- data.frame(
  x = c(exon_coordinates$end[1:(n_exon - 3)] + 50, # Final 3 exons are shown without breaks ((n_exon - 3))
        exon_coordinates$start[2:(n_exon - 2)] - 50))

p3b <- p3 +
  geom_segment(
    data = break_lines,
    aes(x = x, xend = x, y = -0.04, yend = 0.04),  # choose height to match your y-lims
    inherit.aes = FALSE,
    colour = "grey40",
    linewidth = 0.5,
    linetype = "solid")


# Add in y-axis breaks
p4 <- p3b

for (i in 1:(nrow(exon_coordinates) -3)) {
  p4 <- p4 + scale_x_break(c(exon_coordinates$end[i] + 50, exon_coordinates$start[i+1] - 50), 
                           scales=0.90, 
                           ticklabels = NULL)} 

left_boundary_ticks <- exon_coordinates$end[1:(n_exon-2)]

# Re-scale and label the yaxis for the final plot
p5 <- p4 +
  scale_x_continuous(limits = c(48303751, 48481890),
                     breaks = left_boundary_ticks,
                     labels = function(x) format(x, big.mark = ",", scientific = FALSE),
                     expand = expansion(mult = c(0.05, 0.05))) +
  
  theme(axis.line.x = element_line(),
        
        axis.text.x.bottom  = element_text(color="black", size=8, angle = 90, 
                                           hjust = 1, vjust = 1),
        axis.ticks.x.bottom = element_line(color="black", linewidth = 0.5),
        
        
        axis.text.x.top   = element_blank(),
        axis.ticks.x.top  = element_blank(),
        axis.title.x.top  = element_blank(),
        axis.line.x.top   = element_blank())


print(p5)

# Save final plot
ggsave("RB1_lollipop.png", p5, width = 15, height = 6, units = "in", dpi = 600)








