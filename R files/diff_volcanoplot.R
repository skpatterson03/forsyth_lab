setwd("/Users/sarahpatterson/Desktop/Forsyth_Lab/samples_and_codes/pre_processed/volcano_plts")
# Load in the data
data <- read.csv("rna_seq_merged_methyl.csv")

pval_cutoff <- 0.01
log2FC_cutoff <- 1


library(ggplot2)
library(dplyr)
library(ggrepel)

# Prepare data
plot_data <- data %>%
  filter(Total_mod_count > 0) %>%
  mutate(
    mod_location = case_when(
      CodingRegion_mod_count > 0 & Promoter_mod_count > 0 ~ "Both (Coding & Promoter)",
      CodingRegion_mod_count > 0 & Promoter_mod_count == 0 ~ "Coding Region Only",
      CodingRegion_mod_count == 0 & Promoter_mod_count > 0 ~ "Promoter Only"
    ),
    
    neg_log10_FDR = -log10(FDR),
    
    # Only label genes with |log2FC| > 2
    point_label = ifelse(
      log2FC > 2 | log2FC < -2,
      paste0(Gene_ID, " (Δ", Total_mod_count, ")"),
      NA
    )
  )

mod_colors <- c(
  "Coding Region Only"       = "red",  
  "Promoter Only"            = "blue",  
  "Both (Coding & Promoter)" = "purple"   
)
  

ggplot(plot_data, aes(x = log2FC, y = neg_log10_FDR)) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  
  geom_point(
    aes(color = mod_location),
    shape = 16,
    size  = 2.5,
    alpha = 0.7
  ) +
  
  geom_text_repel(
    aes(label = point_label, color = mod_location),
    size          = 3,
    fontface      = "bold",
    box.padding   = 0.4,
    point.padding = 0.3,
    segment.size  = 0.3,
    segment.alpha = 0.6,
    max.overlaps  = Inf,
    show.legend   = FALSE
  ) +
  
  scale_color_manual(values = mod_colors, name = "Modification Location") +
  
  labs(
    title    = "Volcano Plot — ctrl Differential Modifications",
    subtitle = "Labels show Gene ID and number of methylation differences (Δ)",
    x        = expression(log[2]~Fold~Change),
    y        = expression(-log[10]~FDR)
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    legend.position  = "right",
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
