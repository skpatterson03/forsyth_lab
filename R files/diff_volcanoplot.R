# Load necessary libraries
library(ggplot2)
library(dplyr)

# Load in the data
data <- read.csv("../data/DEG_methyl&trans_Loh2021_26695.csv")

pval_cutoff <- 0.01
log2FC_cutoff <- 1

# Prepare the data
# Calculate log2 fold change and ensure it's numeric
data$log2_fold_change <- log2(data$fold)
data$log10_fdr <- -log10(data$FDR)

# Create the volcano plot
volcano_plot <- ggplot(data, aes(x = log2_fold_change, y = log10_fdr)) +
  # Base scatter plot with color coding
  geom_point(aes(color = case_when(
    log2_fold_change > 2 ~ "Upregulated",
    log2_fold_change < -2 ~ "Downregulated",
    log2_fold_change >= -1 & log2_fold_change <= 1 ~ "Not Significant",
    ctrl.differences == 0 & (log2_fold_change > 1 | log2_fold_change < -1) ~ "Marginally Significant (No Diff)",
    ctrl.differences > 0 & (log2_fold_change > 1 | log2_fold_change < -1) ~ "Marginally Significant (With Diff)"
  )), alpha = 0.7) +
  
  # Color scheme
  scale_color_manual(values = c(
    "Upregulated" = "red", 
    "Downregulated" = "blue", 
    "Not Significant" = "grey", 
    "Marginally Significant (No Diff)" = "lightpink", 
    "Marginally Significant (With Diff)" = "orange"
  )) +
  
  # Labels for genes with significant fold change
  geom_text_repel(
    data = subset(data, abs(log2_fold_change) >= 2),
    aes(label = paste0(gene, "\nΔ", ctrl.differences)),
    size = 3,
    box.padding = 0.5,
    max.overlaps = Inf
  ) +
  
  # Thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "gray") +
  
  # Titles and labels
  labs(
    title = "H. pylori Methylation Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10(FDR)",
    color = "Regulation Status"
  ) +
  
  # Theme
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )
print(volcano_plot)
