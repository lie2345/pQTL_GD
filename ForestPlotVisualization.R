# ==============================================================================
# Forest Plot for Multi-Cohort MR Results
# ==============================================================================
# This script creates a publication-quality forest plot comparing Mendelian
# Randomization results from discovery (deCODE) and replication (UKB-PPP) cohorts.
# The plot displays effect estimates (OR and 95% CI) for ferroptosis-related
# proteins associated with Graves' disease risk.
#
# Author: Huanxian Luo
# Date: 2025-01-15
# Paper: DNA Methylation-Mediated Ferroptosis Pathway Promotes Graves' Disease
#        Risk: A Multi-Omics Mendelian Randomization Study
# ==============================================================================

# -----------------------------------------------------------------------------
# Load required packages
# -----------------------------------------------------------------------------
# Install packages if not already installed
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("data.table")) install.packages("data.table")

library(ggplot2)
library(dplyr)
library(data.table)

# Set random seed for reproducibility
set.seed(123)

# -----------------------------------------------------------------------------
# Step 1: Set working directory and file paths
# -----------------------------------------------------------------------------
# NOTE: Modify these paths according to your directory structure
# Or use relative paths if running from project root

# Set working directory (modify as needed)
# setwd("path/to/your/project")

# Define input file paths
decode_result_file <- "results/deCODE_ferroptosis_results.tsv"
ukbppp_result_file <- "results/UKB_PPP_ferroptosis_results.tsv"
decode_full_file <- "results/deCODE_result_all.tsv"
ukbppp_full_file <- "results/UKB_PPP_result_all.tsv"

# Define output directory
output_dir <- "figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# Step 2: Load and filter significant results
# -----------------------------------------------------------------------------
cat("=== Loading MR Results ===\n")

# Load deCODE discovery cohort results
# Filter for significant associations with no evidence of heterogeneity or pleiotropy
mydata1 <- fread(decode_result_file, header = TRUE, sep = "\t") %>% 
  as.data.frame() %>%
  filter(IVW_pvalue < 0.05 &        # Significant association
         Heter.p > 0.05 &            # No significant heterogeneity
         egger_intercept_pval > 0.05) # No horizontal pleiotropy

cat(sprintf("deCODE: %d significant proteins (after filtering)\n", nrow(mydata1)))

# Load UKB-PPP replication cohort results with same filters
mydata2 <- fread(ukbppp_result_file, header = TRUE, sep = "\t") %>% 
  as.data.frame() %>%
  filter(IVW_pvalue < 0.05 & 
         Heter.p > 0.05 & 
         egger_intercept_pval > 0.05)

cat(sprintf("UKB-PPP: %d significant proteins (after filtering)\n", nrow(mydata2)))

# -----------------------------------------------------------------------------
# Step 3: Identify genes significant in either cohort
# -----------------------------------------------------------------------------
# Find union of genes significant in either cohort
gene_union <- union(mydata1$gene, mydata2$gene)

cat(sprintf("\nTotal unique genes significant in at least one cohort: %d\n", 
            length(gene_union)))

# -----------------------------------------------------------------------------
# Step 4: Load complete results for overlapping genes
# -----------------------------------------------------------------------------
# Load complete deCODE results
mydata1 <- fread(decode_full_file, header = TRUE, sep = "\t") %>% 
  as.data.frame()

# Extract gene names from protein IDs
# Expected format: "ID_number_GENE_protein"
mydata1$gene <- sapply(strsplit(mydata1$id, "_"), function(x) x[3])

# Filter for genes of interest and remove problematic variants
mydata1 <- mydata1 %>% 
  filter(gene %in% gene_union) %>%
  na.omit() %>%
  # Remove specific variants with known issues (modify as needed)
  filter(id != "6168_11_TP53_p53_R175H") %>%
  filter(id != "19567_1_EGFR_EGFRvIII") %>%
  mutate(group = "deCODE")

cat(sprintf("deCODE: %d protein-gene combinations for forest plot\n", nrow(mydata1)))

# Load complete UKB-PPP results
mydata2 <- fread(ukbppp_full_file, header = TRUE, sep = "\t") %>% 
  as.data.frame()

# Extract gene names (different ID format)
mydata2$gene <- sapply(strsplit(mydata2$id, "_"), function(x) x[1])

mydata2 <- mydata2 %>% 
  filter(gene %in% gene_union) %>%
  na.omit() %>%
  mutate(group = "UKB-PPP")

cat(sprintf("UKB-PPP: %d protein-gene combinations for forest plot\n", nrow(mydata2)))

# -----------------------------------------------------------------------------
# Step 5: Combine and sort data for visualization
# -----------------------------------------------------------------------------
# Combine both cohorts
mydata <- rbind(mydata1, mydata2) %>%
  # Calculate minimum p-value for each gene (for sorting)
  group_by(gene) %>%
  mutate(min_pvalue = min(IVW_pvalue)) %>%
  ungroup() %>%
  # Sort by: 1) gene-level min p-value, 2) gene name, 3) within-gene p-value
  arrange(min_pvalue, gene, IVW_pvalue) %>%
  select(-min_pvalue) %>%
  # Create gene grouping variables for alternating background colors
  mutate(
    gene_numeric = as.numeric(factor(gene, levels = unique(gene))),
    gene_even = gene_numeric %% 2 == 0
  ) %>%
  # Rename columns for clarity
  rename(pvalue = IVW_pvalue) %>%
  mutate(
    label = gene,
    outcome = "Graves' disease",
    "P value" = pvalue
  )

# Extract simple labels from protein IDs
mydata$label <- sapply(strsplit(mydata$id, "_"), "[", 1)

cat(sprintf("\nTotal rows in combined dataset: %d\n", nrow(mydata)))

# -----------------------------------------------------------------------------
# Step 6: Format values for display
# -----------------------------------------------------------------------------
# Format P-values
mydata$pvalue_display <- ifelse(
  mydata$pvalue < 0.001,
  sprintf("%.2e", mydata$pvalue),  # Scientific notation for very small p-values
  sprintf("%.3f", mydata$pvalue)   # Standard notation otherwise
)

# Format OR and 95% CI as "OR(LCI, UCI)"
mydata$'OR(95%CI)' <- ifelse(
  is.na(mydata$or), "",
  sprintf('%.3f(%.3f, %.3f)',
          mydata$or, mydata$or_lci95, mydata$or_uci95)
)

# Replace NA values with empty strings
mydata[is.na(mydata)] <- " "

# -----------------------------------------------------------------------------
# Step 7: Prepare data for plotting
# -----------------------------------------------------------------------------
mydata <- mydata %>%
  # Show gene name only for first occurrence
  mutate(gene_display = ifelse(duplicated(gene), "", gene)) %>%
  # Reverse order (ggplot builds from bottom up)
  arrange(desc(row_number())) %>%
  # Add row identifiers
  mutate(
    row_id = row_number(),
    # Create alternating background grouping
    row_even = group == "deCODE",
    # Create gene-level grouping for alternating colors
    gene_group = as.numeric(factor(gene)),
    gene_even = gene_group %% 2 == 0
  )

# Store number of rows for plot dimensions
n_rows <- nrow(mydata)

cat(sprintf("Final dataset prepared with %d rows\n", n_rows))

# -----------------------------------------------------------------------------
# Step 8: Create forest plot
# -----------------------------------------------------------------------------
cat("\n=== Creating Forest Plot ===\n")

p <- ggplot(mydata, aes(y = factor(row_id, levels = rev(row_id)))) +
  
  # Add alternating background rectangles by gene
  geom_rect(aes(xmin = -Inf, xmax = 3.5,
                ymin = as.numeric(factor(row_id)) - 0.5,
                ymax = as.numeric(factor(row_id)) + 0.5,
                fill = gene_even),
            alpha = 0.3) +
  
  # Add vertical reference line at OR = 1
  geom_segment(x = 1, xend = 1,
               y = 0.5, yend = n_rows + 0.5,
               linetype = "dashed", color = "#76A2B9",
               linewidth = 0.8) +
  
  # Add tick marks on x-axis
  geom_segment(x = 0, xend = 0, y = 0.5, yend = 0.6,
               color = "grey50", linewidth = 0.5) +
  geom_segment(x = 2, xend = 2, y = 0.5, yend = 0.6,
               color = "grey50", linewidth = 0.5) +
  
  # Add confidence interval error bars
  geom_errorbarh(aes(xmin = or_lci95, xmax = or_uci95, color = or > 1),
                 height = 0.2, linewidth = 1) +
  
  # Add point estimates (OR values)
  geom_point(aes(x = or, color = or > 1), size = 3, shape = 16) +
  
  # Add text annotations - Gene names
  geom_text(aes(x = -1.1, label = gene_display),
            hjust = 0.3, size = 3.2) +
  
  # Add text annotations - Database/Cohort
  geom_text(aes(x = -0.55, label = group),
            hjust = 0.3, size = 3.2) +
  
  # Add text annotations - Number of SNPs
  geom_text(aes(x = 0, label = nsnp),
            hjust = 0, size = 3.2) +
  
  # Add text annotations - OR and 95% CI
  geom_text(aes(x = 1.8, label = `OR(95%CI)`),
            size = 3.2, hjust = 0) +
  
  # Add text annotations - P-values
  geom_text(aes(x = 3.0, label = pvalue_display),
            size = 3.2, hjust = 0) +
  
  # Add column headers
  annotate("text", x = -1.1, y = n_rows + 0.8,
           label = "Gene", size = 3.5, hjust = 0.3, fontface = "bold") +
  annotate("text", x = -0.55, y = n_rows + 0.8,
           label = "Database", size = 3.5, hjust = 0.3, fontface = "bold") +
  annotate("text", x = 0, y = n_rows + 0.8,
           label = "N-SNPs", size = 3.5, hjust = 0.3, fontface = "bold") +
  annotate("text", x = 1.8, y = n_rows + 0.8,
           label = "OR(95%CI)", size = 3.5, hjust = 0, fontface = "bold") +
  annotate("text", x = 3.0, y = n_rows + 0.8,
           label = "P value", size = 3.5, hjust = 0, fontface = "bold") +
  
  # Add top separator line
  geom_segment(x = -1.35, xend = 3.5,
               y = n_rows + 0.5,
               yend = n_rows + 0.5,
               color = "black", linewidth = 0.5) +
  
  # Add bottom separator line
  geom_segment(x = -1.35, xend = 3.5,
               y = 0.5,
               yend = 0.5,
               color = "black", linewidth = 0.5) +
  
  # Set color schemes
  # Point and error bar colors based on OR direction
  scale_color_manual(
    values = c("FALSE" = "#72B6A1",  # OR < 1 (protective)
               "TRUE" = "#4189C8"),   # OR > 1 (risk)
    guide = "none") +
  
  # Background alternating colors
  scale_fill_manual(
    values = c("TRUE" = "#C0C0BF",   # Even gene groups
               "FALSE" = "#FFFFFF"),  # Odd gene groups
    guide = "none") +
  
  # Configure x-axis
  scale_x_continuous(
    limits = c(-1.35, 3.5),
    breaks = c(0, 1, 2),
    labels = c("0", "1", "2"),
    expand = c(0, 0)) +
  
  # Configure y-axis
  scale_y_discrete(
    limits = c(as.character(1:n_rows), "title"),
    labels = c(rep("", n_rows), "")) +
  
  # Add plot title and labels
  labs(
    x = "Odds Ratio",
    y = "",
    title = "Forest Plot: deCODE and UK Biobank pQTL-Graves' Disease Associations",
    subtitle = paste0("Ferroptosis-related proteins (n = ", 
                     length(unique(mydata$gene)), " genes)")) +
  
  # Theme customization
  theme_minimal() +
  theme(
    # Title formatting
    plot.title = element_text(
      hjust = 0.5,
      vjust = 1,
      size = 12,
      face = "bold",
      margin = margin(b = 5)),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10,
      margin = margin(b = 10)),
    
    # Remove grid lines
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid = element_blank(),
    
    # Remove y-axis elements
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    
    # Background colors
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_blank(),
    
    # Legend
    legend.position = "none",
    
    # Axis labels
    axis.title.x = element_text(margin = margin(t = 10), size = 11),
    
    # Plot margins
    plot.margin = margin(t = 15, r = 15, b = 15, l = 15, unit = "pt")
  )

# Display plot
print(p)

# -----------------------------------------------------------------------------
# Step 9: Save forest plot to file
# -----------------------------------------------------------------------------
cat("\n=== Saving Forest Plot ===\n")

# Calculate dynamic height based on number of rows
# Rule of thumb: ~0.3 inches per row, minimum 6 inches
plot_height <- max(6, n_rows * 0.3 + 2)
plot_width <- 10

# Save as high-resolution TIFF
output_tiff <- file.path(output_dir, "forest_plot_multicohort.tif")
ggsave(output_tiff,
       plot = p,
       device = "tiff",
       width = plot_width,
       height = plot_height,
       units = "in",
       dpi = 600,
       compression = "lzw")

cat(sprintf("TIFF saved: %s\n", output_tiff))

# Save as PDF (vector format, good for editing)
output_pdf <- file.path(output_dir, "forest_plot_multicohort.pdf")
ggsave(output_pdf,
       plot = p,
       device = "pdf",
       width = plot_width,
       height = plot_height,
       units = "in")

cat(sprintf("PDF saved: %s\n", output_pdf))

# Save as PNG (for presentations/web)
output_png <- file.path(output_dir, "forest_plot_multicohort.png")
ggsave(output_png,
       plot = p,
       device = "png",
       width = plot_width,
       height = plot_height,
       units = "in",
       dpi = 300)

cat(sprintf("PNG saved: %s\n", output_png))

# -----------------------------------------------------------------------------
# Step 10: Generate summary statistics
# -----------------------------------------------------------------------------
cat("\n=== Summary Statistics ===\n")

# Count significant results by cohort
summary_stats <- mydata %>%
  group_by(group) %>%
  summarise(
    n_proteins = n(),
    n_unique_genes = n_distinct(gene),
    median_or = median(or, na.rm = TRUE),
    median_pvalue = median(pvalue, na.rm = TRUE),
    n_risk = sum(or > 1, na.rm = TRUE),
    n_protective = sum(or < 1, na.rm = TRUE)
  )

print(summary_stats)

# Save summary statistics
summary_file <- file.path(output_dir, "forest_plot_summary.txt")
write.table(summary_stats, summary_file, 
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("\nSummary statistics saved: %s\n", summary_file))

# -----------------------------------------------------------------------------
# Step 11: Save complete plotting data
# -----------------------------------------------------------------------------
# Save the processed data used for plotting (for reproducibility)
data_file <- file.path(output_dir, "forest_plot_data.txt")
write.table(mydata, data_file,
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Plotting data saved: %s\n", data_file))

cat("\n=== Forest plot creation completed successfully ===\n")

# -----------------------------------------------------------------------------
# Additional Notes
# -----------------------------------------------------------------------------
# Color scheme explanation:
# - Green (#72B6A1): Protective effects (OR < 1)
# - Blue (#4189C8): Risk effects (OR > 1)
# - Gray/White alternating: Different gene groups for easier reading
#
# The plot shows:
# - Gene names (only shown once per gene)
# - Database source (deCODE or UKB-PPP)
# - Number of genetic instruments (SNPs)
# - Effect size with 95% confidence interval
# - P-value from IVW MR analysis
#
# Filtering criteria applied:
# - IVW P < 0.05 (significant association)
# - Heterogeneity P > 0.05 (no evidence of heterogeneity)
# - Egger intercept P > 0.05 (no evidence of horizontal pleiotropy)
# ==============================================================================
