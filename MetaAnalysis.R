# ==============================================================================
# Meta-Analysis of MR Results Across Discovery and Replication Cohorts
# ==============================================================================
# This script performs fixed-effects inverse-variance weighted meta-analysis
# to pool effect estimates from deCODE (discovery) and UKB-PPP (replication) cohorts.
#
# Author: Huanxian Luo
# Date: 2025-11-11
# Paper: DNA Methylation-Mediated Ferroptosis Pathway Promotes Graves' Disease
#        Risk: A Multi-Omics Mendelian Randomization Study
# ==============================================================================

# -----------------------------------------------------------------------------
# Load required packages
# -----------------------------------------------------------------------------
# Install packages if not already installed
if (!require("meta")) install.packages("meta")
if (!require("data.table")) install.packages("data.table")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")

library(meta)
library(data.table)
library(dplyr)
library(ggplot2)

# Set random seed for reproducibility
set.seed(123)

# -----------------------------------------------------------------------------
# Step 1: Load and prepare data
# -----------------------------------------------------------------------------
# Set working directory (modify as needed)
# setwd("path/to/your/directory")

# Read input data containing OR and 95% CI from both cohorts
# Expected format:
# id           or          or_lci95    or_uci95
# deCODE       1.14        1.09        1.19
# UKB-PPP      1.10        1.07        1.13

meta_input <- read.table("meta_input.txt", header = TRUE, sep = "\t", 
                         stringsAsFactors = FALSE)

# Display input data
cat("=== Input Data ===\n")
print(meta_input)

# -----------------------------------------------------------------------------
# Step 2: Data transformation for meta-analysis
# -----------------------------------------------------------------------------
# Meta-analysis requires effect estimates and standard errors on the log scale

# Convert OR to log(OR)
lnor <- log(meta_input$or)

# Convert confidence interval bounds to log scale
lnuci <- log(meta_input$or_uci95)
lnlci <- log(meta_input$or_lci95)

# Calculate standard error from confidence interval
# Formula: SE = (ln(UCI) - ln(LCI)) / (2 * 1.96)
# This is based on: 95% CI = estimate ± 1.96 × SE
selnor <- (lnuci - lnlci) / (2 * 1.96)

# -----------------------------------------------------------------------------
# Step 3: Perform meta-analysis
# -----------------------------------------------------------------------------
# Use metagen() function with inverse-variance weighting
# sm = "OR" ensures results are back-transformed to OR scale
meta_result <- metagen(TE = lnor,           # Treatment effect (log OR)
                       seTE = selnor,        # Standard error of log OR
                       sm = "OR",            # Summary measure (Odds Ratio)
                       data = meta_input,
                       studlab = meta_input$id,  # Study labels
                       fixed = TRUE,         # Fixed-effect model
                       random = FALSE,       # Do not compute random-effects model
                       method.tau = "DL",    # DerSimonian-Laird for tau² (if random = TRUE)
                       hakn = FALSE,         # Hartung-Knapp adjustment
                       prediction = FALSE)   # Prediction interval

# Display detailed results
cat("\n=== Meta-Analysis Results ===\n")
summary(meta_result)

# -----------------------------------------------------------------------------
# Step 4: Extract and display key results
# -----------------------------------------------------------------------------
# Extract pooled OR and 95% CI
pooled_or <- exp(meta_result$TE.fixed)
pooled_lci <- exp(meta_result$lower.fixed)
pooled_uci <- exp(meta_result$upper.fixed)
pooled_p <- meta_result$pval.fixed

# Extract heterogeneity statistics
i2 <- meta_result$I2
i2_lci <- meta_result$lower.I2
i2_uci <- meta_result$upper.I2
tau2 <- meta_result$tau2
q_stat <- meta_result$Q
q_p <- meta_result$pval.Q

# Display summary
cat("\n=== Summary Statistics ===\n")
cat(sprintf("Pooled OR: %.2f (95%% CI: %.2f-%.2f)\n", pooled_or, pooled_lci, pooled_uci))
cat(sprintf("P-value: %.2e\n", pooled_p))
cat(sprintf("I² statistic: %.1f%% (95%% CI: %.1f%%-%.1f%%)\n", i2, i2_lci, i2_uci))
cat(sprintf("Cochran's Q: %.2f (P = %.3f)\n", q_stat, q_p))
cat(sprintf("Tau²: %.4f\n", tau2))

# Interpret heterogeneity
if(i2 < 25) {
  cat("Heterogeneity: Low\n")
} else if(i2 < 50) {
  cat("Heterogeneity: Moderate\n")
} else if(i2 < 75) {
  cat("Heterogeneity: Substantial\n")
} else {
  cat("Heterogeneity: Considerable\n")
}

# -----------------------------------------------------------------------------
# Step 5: Create forest plot
# -----------------------------------------------------------------------------
# Generate publication-quality forest plot
pdf("meta_analysis_forest_plot.pdf", width = 10, height = 6)
forest(meta_result,
       sortvar = TE,                    # Sort by effect size
       prediction = FALSE,
       print.tau2 = FALSE,
       print.I2 = TRUE,                 # Show I² statistic
       print.I2.ci = TRUE,              # Show I² confidence interval
       col.square = "navy",             # Color for individual studies
       col.square.lines = "navy",
       col.diamond = "maroon",          # Color for pooled estimate
       col.diamond.lines = "maroon",
       col.predict = "red",
       print.Q = TRUE,                  # Show Cochran's Q
       digits = 2,                      # Decimal places
       fontsize = 12,
       colgap.forest.left = "1cm",
       colgap.forest.right = "1cm",
       xlim = c(0.9, 1.3),              # X-axis limits (adjust as needed)
       xlab = "Odds Ratio (95% CI)")
dev.off()

# Also save as PNG
png("meta_analysis_forest_plot.png", width = 3000, height = 1800, res = 300)
forest(meta_result,
       sortvar = TE,
       prediction = FALSE,
       print.tau2 = FALSE,
       print.I2 = TRUE,
       print.I2.ci = TRUE,
       col.square = "navy",
       col.square.lines = "navy",
       col.diamond = "maroon",
       col.diamond.lines = "maroon",
       col.predict = "red",
       print.Q = TRUE,
       digits = 2,
       fontsize = 12,
       colgap.forest.left = "1cm",
       colgap.forest.right = "1cm",
       xlim = c(0.9, 1.3),
       xlab = "Odds Ratio (95% CI)")
dev.off()

cat("\nForest plots saved as 'meta_analysis_forest_plot.pdf' and '.png'\n")

# -----------------------------------------------------------------------------
# Step 6: Save results to file
# -----------------------------------------------------------------------------
# Create results data frame
results_df <- data.frame(
  Cohort = c(meta_input$id, "Pooled"),
  OR = c(meta_input$or, pooled_or),
  Lower_CI = c(meta_input$or_lci95, pooled_lci),
  Upper_CI = c(meta_input$or_uci95, pooled_uci),
  P_value = c(rep(NA, nrow(meta_input)), pooled_p),
  I2 = c(rep(NA, nrow(meta_input)), i2),
  Q_statistic = c(rep(NA, nrow(meta_input)), q_stat),
  Q_p_value = c(rep(NA, nrow(meta_input)), q_p)
)

# Save to file
write.table(results_df, "meta_analysis_results.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nResults saved to 'meta_analysis_results.txt'\n")

# -----------------------------------------------------------------------------
# Step 7: Sensitivity analysis - Leave-one-out
# -----------------------------------------------------------------------------
if(nrow(meta_input) > 2) {
  cat("\n=== Leave-One-Out Sensitivity Analysis ===\n")
  
  loo_results <- data.frame(
    Excluded = character(),
    Pooled_OR = numeric(),
    Lower_CI = numeric(),
    Upper_CI = numeric(),
    I2 = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(i in 1:nrow(meta_input)) {
    # Perform meta-analysis excluding study i
    temp_lnor <- lnor[-i]
    temp_selnor <- selnor[-i]
    temp_studlab <- meta_input$id[-i]
    
    temp_meta <- metagen(TE = temp_lnor,
                        seTE = temp_selnor,
                        sm = "OR",
                        studlab = temp_studlab,
                        fixed = TRUE,
                        random = FALSE)
    
    # Store results
    loo_results <- rbind(loo_results, data.frame(
      Excluded = meta_input$id[i],
      Pooled_OR = exp(temp_meta$TE.fixed),
      Lower_CI = exp(temp_meta$lower.fixed),
      Upper_CI = exp(temp_meta$upper.fixed),
      I2 = temp_meta$I2
    ))
  }
  
  print(loo_results)
  write.table(loo_results, "leave_one_out_analysis.txt", 
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\nLeave-one-out results saved to 'leave_one_out_analysis.txt'\n")
}

cat("\n=== Meta-analysis completed successfully ===\n")
