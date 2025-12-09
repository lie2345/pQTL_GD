# ==============================================================================  
# Mendelian Randomization Mediation Analysis  
# ==============================================================================  
# This script performs mediation analysis using Mendelian Randomization (MR)   
# summary statistics. It calculates direct effects, indirect (mediated) effects,  
# and bootstrapped confidence intervals for the mediation effect.  
#  
# Author: [Huanxian Luo]  
# Date: [2025-11-11]  
# Paper: [DNA Methylation-Mediated Ferroptosis Pathway Promotes Graves' Disease
# Risk: A Multi-Omics Mendelian Randomization Study]  
# ==============================================================================
# The analysis requires three sets of MR results:  
# 1. Exposure → Outcome (total effect)  
# 2. Exposure → Mediator (path a)  
# 3. Mediator → Outcome (path b)  
# ==============================================================================  

# -----------------------------------------------------------------------------  
# Load required packages  
# -----------------------------------------------------------------------------  
library(data.table)  # For efficient data reading and manipulation  
library(dplyr)       # For data filtering and manipulation  
library(boot)        # For bootstrapping confidence intervals  

# Set random seed for reproducibility  
set.seed(123)  

# -----------------------------------------------------------------------------  
# Step 1: Define the exposure and load data  
# -----------------------------------------------------------------------------  
# Define exposure name for filtering and output file naming  
exposure_name <- "ZDHHC4"  # Name of the exposure variable  

# Load MR results datasets  
# Total effect: Exposure → Outcome results  
data1 <- fread("ExposureOutcomeResult.tsv", header = TRUE, sep = "\t")  

# Path a: Exposure → Mediator results  
data2 <- fread("ExposureMediatorResult.tsv")  

# Path b: Mediator → Outcome results  
data3 <- fread("MediatorOutcomeResult.tsv")  

# -----------------------------------------------------------------------------  
# Step 2: Data preprocessing and validation  
# -----------------------------------------------------------------------------  
# Check if exposure exists in the data  
if(!exposure_name %in% data1$id) {  
  stop(paste("Exposure", exposure_name, "not found in exposure-outcome results"))  
}  

# Filter mediator-outcome results to exclude potentially invalid instruments  
# Remove results with evidence of heterogeneity (p < 0.05) or pleiotropy (Egger p < 0.05)  
data3 <- data3 %>% filter(Heter.p > 0.05 & egger_intercept_pval > 0.05)  

# Format mediator-outcome data for merging  
data3 <- data.frame(  
  id = data3$id,                # Mediator ID  
  nsnp = data3$nsnp,            # Number of SNPs used in the analysis  
  method = "IVW",               # Method used (Inverse Variance Weighted)  
  b = data3$beta,               # Effect size  
  se = data3$se,                # Standard error  
  or = data3$or,                # Odds ratio  
  or_lci95 = data3$or_lci95,    # Lower 95% CI of OR  
  or_uci95 = data3$or_uci95,    # Upper 95% CI of OR  
  pvalue = data3$IVW_pvalue     # P-value  
)  

# -----------------------------------------------------------------------------  
# Step 3: Merge datasets and calculate mediation effects  
# -----------------------------------------------------------------------------  
# Merge exposure-mediator and mediator-outcome data by mediator ID  
md2 <- merge(data2, data3, by = "id")  

# Extract total effect of exposure on outcome  
md2$beta_all <- data1$beta[data1$id == exposure_name]  

# Path a effect: Exposure → Mediator  
md2$betaa <- md2$b.x  # Effect from exposure to mediator (a path)  

# Path b effect: Mediator → Outcome  
md2$betab <- md2$b.y  # Effect from mediator to outcome (b path)  

# Calculate indirect effect (mediated effect)  
# Indirect effect = a*b (product of path coefficients)  
md2$betac <- md2$betaa * md2$betab  

# Calculate direct effect (total effect minus indirect effect)  
md2$beta_dir <- md2$beta_all - md2$betac  

# Calculate standard error for the mediated effect using the delta method  
# Formula: SE(ab) = √(a²*SE(b)² + b²*SE(a)²)  
md2$se <- sqrt(md2$betaa^2 * md2$se.y^2 + md2$betab^2 * md2$se.x^2)  

# -----------------------------------------------------------------------------  
# Step 4: Define bootstrapping functions for confidence intervals  
# -----------------------------------------------------------------------------  
# Function to calculate mediation effect for bootstrapping  
med_effect <- function(data, indices) {  
  # Sample with replacement  
  d <- data[indices, ]  
  
  # Generate random effect estimates from normal distributions  
  # centered at the point estimates with SEs as the standard deviations  
  beta_xz <- rnorm(1, d$b.x, d$se.x)   # Exposure → Mediator effect  
  beta_zy <- rnorm(1, d$b.y, d$se.y)   # Mediator → Outcome effect  
  
  # Check for invalid values (NA or Inf)  
  if(any(is.na(c(beta_xz, beta_zy))) || any(is.infinite(c(beta_xz, beta_zy)))) {  
    return(NA)  
  }  
  
  # Calculate mediation effect (indirect effect)  
  med <- beta_xz * beta_zy  
  return(med)  
}   

# Function to extract confidence intervals from bootstrap results  
get_ci <- function(boot_results, conf = 0.95) {  
  ci_results <- list()  
  
  # Use percentile method for confidence intervals  
  # (Other methods like BCa could be added if needed)  
  types <- c("perc")  
  
  for(type in types) {  
    tryCatch({  
      ci_results[[type]] <- boot.ci(boot_results,   
                                   type = type,   
                                   conf = conf)  
      cat(sprintf("Successfully calculated %s confidence interval\n", type))  
    }, error = function(e) {  
      cat(sprintf("Failed to calculate %s confidence interval: %s\n",   
                 type, e$message))  
      ci_results[[type]] <- NULL  
    })  
  }  
  
  return(ci_results)  
}  

# -----------------------------------------------------------------------------  
# Step 5: Perform bootstrapping for each mediator  
# -----------------------------------------------------------------------------  
# Initialize columns for bootstrap results  
md2$prec_ci_lower <- NA  # Lower bound of CI  
md2$prec_ci_upper <- NA  # Upper bound of CI  
md2$prec_lci_p <- NA     # Lower bound of proportion mediated  
md2$prec_uci_p <- NA     # Upper bound of proportion mediated  
md2$p_value <- NA        # P-value from bootstrap  

# Display information before starting bootstrap  
cat(sprintf("\nPerforming bootstrap analysis for %d potential mediators\n", nrow(md2)))  
cat("This may take some time...\n")  

# Perform bootstrapping for each mediator  
for(i in 1:nrow(md2)) {  
  # Display progress  
  cat(sprintf("Analyzing mediator %d of %d: %s\n", i, nrow(md2), md2$id[i]))  
  
  # Extract current row data  
  current_data <- md2[i, ]  
  
  # Execute bootstrap with 10,000 replications  
  boot_result <- boot(  
    data = current_data,  
    statistic = med_effect,  
    R = 10000  
  )  
  
  # Calculate confidence intervals  
  ci_results <- get_ci(boot_result)   
  
  # Calculate empirical p-value from bootstrap distribution  
  # (two-tailed test based on proportion of estimates crossing zero)  
  p_value <- 2 * min(  
    mean(boot_result$t >= 0, na.rm = TRUE),  
    mean(boot_result$t < 0, na.rm = TRUE)  
  )  
  
  # Store results in the data frame  
  md2$prec_ci_lower[i] <- ci_results[["perc"]][["percent"]][4]  # 2.5% percentile  
  md2$prec_ci_upper[i] <- ci_results[["perc"]][["percent"]][5]  # 97.5% percentile  
  
  # Calculate proportion of total effect that is mediated (with CI)  
  md2$prec_lci_p[i] <- md2$prec_ci_lower[i] / md2$beta_all  
  md2$prec_uci_p[i] <- md2$prec_ci_upper[i] / md2$beta_all  
  md2$p_value[i] <- p_value  
}  

# -----------------------------------------------------------------------------  
# Step 6: Save results and provide summary  
# -----------------------------------------------------------------------------  
# Define output filename  
output_file <- paste0(exposure_name, "MediationResult.tsv")  

# Save results to tab-separated file  
write.table(md2, output_file, sep = "\t", quote = FALSE, row.names = FALSE)  

# Display summary of results  
cat("\n=== Mediation Analysis Results ===\n")  
cat(sprintf("Exposure: %s\n", exposure_name))  
cat(sprintf("Total effect: %.4f\n", md2$beta_all[1]))  
cat(sprintf("Number of potential mediators analyzed: %d\n", nrow(md2)))  
cat(sprintf("Mediators with significant mediation (p < 0.05): %d\n",   
            sum(md2$p_value < 0.05, na.rm = TRUE)))  
cat(sprintf("Results saved to: %s\n", output_file))  
cat("===============================\n")
