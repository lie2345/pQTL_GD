# ============================================================================  
# Mendelian Randomization Analysis Script  
# ============================================================================  
# This script performs Mendelian Randomization (MR) analysis to estimate causal  
# effects between an exposure and outcome using genetic variants as   
# instrumental variables.  
#   
# Author: [Huanxian Luo]  
# Date: [2025-11-11]  
# Paper: [DNA Methylation-Mediated Ferroptosis Pathway Promotes Graves' Disease
# Risk: A Multi-Omics Mendelian Randomization Study]  
# ============================================================================  

# -----------------------------------------------------------------------------  
# Function to check and install required packages if not already installed  
# -----------------------------------------------------------------------------  
install_if_missing <- function(package) {  
  # Check if package is already installed; if not, attempt to install it  
  if (!require(package, character.only = TRUE, quietly = TRUE)) {  
    # First try installing from CRAN  
    if (!requireNamespace("BiocManager", quietly = TRUE))  
      install.packages("BiocManager")  
    
    tryCatch({  
      install.packages(package)  
    }, error = function(e) {  
      # If CRAN installation fails, try from Bioconductor  
      BiocManager::install(package)  
    })  
    
    library(package, character.only = TRUE)  
    cat(paste0("Package installed and loaded: ", package, "\n"))  
  } else {  
    cat(paste0("Package loaded: ", package, "\n"))  
  }  
}  

# Load required packages for Mendelian Randomization analysis  
packages <- c("TwoSampleMR", "data.table", "tidyverse", "ggplot2")  
invisible(sapply(packages, install_if_missing))   

# -----------------------------------------------------------------------------  
# Step 1: Read exposure data (genetic instruments)  
# -----------------------------------------------------------------------------  
# Read the SNPs (genetic variants) to be used as instrumental variables  
# Format must include columns for SNP IDs, effect sizes, standard errors, etc.  
expo_rt <- read_exposure_data(  
  filename = "exposure_data.txt",  
  clump = FALSE,  
  sep = "\t",  
  snp_col = "SNP",          # Column containing SNP IDs  
  beta_col = "beta",        # Column containing effect sizes  
  se_col = "se",            # Column containing standard errors  
  effect_allele_col = "A1", # Column containing effect alleles  
  other_allele_col = "A2",  # Column containing other alleles  
  eaf_col = "eaf",          # Column containing effect allele frequency  
  pval_col = "p",           # Column containing p-values  
  samplesize_col = "n"      # Column containing sample size  
)  

# Filter SNPs to include only genome-wide significant variants (p < 5e-8)  
expo_rt <- subset(expo_rt, pval.exposure < 5e-8)  

# Perform LD clumping to ensure independent instruments  
# Parameters: 10,000kb window, r² threshold of 0.1, European population  
expo_rt <- clump_data(expo_rt, clump_kb = 10000, clump_r2 = 0.1, clump_p1 = 1, clump_p2 = 1, pop = "EUR")  

print(paste("Number of SNPs in exposure data:", nrow(expo_rt)))   

# -----------------------------------------------------------------------------  
# Step 2: Read outcome data  
# -----------------------------------------------------------------------------  
# Define column mapping for outcome data file  
col_map <- c(  
  SNP = "rsids",          # SNP identifier  
  beta = "beta",          # Effect size  
  se = "sebeta",          # Standard error  
  eaf = "af_alt",         # Effect allele frequency  
  effect_allele = "alt",  # Effect allele  
  other_allele = "ref",   # Other allele  
  pval = "pval",          # P-value  
  chr = "Chrom",          # Chromosome  
  pos = "Pos"             # Position  
)  
   
# Read outcome data and filter for SNPs present in exposure data  
# This ensures we only analyze SNPs that were identified as instruments  
outcome_dat <- fread(  
  file = "finngen_R12_F5_OCD.gz",  
  select = unname(col_map),           
  col.names = names(col_map)        
)[SNP %in% expo_rt$SNP]   
  
# Format outcome data for MR analysis using TwoSampleMR format  
outc_rt <- format_data(  
  as.data.frame(outcome_dat),   
  type = "outcome",   
  snps = NULL,   
  phenotype_col = "Phenotype",   
  snp_col = "SNP",   
  beta_col = "beta",   
  se_col = "se",   
  eaf_col = "eaf",   
  effect_allele_col = "effect_allele",   
  other_allele_col = "other_allele",   
  pval_col = "pval",   
  units_col = "units",   
  ncase_col = "ncase",   
  ncontrol_col = "ncontrol",   
  samplesize_col = "Samplesize",   
  gene_col = "gene",   
  id_col = "id",   
  min_pval = 1e-200,   
  log_pval = FALSE,   
  chr_col = "chr",   
  pos_col = "pos"  
)  
outc_rt$data_source.outcome <- "textfile"  
  
# Remove SNPs that are genome-wide significant in the outcome  
# This helps avoid potential reverse causation or horizontal pleiotropy  
outc_rt = outc_rt[outc_rt$pval.outcome > 5e-8, ]  
outc_rt$samplesize.outcome = 447227  # Set sample size for outcome dataset  

# -----------------------------------------------------------------------------  
# Step 3: Harmonize exposure and outcome data  
# -----------------------------------------------------------------------------  
# Align effect alleles between exposure and outcome datasets  
# action=2 means try to infer forward strand alleles, remove palindromic SNPs   
# with ambiguous allele frequencies (A/T or G/C with MAF~0.5)  
harm_rt <- harmonise_data(  
  exposure_dat = expo_rt,   
  outcome_dat = outc_rt,  
  action = 2  
)  
    
# -----------------------------------------------------------------------------  
# Step 4: Calculate instrument strength (F-statistic)  
# -----------------------------------------------------------------------------  
# F-statistic is used to assess instrument strength. Weak instruments (F<10)  
# can lead to bias. We calculate F-statistics with or without EAF depending on data.  
if (all(is.na(expo_rt$eaf.exposure))) {  
  # If effect allele frequency is missing, calculate F-statistic without EAF  
  harm_rt <- subset(harm_rt, mr_keep == TRUE)  
  # Calculate variance explained (R²) for each SNP  
  harm_rt$R2 <- (2 * (harm_rt$beta.exposure^2)) /  
    (2 * (harm_rt$beta.exposure^2) +  
       2 * harm_rt$samplesize.exposure * harm_rt$se.exposure^2)  
  # Calculate F-statistic for each SNP  
  harm_rt$f <- harm_rt$R2 * (harm_rt$samplesize.exposure - 2) / (1 - harm_rt$R2)  
  harm_rt$meanf <- mean(harm_rt$f)  
} else {    
  # Calculate with EAF - provides better estimate of R² and F-statistic  
  harm_rt$R2 <- (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure)) /  
    (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) +  
       2 * harm_rt$samplesize.exposure * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) * harm_rt$se.exposure^2)  
  harm_rt$f <- harm_rt$R2 * (harm_rt$samplesize.exposure - 2) / (1 - harm_rt$R2)  
  harm_rt$meanf <- mean(harm_rt$f)  
}  
    
# Filter SNPs to keep only strong instruments (F > 10) to avoid weak instrument bias  
harm_rt <- harm_rt[harm_rt$f > 10, ]  
    
# -----------------------------------------------------------------------------  
# Step 5: Steiger filtering - ensure instruments affect exposure more than outcome  
# -----------------------------------------------------------------------------  
# This removes SNPs that explain more variance in the outcome than the exposure  
# (which would suggest reverse causation or invalid instruments)  
harm_rt = steiger_filtering(harm_rt)  

# -----------------------------------------------------------------------------  
# Step 6: Perform MR analyses using various methods  
# -----------------------------------------------------------------------------  
# Run main MR analysis (includes Inverse Variance Weighted, MR Egger, Weighted Median)  
# Different methods have different assumptions about pleiotropy  
mr_result <- mr(harm_rt)  

# Convert results to odds ratios (for binary outcomes)  
# This transforms beta coefficients to odds ratios for easier interpretation  
result_or = generate_odds_ratios(mr_result)  

# -----------------------------------------------------------------------------  
# Step 7: Save harmonized data and main results  
# -----------------------------------------------------------------------------  
write.table(harm_rt, file = "harmonise.tsv", row.names = FALSE, sep = "\t", quote = FALSE)  
write.table(result_or[, 5:ncol(result_or)], file = "OR.tsv", row.names = FALSE, sep = "\t", quote = FALSE)  

# -----------------------------------------------------------------------------  
# Step 8: MR sensitivity analyses  
# -----------------------------------------------------------------------------  
# Test for horizontal pleiotropy using MR Egger intercept  
# A significant intercept suggests directional pleiotropy  
pleiotropy = mr_pleiotropy_test(harm_rt)  
write.table(pleiotropy, file = "pleiotropy.tsv", sep = "\t", quote = FALSE, row.names = FALSE)  

# Test for heterogeneity using Cochran's Q and calculate I² statistic  
# High heterogeneity suggests potential pleiotropy or instrument invalidity  
heterogeneity = mr_heterogeneity(harm_rt)  
heterogeneity$I2 <- ifelse(heterogeneity$Q > 0, (heterogeneity$Q - heterogeneity$Q_df) / heterogeneity$Q, 0)  

write.table(heterogeneity, file = "heterogeneity.tsv", sep = "\t", quote = FALSE, row.names = FALSE)  

# Perform single-SNP analysis (leave-one-in)  
# This shows the causal effect estimate when using each SNP individually  
singlesnp_res <- mr_singlesnp(harm_rt)  
singlesnpOR = generate_odds_ratios(singlesnp_res)  
write.table(singlesnpOR, file = "singlesnpOR.txt", row.names = FALSE, sep = "\t", quote = FALSE)  

# Create forest plot of single-SNP results  
# This visualizes the consistency of effects across individual instruments  
p2 <- mr_forest_plot(singlesnp_res)  
ggsave(p2[[1]], file = "forest.pdf", width = 8, height = 8)  

# Perform leave-one-out analysis to identify influential outliers  
# This helps detect if any single SNP drives the overall result  
sen_res <- mr_leaveoneout(harm_rt)  
p3 <- mr_leaveoneout_plot(sen_res)  
ggsave(p3[[1]], file = "sensitivity-analysis.pdf", width = 8, height = 8)  

# Create funnel plot to visualize potential directional pleiotropy  
# Asymmetry suggests directional pleiotropy  
p4 <- mr_funnel_plot(singlesnp_res)  
ggsave(p4[[1]], file = "funnelplot.pdf", width = 8, height = 8)  

# Run MR-PRESSO to detect and correct for horizontal pleiotropy  
# This method identifies and removes outlier SNPs that may exhibit pleiotropy  
presso = run_mr_presso(harm_rt, NbDistribution = 1000)  
capture.output(presso, file = "presso.txt")
