# pQTL_GD
MR analysis of pQTL and Graves‘ Disease

## Overview

This repository contains a collection of R scripts for performing different types of Mendelian Randomization (MR) analyses. These scripts enable causal inference between exposures and outcomes using genetic variants as instrumental variables. The toolkit includes standard MR analysis and mediation analysis.

## Scripts

### 1. MR.R
Standard Mendelian Randomization analysis script that implements multiple MR methods using summary statistics.

**Key features:**
- Multiple MR methods (Inverse Variance Weighted, Weighted Median, MR-Egger, etc.)
- Horizontal pleiotropy assessment
- Sensitivity analyses
- Visualization of MR results

### 2. MediationAnalysis.R
Performs mediation analysis using MR summary statistics to quantify direct and indirect (mediated) effects of genetic exposures on outcomes through potential mediators.

**Key features:**
- Calculates direct and indirect effects
- Bootstrap confidence intervals
- Proportion of effect mediated
- Filters for valid instruments

## Requirements

### R Packages
```r
# Common packages (available on CRAN)  
install.packages(c("data.table", "dplyr", "boot", "ggplot2", "meta"))  

# For MR.R 
install.packages("remotes")  
remotes::install_github("MRCIEU/TwoSampleMR")  

```

## Input Data

### MR.R
Requires exposure and outcome summary statistics with the following columns:
- SNP ID (rsid)
- Effect allele
- Other allele
- Effect allele frequency
- Beta coefficient
- Standard error
- P-value
- Sample size


### MediationAnalysis.R
Requires three sets of MR summary statistics:
1. **Exposure → Outcome results** (`ExposureOutcomeResult.tsv`)
2. **Exposure → Mediator results** (`ExposureMediatorResult.tsv`)
3. **Mediator → Outcome results** (`MediatorOutcomeResult.tsv`)

## Usage

### MR.R
```r
# Edit the input file paths and parameters at the beginning of the script
# Then run:
source("MR.R")
```

### MediationAnalysis.R
```r
# Edit the exposure name and file paths:
exposure_name <- "YourExposureName"

# Then run:
source("MediationAnalysis.R")
```

## Output

### MR.R
- `harmonise.tsv`: Harmonized exposure and outcome data
- `OR.tsv`: Odds ratios from different MR methods
- `pleiotropy.tsv`: MR-Egger intercept test for directional pleiotropy
- `heterogeneity.tsv`: Heterogeneity statistics including Cochran's Q and I²
- `singlesnpOR.txt`: Results from single-SNP analysis
- `forest.pdf`: Forest plot showing effect estimates for each SNP
- `sensitivity-analysis.pdf`: Leave-one-out analysis plot
- `funnelplot.pdf`: Funnel plot for detecting directional pleiotropy
- `presso.txt`: MR-PRESSO results for outlier detection and correction

### MediationAnalysis.R
- `[exposure_name]MediationResult.tsv`: Mediation analysis results including:
  - Direct and indirect effects
  - Bootstrap confidence intervals
  - Proportion of effect mediated
  - P-values

## Interpretation of Results

### MR Results
- **Causal estimate (beta)**: The effect of the exposure on the outcome
- **OR**: Odds ratio (for binary outcomes)
- **P-value**: Statistical significance of the causal effect
- **Heterogeneity tests**: Assess whether all SNPs provide consistent estimates (Q and I² statistics)
- **Pleiotropy tests**: MR-Egger intercept for assessing directional pleiotropy
- **Sensitivity analyses**: Leave-one-out and single-SNP analyses to identify influential outliers

### Mediation Results
- **Indirect effect**: Product of path a and path b coefficients
- **Direct effect**: Total effect minus indirect effect
- **Proportion mediated**: Percentage of total effect explained by the mediator
- **Bootstrap CI**: Confidence intervals for indirect effect
- **P-value**: Statistical significance of mediation

## References

1. Zhu, Z., Zheng, Z., Zhang, F., Wu, Y., Trzaskowski, M., Maier, R., Robinson, M. R., McGrath, J. J., Visscher, P. M., Wray, N. R., & Yang, J. (2018). Causal associations between risk factors and common diseases inferred from GWAS summary data. Nature Communications, 9, 224. https://doi.org/10.1038/s41467-017-02317-2.

2. Xue, A., Jiang, L., Makalic, E., Jia, T., Shadrin, A., Landén, M., Robinson, M. R., & Zhu, Z. (2024). Unravelling the complex causal effects of substance use behaviours on common diseases. Communications Medicine, 4(1), 43. https://doi.org/10.1038/s43856-024-00411-3.

## License

This project is licensed under the MIT License.

## Citation

If you use these scripts in your research, please cite:

```
Huanxian Luo. (2025). Mendelian Randomization Analysis Code. 
GitHub Repository. https://github.com/lie2345/pQTL_GD
```

## Contributing

Contributions to improve the code are welcome. Please feel free to submit issues or pull requests.

## Contact

Huanxian Luo - luohuanx@alumni.sysu.edu.cn

Project Link: [https://github.com/lie2345/pQTL_GD](https://github.com/lie2345/pQTL_GD)
