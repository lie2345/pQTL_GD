# pQTL_GD

Multi-Omics Mendelian Randomization Analysis of Ferroptosis-Related Proteins and Graves' Disease

## Overview

This repository contains a comprehensive collection of R scripts for performing Mendelian Randomization (MR) analyses investigating the causal relationships between ferroptosis-related proteins, DNA methylation, and Graves' disease risk. The toolkit implements standard two-sample MR, meta-analysis across multiple cohorts, mediation analysis, and publication-quality data visualization.


## Scripts

### 1. MR.R
Core Mendelian Randomization analysis script implementing multiple MR methods using summary statistics from pQTL and GWAS data.

**Key features:**
- Multiple MR methods (IVW, Weighted Median, MR-Egger, Weighted Mode)
- Horizontal pleiotropy assessment (MR-Egger intercept test)
- Heterogeneity testing (Cochran's Q, I² statistics)
- Outlier detection and correction (MR-PRESSO)
- Comprehensive sensitivity analyses (leave-one-out, single-SNP)
- Automated visualization (forest plots, funnel plots, scatter plots)

**Application:**
- Phase 1: Proteome-wide MR screening (ferroptosis proteins → Graves' disease)
- Discovery cohort: deCODE pQTL data
- Replication cohort: UK Biobank Pharma Proteomics Project

---

### 2. MetaAnalysis.R
Fixed-effects meta-analysis combining MR results from discovery and replication cohorts.

**Key features:**
- Inverse-variance weighted meta-analysis
- Heterogeneity assessment (I², Cochran's Q)
- Leave-one-out sensitivity analysis
- Publication-quality forest plots
- Automatic conversion between OR and log(OR) scales

**Application:**
- Pooling effect estimates from deCODE and UKB-PPP cohorts
- Evaluating consistency across independent populations
- Generating combined evidence for causal relationships

---

### 3. ForestPlotVisualization.R
Creates publication-ready forest plots comparing MR results across multiple cohorts and proteins.

**Key features:**
- Multi-cohort comparison visualization
- Gene-based grouping with alternating backgrounds
- Automatic result filtering (significance, heterogeneity, pleiotropy)
- Dynamic plot sizing based on number of results
- Multiple output formats (TIFF, PDF, PNG)
- Customizable color schemes for risk/protective effects

**Application:**
- Visual presentation of proteome-wide MR results
- Comparison of effect estimates between deCODE and UKB-PPP
- Identifying consistent signals across cohorts

---

### 4. MediationAnalysis.R
Two-step MR mediation analysis to quantify direct and indirect effects of exposures through mediators.

**Key features:**
- Product-of-coefficients mediation estimation
- Bootstrap confidence intervals (10,000 iterations)
- Mediation proportion calculation
- Automated filtering for valid mediators
- Empirical p-value calculation

**Application:**
- Phase 2: Testing DNA methylation → protein → Graves' disease pathways
- Quantifying proportion of effect mediated through epigenetic regulation
- Identifying regulatory CpG sites for ferroptosis proteins

---

## Requirements

### R Version
- R ≥ 4.4.3

### Required Packages

```r
# Install from CRAN
install.packages(c(
  "data.table",    # Fast data manipulation
  "dplyr",         # Data wrangling
  "ggplot2",       # Data visualization
  "boot",          # Bootstrap methods
  "meta"           # Meta-analysis
))

# Install TwoSampleMR (for MR.R)
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")

# Install MR-PRESSO (for MR.R)
install.packages("MRPRESSO")
```

---

## Data Sources

This project uses publicly available GWAS summary statistics:

| Data Type | Source | Description | Access |
|-----------|--------|-------------|--------|
| **pQTL (Discovery)** | deCODE Genetics | 35,559 Icelandic individuals, 4,907 proteins | https://www.decode.com/summarydata/ |
| **pQTL (Replication)** | UKB-PPP | 54,219 UK Biobank participants, 1,536 proteins | https://www.ukbiobank.ac.uk/ |
| **mQTL** | GoDMC | 27,000+ individuals, 420,509 CpG sites | https://www.godmc.org.uk/ |
| **Outcome (Graves' disease)** | FinnGen R12 | 3,012 cases, 497,336 controls | https://www.finngen.fi/en |
| **Ferroptosis genes** | FerrDb V2 | 324 experimentally validated genes | http://www.zhounan.org/ferrdb/ |

---

## Input Data Format

### MR.R
Exposure and outcome summary statistics with the following columns:

**Required columns:**
- `SNP`: rsID or variant identifier
- `effect_allele`: Effect/alternative allele
- `other_allele`: Reference allele
- `eaf`: Effect allele frequency
- `beta`: Effect size (log odds ratio or beta coefficient)
- `se`: Standard error of beta
- `pval`: P-value
- `samplesize`: Sample size (optional but recommended)

**Example:**
```
SNP           effect_allele  other_allele  eaf    beta      se        pval
rs1234567     A              G             0.35   0.045     0.012     0.0002
rs7654321     T              C             0.52  -0.023     0.015     0.125
```

---

### MetaAnalysis.R
Tab-separated file with meta-analysis input:

**Required columns:**
- `id`: Study/cohort identifier
- `or`: Odds ratio
- `or_lci95`: Lower 95% confidence interval
- `or_uci95`: Upper 95% confidence interval

**Example (`meta_input.txt`):**
```
id        or      or_lci95  or_uci95
deCODE    1.14    1.09      1.19
UKB-PPP   1.10    1.07      1.13
```

---

### ForestPlotVisualization.R
MR results files from both cohorts:

**Required columns:**
- `id`: Protein/gene identifier
- `gene`: Gene symbol
- `nsnp`: Number of SNPs used as instruments
- `IVW_pvalue`: P-value from IVW method
- `Heter.p`: Heterogeneity test p-value
- `egger_intercept_pval`: MR-Egger intercept p-value
- `or`: Odds ratio
- `or_lci95`: Lower 95% CI
- `or_uci95`: Upper 95% CI

---

### MediationAnalysis.R
Three separate MR result files:

1. **ExposureOutcomeResult.tsv**: Total effect (e.g., methylation → GD)
2. **ExposureMediatorResult.tsv**: Path a (e.g., methylation → protein)
3. **MediatorOutcomeResult.tsv**: Path b (e.g., protein → GD)

Each file should contain standard MR output columns (`id`, `beta`, `se`, `pval`, etc.)

---

## Usage

### 1. Standard MR Analysis

```r
# Edit input file paths and parameters
source("MR.R")
```

**Workflow:**
1. Load exposure and outcome data
2. Select genetic instruments (P < 5×10⁻⁸, LD r² < 0.1)
3. Harmonize datasets
4. Perform MR analyses
5. Generate diagnostic plots
6. Save results

---

### 2. Meta-Analysis

```r
# Prepare input file (meta_input.txt)
# Run meta-analysis
source("MetaAnalysis.R")
```

**Workflow:**
1. Reads pooled effect estimates from multiple cohorts
2. Converts OR to log scale
3. Calculates standard errors from 95% CI
4. Performs fixed-effects IVW meta-analysis
5. Assesses heterogeneity (I², Q statistic)
6. Generates forest plot

---

### 3. Forest Plot Generation

```r
# Update file paths to your MR results
# Modify filtering criteria if needed (lines 47-52)
source("ForestPlotVisualization.R")
```

**Customization options:**
- Adjust filtering thresholds (significance, heterogeneity, pleiotropy)
- Modify color schemes (`scale_color_manual`, `scale_fill_manual`)
- Change plot dimensions (`plot_width`, `plot_height`)
- Customize column positions (x-coordinates in `geom_text`)

---

### 4. Mediation Analysis

```r
# Set exposure name
exposure_name <- "DPEP1"

# Update file paths to your MR results
source("MediationAnalysis.R")
```

**Workflow:**
1. Loads three sets of MR results
2. Filters mediators based on validity criteria
3. Calculates indirect effect (a × b)
4. Performs bootstrap (10,000 replicates)
5. Calculates mediation proportion
6. Saves results with confidence intervals

---

## Output Files

### MR.R
| File | Description |
|------|-------------|
| `harmonise.txt` | Harmonized exposure-outcome data |
| `OR.txt` | Effect estimates from all MR methods |
| `pleiotropy.txt` | MR-Egger intercept test results |
| `heterogeneity.txt` | Cochran's Q and I² statistics |
| `singlesnpOR.txt` | Individual SNP effect estimates |
| `forest.pdf` | Forest plot of SNP-specific effects |
| `sensitivity-analysis.pdf` | Leave-one-out analysis plot |
| `funnelplot.pdf` | Funnel plot for pleiotropy detection |
| `presso.txt` | MR-PRESSO outlier test results |

---

### MetaAnalysis.R
| File | Description |
|------|-------------|
| `meta_analysis_forest_plot.pdf` | Publication-quality forest plot |
| `meta_analysis_results.txt` | Pooled estimates and heterogeneity stats |
| `leave_one_out_analysis.txt` | Sensitivity analysis results |

---

### ForestPlotVisualization.R
| File | Description |
|------|-------------|
| `forest_plot_multicohort.tif` | High-resolution TIFF (600 DPI) |
| `forest_plot_multicohort.pdf` | Vector format for editing |
| `forest_plot_multicohort.png` | Web-ready PNG (300 DPI) |
| `forest_plot_summary.txt` | Summary statistics by cohort |
| `forest_plot_data.txt` | Complete plotting dataset |

---

### MediationAnalysis.R
| File | Description |
|------|-------------|
| `[exposure_name]MediationResult.tsv` | Complete mediation analysis results |

**Output columns:**
- `id`: Mediator identifier
- `betaa`: Path a effect (exposure → mediator)
- `betab`: Path b effect (mediator → outcome)
- `betac`: Indirect effect (a × b)
- `beta_dir`: Direct effect
- `beta_all`: Total effect
- `prec_ci_lower/upper`: Bootstrap 95% CI for indirect effect
- `prec_lci_p/uci_p`: 95% CI for mediation proportion
- `p_value`: Empirical p-value from bootstrap

---

## Interpretation Guide

### MR Results

**Causal Effect:**
- **Beta/OR**: Magnitude of causal effect
- **95% CI**: Precision of estimate
- **P-value**: Statistical significance

**Validity Checks:**
- **Heterogeneity (I²)**: 
  - < 25%: Low heterogeneity
  - 25-50%: Moderate
  - 50-75%: Substantial
  - \> 75%: Considerable
- **MR-Egger intercept**: P > 0.05 suggests no directional pleiotropy
- **Leave-one-out**: Identifies influential SNPs

**Decision criteria:**
- Significant if: IVW P < threshold AND Heter.p > 0.05 AND Egger intercept P > 0.05

---

### Meta-Analysis Results

**Pooled estimate**: Combined effect across cohorts
**I² statistic**: Between-study heterogeneity
- I² < 50%: Acceptable consistency for fixed-effects model
- I² ≥ 50%: Consider random-effects or investigate sources

**Interpretation:**
- If 95% CI excludes null → significant pooled association
- Low I² + significant pooled effect → robust evidence

---

### Mediation Analysis Results

**Effects:**
- **Total effect (c)**: Overall exposure → outcome effect
- **Indirect effect (a×b)**: Effect mediated through the mediator
- **Direct effect (c')**: Exposure → outcome effect not through mediator

**Mediation proportion:**
- (Indirect effect / Total effect) × 100%
- Indicates percentage of total effect explained by mediator

**Significance:**
- Bootstrap CI not including zero → significant mediation
- P < 0.05 from empirical distribution

---

## Workflow Example

Complete analysis pipeline for a single protein:

```r
# Step 1: Discovery cohort MR
# (Set exposure = DPEP1, outcome = Graves' disease)
source("MR.R")
# Result: DPEP1 associated with GD (OR = 1.14, P = 2.3×10⁻⁵)

# Step 2: Replication cohort MR
# (Repeat with UKB-PPP data)
source("MR.R")
# Result: DPEP1 replicated (OR = 1.10, P = 0.003)

# Step 3: Meta-analysis
# (Combine results from both cohorts)
source("MetaAnalysis.R")
# Result: Pooled OR = 1.12 (95% CI: 1.08-1.16), I² = 15%

# Step 4: Forest plot
source("ForestPlotVisualization.R")
# Output: Visual comparison of both cohorts

# Step 5: Mediation analysis
# (Test if DNA methylation → DPEP1 → GD pathway exists)
exposure_name <- "DPEP1"
source("MediationAnalysis.R")
# Result: cg12345678 mediates 35% of total effect
```

---

## Quality Control

### Instrument Selection
- Genome-wide significance: P < 5×10⁻⁸
- LD clumping: r² < 0.1, window = 10 Mb
- F-statistic > 10 (avoid weak instrument bias)
- Cis-acting only (±1 Mb from gene)

### MR Assumptions
1. **Relevance**: Instruments strongly associated with exposure (✓ F > 10)
2. **Independence**: No confounding of instrument-outcome (✓ germline randomization)
3. **Exclusion restriction**: Instruments affect outcome only through exposure (✓ tested via pleiotropy)

### Sensitivity Thresholds
- Heterogeneity: P > 0.05
- Pleiotropy: MR-Egger intercept P > 0.05
- Outliers: MR-PRESSO detection and removal

---

## Project Structure

```
pQTL_GD/
├── MR.R                          # Core MR analysis
├── MetaAnalysis.R                # Meta-analysis across cohorts
├── ForestPlotVisualization.R     # Forest plot generation
├── MediationAnalysis.R           # Mediation analysis
├── README.md                     # This file
├── LICENSE                       # MIT License
├── data/                         # Input data (not tracked)
│   ├── exposure/
│   ├── outcome/
│   └── mediator/
├── results/                      # Analysis outputs
│   ├── discovery/
│   ├── replication/
│   └── meta/
└── figures/                      # Publication figures
    ├── forest_plots/
    ├── sensitivity/
    └── mediation/
```

---

## Citation

If you use these scripts in your research, please cite:

```bibtex
@article{luo2025methylation,
  title={DNA Methylation-Mediated Ferroptosis Pathway Promotes Graves' Disease Risk: A Multi-Omics Mendelian Randomization Study},
  author={Luo, Huanxian and [Co-authors]},
  journal={[Journal Name]},
  year={2025},
  note={In preparation}
}

@software{luo2025mrcode,
  author={Luo, Huanxian},
  title={Multi-Omics Mendelian Randomization Analysis Code},
  year={2025},
  url={https://github.com/lie2345/pQTL_GD}
}
```

---

## References

**Methodological:**

1. **TwoSampleMR**: Hemani G, et al. (2018). The MR-Base platform supports systematic causal inference across the human phenome. *eLife*, 7:e34408. https://doi.org/10.7554/eLife.34408

2. **MR-PRESSO**: Verbanck M, et al. (2018). Detection of widespread horizontal pleiotropy in causal relationships inferred from Mendelian randomization between complex traits and diseases. *Nat Genet*, 50:693-698. https://doi.org/10.1038/s41588-018-0099-7

3. **Mediation Analysis**: Relton CL & Davey Smith G (2012). Two-step epigenetic Mendelian randomization: a strategy for establishing the causal role of epigenetic processes in pathways to disease. *Int J Epidemiol*, 41:161-176. https://doi.org/10.1093/ije/dyr233

**Data Sources:**

4. **deCODE pQTL**: Ferkingstad E, et al. (2021). Large-scale integration of the plasma proteome with genetics and disease. *Nat Genet*, 53:1712-1721. https://doi.org/10.1038/s41588-021-00978-w

5. **UKB-PPP**: Sun BB, et al. (2023). Plasma proteomic associations with genetics and health in the UK Biobank. *Nature*, 622:329-338. https://doi.org/10.1038/s41586-023-06592-6

6. **GoDMC**: Min JL, et al. (2021). Genomic and phenotypic insights from an atlas of genetic effects on DNA methylation. *Nat Genet*, 53:1311-1321. https://doi.org/10.1038/s41588-021-00923-x

7. **FinnGen**: Kurki MI, et al. (2023). FinnGen provides genetic insights from a well-phenotyped isolated population. *Nature*, 613:508-518. https://doi.org/10.1038/s41586-022-05473-8

8. **FerrDb**: Zhou N & Bao J (2020). FerrDb: a manually curated resource for regulators and markers of ferroptosis and ferroptosis-disease associations. *Database*, 2020:baaa021. https://doi.org/10.1093/database/baaa021

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Contributing

Contributions are welcome! Please feel free to:
- Report bugs or issues
- Suggest new features
- Submit pull requests
- Improve documentation

**Guidelines:**
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/NewFeature`)
3. Commit changes (`git commit -m 'Add NewFeature'`)
4. Push to branch (`git push origin feature/NewFeature`)
5. Open a Pull Request

---

## Contact

**Author**: Huanxian Luo
**Email**: luohuanx@alumni.sysu.edu.cn
**Institution**: Sun Yat-sen University

**Project Repository**: https://github.com/lie2345/pQTL_GD

---

## Acknowledgments

- **Data providers**: FinnGen, UK Biobank, deCODE Genetics, GoDMC Consortium
- **Method developers**: MR-Base team, TwoSampleMR contributors
- **Funding**: [Add your funding sources if applicable]

---

## Version History

- **v1.0.0** (2025-01-15): Initial release
  - Core MR analysis
  - Meta-analysis functionality
  - Forest plot visualization
  - Mediation analysis

---

## Troubleshooting

**Common issues:**

1. **TwoSampleMR installation fails**
   ```r
   # Try alternative installation
   devtools::install_github("MRCIEU/TwoSampleMR")
   ```

2. **Memory errors with large datasets**
   ```r
   # Increase memory limit
   memory.limit(size = 16000)  # Windows
   ```

3. **LD clumping fails**
   - Ensure 1000 Genomes reference panel is accessible
   - Check TwoSampleMR version is up to date

4. **Forest plot text overlap**
   - Adjust x-coordinates in `geom_text()` calls
   - Increase plot width/height

For other issues, please open a GitHub issue with:
- R version and package versions
- Error message
- Minimal reproducible example

---

**Last updated**: January 15, 2025
