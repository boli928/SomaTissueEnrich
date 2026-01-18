# SomaTissueEnrich

**SomaEnrich** is an R package for performing tissue enrichment analysis specifically optimized for SomaScan data, using the Human Protein Atlas (HPA).

## Installation

You can install the development version from GitHub:

```r
install.packages("devtools")
devtools::install_github("boli928/SomaEnrich")
```

## Quick Start

```r
library(SomaEnrich)

# 1. Load your data
my_de_genes <- c("ALB", "CRP", "INS") # Your DE list
my_background <- c("ALB", "CRP", "INS", "TP53", "EGFR", ...) # Your full SomaScan menu

# 2. Run Enrichment
results <- run_tissue_enrichment(
    input_genes = my_de_genes, 
    background_genes = my_background
)

# 3. View Results
head(results)
