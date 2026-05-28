# Naïve Bayes Classification of Skin Diseases Using Cytokine-Response Gene Signatures

A machine learning pipeline for multi-class disease classification using keratinocyte cytokine-response gene signatures derived from RNA sequencing data.

**Keywords:** Naïve Bayes | Gene Signatures | RNA-seq | Disease Classification | Leave-One-Out Cross-Validation

## Overview

This repository implements a Naïve Bayes classifier trained on cytokine-response gene signatures to discriminate between skin diseases. The approach uses gene expression profiles from primary human keratinocytes stimulated with various inflammatory cytokines to build a multi-class disease predictor with unbiased performance assessment via leave-one-out cross-validation.

## Key Features

- **Cytokine-response gene signatures:** Expression profiles from 10 cytokine stimulations
- **Efficient feature selection:** Top 15 genes per cytokine, ranked by statistical significance
- **Multi-class classification:** Simultaneous discrimination between 24+ diseases
- **Rigorous validation:** Leave-one-out cross-validation for unbiased assessment
- **Robust implementation:** Handles zero-variance genes and numerical stability
- **ROC curve analysis:** One-vs-rest evaluation for each disease with AUC quantification

## Data & Methods

### Cytokine Stimulation Protocol

Primary human keratinocytes were stimulated for 8 hours with:
- IL-4
- IL-13
- IL-17A
- IL-17A + TNF-α
- TNF-α
- IFN-α
- IFN-γ
- IL-36A
- IL-36B
- IL-36G

Unstimulated controls were processed in parallel for differential expression analysis.

### Gene Signature Derivation

1. **Differential expression analysis:** Identify cytokine-induced genes vs. unstimulated controls
2. **Selection criteria:** log2 fold-change ≥ 1
3. **Feature selection:** Top 15 genes per cytokine, ranked by adjusted p-value
4. **Combined signature:** ~150 genes total across all cytokines

### Naïve Bayes Classifier

**Model assumptions:**
- Conditional independence of gene expression given disease class
- Normal (Gaussian) distribution of gene expression within each class
- Uniform class priors: P(disease) = 1 / (number of diseases)

**Training:**
- For each disease and each gene: calculate mean and variance from training samples
- Exclude genes with zero variance in any disease class (prevents numerical errors)

**Prediction:**
- Compute posterior probability: P(disease | gene expression) using Bayes' theorem
- Use log-space calculations for numerical stability
- Normalize posterior probabilities to sum to 1

### Leave-One-Out Cross-Validation (LOOCV)

- Withhold one sample at a time during testing
- Retrain classifier on remaining n-1 samples
- Predict held-out sample's disease probability
- Repeat for all n samples
- Provides unbiased estimate of generalization performance

### Sample Filtering

- Only diseases with ≥5 samples included in analysis
- Rationale: Ensures adequate training data and prevents overfitting

## Project Structure

```
.
├── README.md                          # This file
├── Classifier.R                       # MA.classification.loo() function
├── analysis_example.R                 # Example workflow
├── figures/
│   └── roc_curves_skin.pdf           # ROC curves for all diseases
├── data/
│   ├── cytokine_signatures/          # Differential expression results per cytokine
│   ├── sample_metadata.csv           # Sample-disease annotations
│   └── expression_matrix.Rds         # Normalized expression data
├── results/
│   ├── classifier_predictions.csv    # Sample × disease posterior probabilities
│   ├── disease_rankings.csv          # Ranked predictions per sample
│   └── performance_metrics.txt       # AUC and accuracy statistics
└── LICENSE
```

## Usage

### Core Function: `MA.classification.loo()`

Located in `Classifier.R`

```r
res.naivebayes <- MA.classification.loo(
  tempinput.pheno,  # 2-column phenotype (sample names, disease labels)
  tempinput.exp,    # Expression matrix (genes × samples)
  tempinput.genes   # Character vector of signature genes
)
```

**Arguments:**
- `tempinput.pheno`: Data frame with 2 columns
  - Column 1: Sample names (character)
  - Column 2: Disease labels (character/factor)
- `tempinput.exp`: Numeric matrix
  - Rows: Gene names (should match `tempinput.genes`)
  - Columns: Sample names
  - Values: Log-normalized expression values
- `tempinput.genes`: Character vector of gene names to use as features

**Returns:**
- Numeric matrix of posterior probabilities
  - Rows: Sample names
  - Columns: Disease names
  - Values: Posterior probability (0–1) of each sample belonging to each disease

### Example Workflow

```r
# Load required libraries
library(tidyverse)
library(pROC)

# Source the classifier function
source("Classifier.R")

# Filter for diseases with ≥5 samples
s <- "skin"
tempp <- table(site.anno[[s]]$Disease)
tempdiseases <- names(tempp)[tempp >= 5]

# Extract filtered expression and phenotype data
temp.data <- site.DGEList_scale_voom_invnorm[[s]][,
  site.anno[[s]]$Disease %in% tempdiseases]
temp.pheno <- cbind(
  colnames(temp.data),
  as.character(site.anno[[s]]$Disease[site.anno[[s]]$Disease %in% tempdiseases])
)

# Select signature genes: top 15 per cytokine
tempcyt <- lapply(names(cytokineinduced), function(x) {
  data.frame(cytokineinduced[[x]]) %>%
    mutate(cytokine = x, gene = rownames(cytokineinduced[[x]]))
})

temp.cyt.upgenes <- tibble(as.data.frame(do.call(rbind, tempcyt))) %>%
  filter(lfcMLE >= 1) %>%
  arrange(pvalue) %>%
  group_by(cytokine) %>%
  slice(1:15)

# Run leave-one-out cross-validation classification
res.naivebayes <- MA.classification.loo(
  temp.pheno,
  temp.data,
  intersect(rownames(temp.data), unique(temp.cyt.upgenes$gene))
)

# Calculate disease rankings for each sample
res.naivebayes.rank <- t(apply(res.naivebayes, 1, function(x) rank(-x)))

# Determine rank of true disease for each sample
res.naivebayes.rank.truepheno <- character(dim(temp.pheno)[1])
for (i in 1:dim(temp.pheno)[1]) {
  res.naivebayes.rank.truepheno[i] <- res.naivebayes.rank[i, temp.pheno[i, 2]]
}

# Generate ROC curves for each disease
pdf("roc_curves_skin.pdf", width = 11, height = 11)
par(mfrow = c(5, 5))
for (d in colnames(res.naivebayes)) {
  true_labels <- as.numeric(temp.pheno[, 2] == d)
  plot.roc(true_labels, res.naivebayes[, d],
    percent = TRUE, main = d, print.auc = TRUE)
}
dev.off()
```

### Output Interpretation

**`res.naivebayes`** - Posterior Probabilities Matrix
```
                   Disease_A Disease_B Disease_C
Sample_001_disease_A  0.85      0.10      0.05
Sample_002_disease_B  0.15      0.75      0.10
Sample_003_disease_C  0.30      0.20      0.50
```
Each row represents a test sample; posterior probabilities sum to 1 across diseases.

**`res.naivebayes.rank`** - Disease Rankings Matrix
```
                   Disease_A Disease_B Disease_C
Sample_001_disease_A        1         2         3
Sample_002_disease_B        2         1         3
Sample_003_disease_C        3         2         1
```
Rank 1 = highest posterior probability for that sample.

**`res.naivebayes.rank.truepheno`** - True Disease Rank
```
Sample_001_disease_A: 1  (correctly ranked first)
Sample_002_disease_B: 1  (correctly ranked first)
Sample_003_disease_C: 1  (correctly ranked first)
```
Values of 1 indicate correct top-ranked classification.

## Technical Implementation Details

### Dynamic Gene Selection in LOOCV

During each cross-validation iteration:

1. **Calculate per-gene parameters** (mean, variance) for each disease using training samples only (n-1 samples)
2. **Identify problematic genes:** Genes with zero variance in any disease class
3. **Exclude problematic genes** from likelihood calculation (prevents division by zero and numerical instability)
4. **Predict held-out sample** using only the remaining genes with valid variance estimates

This adaptive approach ensures robustness across different training subsets.

### Numerical Stability

All posterior probability calculations are performed in log space to prevent numerical underflow:

```
log(posterior) = log(likelihood) + log(prior)
posterior = exp(log_posterior) / sum(exp(log_posterior))
```

This approach maintains numerical precision even with small probabilities across many genes and diseases.

## Requirements

- R ≥ 4.0
- tidyverse
- pROC (for ROC curve visualization)

### Installation

```r
install.packages(c("tidyverse", "pROC"))
```

If performing differential expression analysis to derive signatures:
```r
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
```

## Citation

If you use this classifier or methodology in your research, please cite:

[Your manuscript title]  
[Authors]  
[Journal] (Year)  
DOI: [your DOI]

## License

This project is licensed under the MIT License. See the LICENSE file for details.
