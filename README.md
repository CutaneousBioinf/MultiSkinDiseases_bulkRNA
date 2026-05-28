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
- **ROC curve analysis:** One-vs-rest evaluation for each disease

## Data & Methods

### Cytokine Stimulation Protocol

Primary human keratinocytes stimulated for 8 hours with:
- IL-4, IL-13, IL-17A, IL-17A+TNF-α, TNF-α, IFN-α, IFN-γ, IL-36A, IL-36B, IL-36G

### Gene Signature Derivation

1. Differential expression analysis (stimulated vs. unstimulated)
2. Selection: log2 fold-change ≥ 1
3. Feature selection: Top 15 genes per cytokine by p-value
4. Combined signature: ~150 genes

### Naïve Bayes Classifier

**Model assumptions:**
- Conditional independence of genes given disease
- Normal distribution of gene expression per class
- Uniform class priors

**Training:**
- Calculate mean and variance per gene per disease
- Exclude genes with zero variance (numerical stability)

**Prediction:**
- Compute posterior probability: P(disease | expression)
- Use log-space for numerical stability

### Leave-One-Out Cross-Validation (LOOCV)

- Withhold one sample at a time
- Retrain classifier on n-1 samples
- Predict held-out sample
- Repeat for all n samples
- Provides unbiased generalization estimate

### Sample Filtering

Only diseases with ≥5 samples included (prevents overfitting).

## Project Structure

```
.
├── README.md                          # This file
├── naïve_bayes_classifier.R          # Main classification function
├── analysis_example.R                 # Example workflow
├── figures/
│   └── roc_curves_skin.pdf           # ROC curves for all diseases
├── data/
│   ├── cytokine_signatures/          # Differential expression results
│   ├── sample_metadata.csv           # Sample-disease annotations
│   └── expression_matrix.Rds         # Normalized expression data
├── results/
│   ├── classifier_predictions.csv    # Posterior probabilities
│   ├── disease_rankings.csv          # Ranked predictions
│   └── performance_metrics.txt       # AUC and accuracy
└── LICENSE
```

## Usage

### Core Function: `MA.classification.loo()`

```r
res.naivebayes <- MA.classification.loo(
  tempinput.pheno,  # 2-column phenotype (sample names, disease labels)
  tempinput.exp,    # Expression matrix (genes × samples)
  tempinput.genes   # Character vector of signature genes
)
```

**Arguments:**
- `tempinput.pheno`: Data frame with 2 columns (sample names, disease labels)
- `tempinput.exp`: Expression matrix (log-normalized values)
- `tempinput.genes`: Gene names to use as features

**Returns:**
- Matrix of posterior probabilities (samples × diseases)

### Example Workflow

```r
library(tidyverse)
library(pROC)
source("naïve_bayes_classifier.R")

# Filter for diseases with ≥5 samples
s <- "skin"
tempp <- table(site.anno[[s]]$Disease)
tempdiseases <- names(tempp)[tempp >= 5]

# Extract data
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

# Run LOOCV classifier
res.naivebayes <- MA.classification.loo(
  temp.pheno,
  temp.data,
  intersect(rownames(temp.data), unique(temp.cyt.upgenes$gene))
)

# Calculate disease rankings
res.naivebayes.rank <- t(apply(res.naivebayes, 1, function(x) rank(-x)))

# Get rank of true disease
res.naivebayes.rank.truepheno <- character(dim(temp.pheno)[1])
for (i in 1:dim(temp.pheno)[1]) {
  res.naivebayes.rank.truepheno[i] <- res.naivebayes.rank[i, temp.pheno[i, 2]]
}

# Generate ROC curves
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

**Posterior probabilities** - sample × disease matrix where each row sums to 1
**Disease rankings** - rank 1 = highest probability for that sample
**True disease rank** - should be 1 for correctly classified samples

## Technical Implementation

### Gene Selection During LOOCV

Each iteration:
1. Calculate mean and variance per gene per disease (training set only)
2. Identify genes with zero variance in any class
3. Exclude problematic genes from prediction
4. Predict held-out sample using remaining genes

### Numerical Stability

Computations in log space prevent underflow:
```
log(posterior) = log(likelihood) + log(prior)
posterior = exp(log_posterior) / sum(exp(log_posterior))
```

## Requirements

- R ≥ 4.0
- tidyverse
- pROC

Install:
```r
install.packages(c("tidyverse", "pROC"))
```

## Citation

[Your manuscript citation]

## Authors

Alex Tsoi, PhD  
University of Michigan School of Public Health

## License

MIT License - see LICENSE file

## Contact

For questions: alextsoi@umich.edu

