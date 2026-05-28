# Naïve Bayes Classification of Skin Diseases Using Cytokine-Response Gene Signatures

A machine learning pipeline for multi-class disease classification using keratinocyte cytokine-response gene signatures derived from RNA sequencing data.

**Keywords:** Naïve Bayes | Gene Signatures | RNA-seq | Disease Classification | Leave-One-Out Cross-Validation

## Overview

This repository implements a Naïve Bayes classifier trained on cytokine-response gene signatures to discriminate between skin diseases. The approach uses gene expression profiles from primary human keratinocytes stimulated with various inflammatory cytokines to build a multi-class disease predictor with unbiased performance assessment via leave-one-out cross-validation.


## Data & Methods (data needs to processed according to the steps below)

### Gene Signature Derivation
e.g. use up-regulated genes from cytokine induced keratinocytes

1. **Differential expression analysis:** DESeq2 to compute summary statistics
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


### Data Structure

Your data should be organized as follows:

**`site.anno`** - List of phenotype annotations (one element per tissue type)
```r
site.anno[[tissue_name]]  # e.g., site.anno[["skin"]]
# Should be a data frame with columns:
# - Sample_ID (or equivalent)
# - Disease (disease/phenotype labels)
# - Any other clinical/metadata columns
```

**`site.DGEList_scale_voom_invnorm`** - List of expression matrices (one per tissue type)
```r
site.DGEList_scale_voom_invnorm[[tissue_name]]  # e.g., site.DGEList_scale_voom_invnorm[["skin"]]
# Should be a numeric matrix with:
# - Rows: Gene names
# - Columns: Sample names (matching site.anno[[tissue_name]])
# - Values: Log-normalized expression (voom-transformed, inverse normal transformed)
```

**`cytokineinduced`** - List of differential expression results (one per cytokine)
```r
cytokineinduced[["IL-4"]]   # Example for one cytokine
cytokineinduced[["IL-13"]]
cytokineinduced[["IL-17A"]]
# ... and so on for all 10 cytokines
# Each element should be a matrix or data frame with columns:
# - baseMean (or equivalent expression level)
# - lfcMLE (log2 fold-change)
# - pvalue (p-value)
# - padj (adjusted p-value)
# - Rownames: Gene names
```

### Complete Example Workflow

```r
# ========================================================================
# Load required libraries
# ========================================================================
library(tidyverse)
library(pROC)

# ========================================================================
# Source the classifier function
# ========================================================================
source("Classifier.R")

# ========================================================================
# Define data objects (you must load/create these before proceeding)
# ========================================================================
# site.anno: List of phenotype annotations per tissue
#   structure: site.anno[["tissue_name"]]$Disease contains disease labels
# site.DGEList_scale_voom_invnorm: List of expression matrices per tissue
#   structure: site.DGEList_scale_voom_invnorm[["tissue_name"]] is genes × samples matrix
# cytokineinduced: List of differential expression results per cytokine
#   structure: cytokineinduced[[cytokine_name]] with columns: baseMean, lfcMLE, pvalue, padj
#   rownames: Gene names

# Example of data loading (modify paths/objects to match your setup):
# load("path/to/site.anno.RData")
# load("path/to/site.DGEList_scale_voom_invnorm.RData")
# load("path/to/cytokineinduced.RData")

# ========================================================================
# Step 1: Select tissue type
# ========================================================================
tissue <- "skin"  # Change to "oral", "joint", etc. as needed

# ========================================================================
# Step 2: Filter for diseases with ≥5 samples
# ========================================================================
disease_counts <- table(site.anno[[tissue]]$Disease)
diseases_filtered <- names(disease_counts)[disease_counts >= 5]

cat("Number of diseases with ≥5 samples:", length(diseases_filtered), "\n")
cat("Diseases included:", paste(diseases_filtered, collapse = ", "), "\n")

# ========================================================================
# Step 3: Extract filtered expression and phenotype data
# ========================================================================
# Get column indices matching filtered diseases
sample_indices <- site.anno[[tissue]]$Disease %in% diseases_filtered

# Extract expression matrix and phenotype for filtered samples
temp.data <- site.DGEList_scale_voom_invnorm[[tissue]][, sample_indices]
temp.pheno <- cbind(
  colnames(temp.data),
  as.character(site.anno[[tissue]]$Disease[sample_indices])
)
colnames(temp.pheno) <- c("Sample_ID", "Disease")

cat("Expression matrix dimensions:", nrow(temp.data), "genes ×", 
    ncol(temp.data), "samples\n")

# ========================================================================
# Step 4: Select signature genes - top 15 per cytokine
# ========================================================================
# Reorganize cytokine DEG results into a single tibble
tempcyt <- lapply(names(cytokineinduced), function(cytokine_name) {
  deg_matrix <- cytokineinduced[[cytokine_name]]
  data.frame(deg_matrix) %>%
    mutate(
      cytokine = cytokine_name,
      gene = rownames(deg_matrix)
    )
})

# Filter and select genes
temp.cyt.upgenes <- tibble(as.data.frame(do.call(rbind, tempcyt))) %>%
  filter(lfcMLE >= 1) %>%           # log2 fold-change ≥ 1
  arrange(pvalue) %>%               # rank by p-value (most significant first)
  group_by(cytokine) %>%
  slice(1:15)                        # top 15 per cytokine

# Get unique gene names that are present in expression matrix
signature_genes <- intersect(
  rownames(temp.data),
  unique(temp.cyt.upgenes$gene)
)

cat("Number of signature genes:", length(signature_genes), "\n")
cat("Genes per cytokine (target: 15):\n")
print(table(temp.cyt.upgenes$cytokine))

# ========================================================================
# Step 5: Run leave-one-out cross-validation classification
# ========================================================================
res.naivebayes <- MA.classification.loo(
  temp.pheno,
  temp.data,
  signature_genes
)

cat("Classification complete. Results dimensions:", 
    nrow(res.naivebayes), "samples ×", 
    ncol(res.naivebayes), "diseases\n")

# ========================================================================
# Step 6: Calculate disease rankings for each sample
# ========================================================================
# For each sample, rank diseases by posterior probability (1 = highest)
res.naivebayes.rank <- t(apply(res.naivebayes, 1, function(x) rank(-x)))

# Get the rank of the true disease for each sample
res.naivebayes.rank.truepheno <- character(nrow(temp.pheno))
for (i in 1:nrow(temp.pheno)) {
  true_disease <- temp.pheno[i, 2]
  res.naivebayes.rank.truepheno[i] <- res.naivebayes.rank[i, true_disease]
}

# Summary statistics
cat("\nClassification performance summary:\n")
cat("Samples with true disease ranked #1:",
    sum(res.naivebayes.rank.truepheno == 1), "/", nrow(temp.pheno), "\n")
cat("Mean rank of true disease:",
    round(mean(as.numeric(res.naivebayes.rank.truepheno)), 2), "\n")

# ========================================================================
# Step 7: Generate ROC curves for each disease
# ========================================================================
pdf("roc_curves_skin.pdf", width = 11, height = 11)
par(mfrow = c(5, 5))

for (disease in colnames(res.naivebayes)) {
  # Create binary labels: 1 = this disease, 0 = all other diseases
  true_labels <- as.numeric(temp.pheno[, 2] == disease)
  
  # Get predicted probabilities for this disease
  predictions <- res.naivebayes[, disease]
  
  # Plot ROC curve
  roc_obj <- plot.roc(true_labels, predictions,
    percent = TRUE,
    main = disease,
    print.auc = TRUE
  )
}
dev.off()

cat("ROC curves saved to: roc_curves_skin.pdf\n")

# ========================================================================
# Step 8: Save results
# ========================================================================
# Save posterior probabilities
write.csv(res.naivebayes, "classifier_predictions.csv")

# Save disease rankings
write.csv(res.naivebayes.rank, "disease_rankings.csv")

# Save true disease ranks
results_summary <- data.frame(
  Sample = rownames(res.naivebayes),
  True_Disease = temp.pheno[, 2],
  True_Disease_Rank = as.numeric(res.naivebayes.rank.truepheno)
)
write.csv(results_summary, "classification_summary.csv", row.names = FALSE)

cat("Results saved successfully!\n")
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


### Disclaimer
While we conducted vigorous normalization, leave one out cross validation, on the data preparation and processing, etc, we acknowledge there can be bias caused by confounding factors including batch effect, biopsy sampling, and the nature of bulk RNA-seq. 
