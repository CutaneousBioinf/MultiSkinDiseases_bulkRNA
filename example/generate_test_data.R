# Create Test Dataset for Naïve Bayes Classifier
# This script generates synthetic expression data with realistic disease patterns

# ========================================================================
# Load libraries
# ========================================================================
library(tidyverse)
set.seed(42)  # For reproducibility

# ========================================================================
# Parameters
# ========================================================================
n_samples <- 150                    # Total number of samples
n_diseases <- 8                     # Number of diseases
n_genes_total <- 2000              # Total number of genes in the matrix
n_signature_genes <- 150           # Number of signature genes (~15 per cytokine)
n_cytokines <- 10                  # Number of cytokines

# Disease distribution (ensuring ≥5 samples per disease)
samples_per_disease <- rep(floor(n_samples / n_diseases), n_diseases)
remainder <- n_samples - sum(samples_per_disease)
samples_per_disease[1:remainder] <- samples_per_disease[1:remainder] + 1

# ========================================================================
# Step 1: Create disease labels and sample metadata
# ========================================================================
disease_names <- paste0("Disease_", 1:n_diseases)
sample_names <- paste0("Sample_", sprintf("%03d", 1:n_samples))

# Create phenotype data
disease_labels <- rep(disease_names, samples_per_disease)
site.anno <- list(
  skin = data.frame(
    Sample_ID = sample_names,
    Disease = disease_labels,
    stringsAsFactors = FALSE
  )
)

cat("Sample distribution:\n")
print(table(site.anno$skin$Disease))

# ========================================================================
# Step 2: Create expression matrix with disease-specific patterns
# ========================================================================
# Initialize expression matrix with background noise
expr_matrix <- matrix(
  rnorm(n_genes_total * n_samples, mean = 0, sd = 1),
  nrow = n_genes_total,
  ncol = n_samples
)

gene_names <- paste0("GENE_", sprintf("%04d", 1:n_genes_total))
rownames(expr_matrix) <- gene_names
colnames(expr_matrix) <- sample_names

# ========================================================================
# Step 3: Add disease-specific signal to signature genes
# ========================================================================
# Select genes that will be signature genes
signature_gene_indices <- sample(1:n_genes_total, n_signature_genes, replace = FALSE)
signature_genes <- gene_names[signature_gene_indices]

# For each disease, create upregulated genes
genes_per_disease <- n_signature_genes / n_diseases

for (d in 1:n_diseases) {
  disease_name <- disease_names[d]
  
  # Select which genes are upregulated in this disease
  start_idx <- ((d - 1) * genes_per_disease) + 1
  end_idx <- min(d * genes_per_disease, n_signature_genes)
  disease_gene_indices <- signature_gene_indices[start_idx:end_idx]
  
  # Get sample indices for this disease
  disease_sample_indices <- which(site.anno$skin$Disease == disease_name)
  
  # Add disease-specific signal (higher mean expression for these genes in this disease)
  for (gene_idx in disease_gene_indices) {
    expr_matrix[gene_idx, disease_sample_indices] <- 
      expr_matrix[gene_idx, disease_sample_indices] + rnorm(
        length(disease_sample_indices),
        mean = 2.5,      # Mean increase in log-scale
        sd = 0.8
      )
  }
  
  cat("Disease", disease_name, "upregulates genes",
      gene_names[disease_gene_indices[1]], "to",
      gene_names[disease_gene_indices[length(disease_gene_indices)]], "\n")
}

# ========================================================================
# Step 4: Create differential expression results for each cytokine
# ========================================================================
# For each cytokine, select ~15 genes that are "induced"
# (These don't have to match disease-specific genes, but we'll make them related)

cytokineinduced <- list()
cytokine_names <- c("IL-4", "IL-13", "IL-17A", "IL-17A+TNF-alpha", "TNF-alpha",
                    "IFN-alpha", "IFN-gamma", "IL-36A", "IL-36B", "IL-36G")

for (cyt in cytokine_names) {
  # Select genes induced by this cytokine (some overlap with signature genes)
  n_induced <- sample(15:20, 1)  # 15-20 induced genes per cytokine
  
  # Bias towards signature genes but include some background genes
  induced_gene_indices <- c(
    sample(signature_gene_indices, min(10, length(signature_gene_indices)), replace = FALSE),
    sample(setdiff(1:n_genes_total, signature_gene_indices), 
           n_induced - 10, replace = FALSE)
  )
  induced_gene_indices <- induced_gene_indices[1:n_induced]
  
  # Create DEG table
  deg_table <- data.frame(
    baseMean = abs(rnorm(n_induced, mean = 5, sd = 2)),
    lfcMLE = rnorm(n_induced, mean = 1.5, sd = 0.5),  # Mean log2 FC ≥ 1
    pvalue = runif(n_induced, 0, 0.05),
    padj = runif(n_induced, 0, 0.05),
    row.names = gene_names[induced_gene_indices]
  )
  
  # Ensure lfcMLE ≥ 1 (to match filter criteria)
  deg_table$lfcMLE[deg_table$lfcMLE < 1] <- 1 + runif(sum(deg_table$lfcMLE < 1), 0, 0.5)
  
  cytokineinduced[[cyt]] <- deg_table
}

cat("\nCytokine-induced genes summary:\n")
for (cyt in names(cytokineinduced)) {
  cat(cyt, ":", nrow(cytokineinduced[[cyt]]), "genes\n")
}

# ========================================================================
# Step 5: Create normalized expression matrix (voom + inverse normal)
# ========================================================================
# Simulate voom transformation: log2(counts + 0.5)
expr_matrix_voom <- log2(pmax(expr_matrix, 0) + 0.5)

# Inverse normal transformation (rankbased)
expr_matrix_invnorm <- apply(expr_matrix_voom, 1, function(x) {
  qnorm((rank(x) - 0.5) / length(x))
})
expr_matrix_invnorm <- t(expr_matrix_invnorm)

# Store in list structure matching your data format
site.DGEList_scale_voom_invnorm <- list(
  skin = expr_matrix_invnorm
)

cat("\nExpression matrix dimensions:",
    nrow(site.DGEList_scale_voom_invnorm$skin), "genes ×",
    ncol(site.DGEList_scale_voom_invnorm$skin), "samples\n")

# ========================================================================
# Step 6: Verify data structure
# ========================================================================
cat("\n--- Data Structure Verification ---\n")
cat("site.anno$skin structure:\n")
print(head(site.anno$skin))

cat("\nsite.DGEList_scale_voom_invnorm$skin dimensions:\n")
print(dim(site.DGEList_scale_voom_invnorm$skin))
print("First few rows and columns:")
print(site.DGEList_scale_voom_invnorm$skin[1:5, 1:5])

cat("\ncytokineinduced[[1]] structure (IL-4):\n")
print(head(cytokineinduced[[1]]))

# ========================================================================
# Step 7: Save test data
# ========================================================================
# Save as RData file
save(site.anno, site.DGEList_scale_voom_invnorm, cytokineinduced,
     file = "test_data.RData")

cat("\n✓ Test data saved to 'test_data.RData'\n")

# ========================================================================
# Step 8: Create a simple test script for users
# ========================================================================
test_script <- '
# ========================================================================
# Test Script: Run classifier on test dataset
# ========================================================================

# Load test data
load("test_data.RData")

# Source classifier
source("Classifier.R")

# Load libraries
library(tidyverse)
library(pROC)

# Run classification workflow
tissue <- "skin"

# Step 1: Filter for diseases with >= 5 samples
disease_counts <- table(site.anno[[tissue]]$Disease)
diseases_filtered <- names(disease_counts)[disease_counts >= 5]

# Step 2: Extract data
sample_indices <- site.anno[[tissue]]$Disease %in% diseases_filtered
temp.data <- site.DGEList_scale_voom_invnorm[[tissue]][, sample_indices]
temp.pheno <- cbind(
  colnames(temp.data),
  as.character(site.anno[[tissue]]$Disease[sample_indices])
)

# Step 3: Select signature genes
tempcyt <- lapply(names(cytokineinduced), function(x) {
  data.frame(cytokineinduced[[x]]) %>%
    mutate(cytokine = x, gene = rownames(cytokineinduced[[x]]))
})

temp.cyt.upgenes <- tibble(as.data.frame(do.call(rbind, tempcyt))) %>%
  filter(lfcMLE >= 1) %>%
  arrange(pvalue) %>%
  group_by(cytokine) %>%
  slice(1:15)

signature_genes <- intersect(
  rownames(temp.data),
  unique(temp.cyt.upgenes$gene)
)

# Step 4: Run classifier
res.naivebayes <- MA.classification.loo(temp.pheno, temp.data, signature_genes)

# Step 5: Evaluate
res.naivebayes.rank <- t(apply(res.naivebayes, 1, function(x) rank(-x)))
res.naivebayes.rank.truepheno <- character(nrow(temp.pheno))
for (i in 1:nrow(temp.pheno)) {
  res.naivebayes.rank.truepheno[i] <- res.naivebayes.rank[i, temp.pheno[i, 2]]
}

cat("Classification Performance:\n")
cat("Correct predictions (rank #1):", sum(res.naivebayes.rank.truepheno == 1),
    "/", nrow(temp.pheno), "\n")
cat("Accuracy:", round(100 * sum(res.naivebayes.rank.truepheno == 1) / nrow(temp.pheno), 2), "%\n")

# Step 6: Generate ROC curves
pdf("test_roc_curves.pdf", width = 11, height = 11)
par(mfrow = c(3, 3))
for (disease in colnames(res.naivebayes)) {
  true_labels <- as.numeric(temp.pheno[, 2] == disease)
  plot.roc(true_labels, res.naivebayes[, disease],
    percent = TRUE, main = disease, print.auc = TRUE)
}
dev.off()

cat("ROC curves saved to test_roc_curves.pdf\n")
'

writeLines(test_script, "test_classifier.R")
cat("✓ Test script saved to 'test_classifier.R'\n")

# ========================================================================
# Summary
# ========================================================================
cat("\n========== SUMMARY ==========\n")
cat("Test dataset created successfully!\n\n")
cat("Files generated:\n")
cat("  1. test_data.RData - Contains site.anno, site.DGEList_scale_voom_invnorm, cytokineinduced\n")
cat("  2. test_classifier.R - Complete test workflow script\n\n")
cat("Dataset characteristics:\n")
cat("  - Samples:", n_samples, "\n")
cat("  - Diseases:", n_diseases, "with", min(samples_per_disease), "-", max(samples_per_disease), "samples each\n")
cat("  - Total genes:", n_genes_total, "\n")
cat("  - Signature genes:", n_signature_genes, "\n")
cat("  - Cytokines:", n_cytokines, "\n\n")
cat("Next steps:\n")
cat("  1. Run: source('test_classifier.R')\n")
cat("  2. Check output for classification accuracy\n")
cat("  3. Review generated ROC curves in test_roc_curves.pdf\n")
