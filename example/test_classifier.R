
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

cat("Classification Performance:
")
cat("Correct predictions (rank #1):", sum(res.naivebayes.rank.truepheno == 1),
    "/", nrow(temp.pheno), "
")
cat("Accuracy:", round(100 * sum(res.naivebayes.rank.truepheno == 1) / nrow(temp.pheno), 2), "%
")

# Step 6: Generate ROC curves
pdf("test_roc_curves.pdf", width = 11, height = 11)
par(mfrow = c(3, 3))
for (disease in colnames(res.naivebayes)) {
  true_labels <- as.numeric(temp.pheno[, 2] == disease)
  plot.roc(true_labels, res.naivebayes[, disease],
    percent = TRUE, main = disease, print.auc = TRUE)
}
dev.off()

cat("ROC curves saved to test_roc_curves.pdf
")

