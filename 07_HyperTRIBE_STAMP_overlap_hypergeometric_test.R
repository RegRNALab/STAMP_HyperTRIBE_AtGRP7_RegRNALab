
#!/usr/bin/env Rscript

# ============================================================
# Hypergeometric analysis of HyperTRIBE–STAMP target-gene overlap
#
# This script tests whether the observed gene-level overlap
# between HyperTRIBE and STAMP candidate AtGRP7 targets is
# greater than expected by chance.
#
# Tissue-specific expressed genes are used as the background.
#
# Expressed gene criterion:
#   CPM >= 1 in at least 3 libraries
#
# Hypergeometric test:
#   Right-tailed test for enrichment
#   P(X >= observed overlap)
#
# Analyses are performed separately for:
#   1. Root
#   2. Cotyledon
#
# Input:
#   featureCounts gene-count matrix
#
# Output:
#   Tissue-specific expressed-gene backgrounds
#   Hypergeometric test results
# ============================================================


# ============================================================
# 1. Load packages
# ============================================================

library(data.table)
library(edgeR)
library(dplyr)


# ============================================================
# 2. Input and output paths
# ============================================================

count_file <-
  "/PATH/gene_counts_24samples.txt"

output_dir <-
  "/PATH/OUTPUT"


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. Read featureCounts output
# ============================================================

fc <- fread(
  count_file,
  data.table = FALSE,
  check.names = FALSE
)


cat("\n============================================\n")
cat("featureCounts file loaded\n")
cat("============================================\n")

cat(
  "Rows:",
  nrow(fc),
  "\n"
)

cat(
  "Columns:",
  ncol(fc),
  "\n\n"
)

cat("Column names:\n")

print(
  colnames(fc)
)


# ============================================================
# 4. Extract gene IDs
#
# Standard featureCounts output contains:
#
# Geneid
# Chr
# Start
# End
# Strand
# Length
# sample1
# sample2
# ...
# ============================================================

if (!"Geneid" %in% colnames(fc)) {
  stop(
    "Cannot find the 'Geneid' column. Check colnames(fc)."
  )
}


gene_id <- fc$Geneid


# ============================================================
# 5. Extract sample count columns
# ============================================================

annotation_columns <- c(
  "Geneid",
  "Chr",
  "Start",
  "End",
  "Strand",
  "Length"
)


sample_columns <- setdiff(
  colnames(fc),
  annotation_columns
)


cat("\n============================================\n")
cat("Detected sample columns\n")
cat("============================================\n")

print(
  sample_columns
)


cat(
  "\nNumber of libraries detected:",
  length(sample_columns),
  "\n"
)


# ============================================================
# 6. Create raw count matrix
# ============================================================

count_matrix <- as.matrix(
  fc[
    ,
    sample_columns,
    drop = FALSE
  ]
)


rownames(count_matrix) <- gene_id


storage.mode(
  count_matrix
) <- "numeric"


cat(
  "\nCount matrix dimensions:\n"
)

print(
  dim(count_matrix)
)


# ============================================================
# 7. Convert raw counts to CPM
# ============================================================

cpm_matrix <- edgeR::cpm(
  count_matrix
)


cat(
  "\nCPM matrix generated.\n"
)


# ============================================================
# 8. Identify root and cotyledon libraries
#
# Change the patterns below if sample names use different
# naming conventions.
# ============================================================


# ----------------------------
# ROOT
# ----------------------------

root_samples <- sample_columns[
  grepl(
    "Root|root|_R_",
    sample_columns
  )
]


# ----------------------------
# COTYLEDON
# ----------------------------

cotyledon_samples <- sample_columns[
  grepl(
    "Cotyledon|cotyledon|Shoot|shoot|_S_",
    sample_columns
  )
]


# ============================================================
# 9. Check detected libraries
# ============================================================

cat("\n============================================\n")
cat("ROOT LIBRARIES\n")
cat("============================================\n")

print(
  root_samples
)


cat(
  "\nNumber of root libraries:",
  length(root_samples),
  "\n"
)


cat("\n============================================\n")
cat("COTYLEDON LIBRARIES\n")
cat("============================================\n")

print(
  cotyledon_samples
)


cat(
  "\nNumber of cotyledon libraries:",
  length(cotyledon_samples),
  "\n"
)


# Safety checks
if (length(root_samples) == 0) {
  stop(
    "No root libraries detected. Check the sample-name pattern."
  )
}


if (length(cotyledon_samples) == 0) {
  stop(
    "No cotyledon libraries detected. Check the sample-name pattern."
  )
}


# ============================================================
# 10. Extract tissue-specific CPM matrices
# ============================================================

root_cpm <- cpm_matrix[
  ,
  root_samples,
  drop = FALSE
]


cotyledon_cpm <- cpm_matrix[
  ,
  cotyledon_samples,
  drop = FALSE
]


# ============================================================
# 11. Define expressed genes
#
# CPM >= 1 in at least 3 libraries
# ============================================================

root_expressed <- rowSums(
  root_cpm >= 1
) >= 3


cotyledon_expressed <- rowSums(
  cotyledon_cpm >= 1
) >= 3


# ============================================================
# 12. Calculate tissue-specific backgrounds
# ============================================================

N_root <- sum(
  root_expressed
)


N_cotyledon <- sum(
  cotyledon_expressed
)


cat("\n============================================\n")
cat("EXPRESSED-GENE BACKGROUNDS\n")
cat("CPM >= 1 in >= 3 libraries\n")
cat("============================================\n")


cat(
  "Root background:",
  format(
    N_root,
    big.mark = ","
  ),
  "genes\n"
)


cat(
  "Cotyledon background:",
  format(
    N_cotyledon,
    big.mark = ","
  ),
  "genes\n"
)


# ============================================================
# 13. Save expressed-gene backgrounds
# ============================================================

root_background_genes <- data.frame(
  Gene_ID = rownames(cpm_matrix)[
    root_expressed
  ]
)


cotyledon_background_genes <- data.frame(
  Gene_ID = rownames(cpm_matrix)[
    cotyledon_expressed
  ]
)


write.csv(
  root_background_genes,
  file.path(
    output_dir,
    "Root_expressed_gene_background.csv"
  ),
  row.names = FALSE
)


write.csv(
  cotyledon_background_genes,
  file.path(
    output_dir,
    "Cotyledon_expressed_gene_background.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. Hypergeometric test function
# ============================================================

run_hypergeometric <- function(
    tissue,
    background,
    HB_total,
    STAMP_total,
    overlap
) {
  
  # Expected overlap under random sampling
  expected <- (
    HB_total * STAMP_total
  ) / background
  
  
  # Fold enrichment
  fold_enrichment <-
    overlap / expected
  
  
  # Right-tailed hypergeometric P-value:
  # probability of observing >= the observed overlap
  p_value <- phyper(
    q = overlap - 1,
    m = HB_total,
    n = background - HB_total,
    k = STAMP_total,
    lower.tail = FALSE
  )
  
  
  data.frame(
    Tissue = tissue,
    Background = background,
    HyperTRIBE_targets = HB_total,
    STAMP_targets = STAMP_total,
    Observed_overlap = overlap,
    Expected_overlap = expected,
    Fold_enrichment = fold_enrichment,
    Hypergeometric_P = p_value
  )
}


# ============================================================
# 15. Root hypergeometric test
#
# HyperTRIBE = 7,661
# STAMP      = 8,234
# Shared     = 4,371
# ============================================================

root_result <- run_hypergeometric(
  tissue = "Root",
  background = N_root,
  HB_total = 7661,
  STAMP_total = 8234,
  overlap = 4371
)


# ============================================================
# 16. Cotyledon hypergeometric test
#
# HyperTRIBE = 6,852
# STAMP      = 7,616
# Shared     = 3,890
# ============================================================

cotyledon_result <- run_hypergeometric(
  tissue = "Cotyledon",
  background = N_cotyledon,
  HB_total = 6852,
  STAMP_total = 7616,
  overlap = 3890
)


# ============================================================
# 17. Combine results
# ============================================================

final_results <- bind_rows(
  root_result,
  cotyledon_result
)


# Keep full P-value while rounding expected overlap
# and fold enrichment for readability.

final_results$Expected_overlap <- round(
  final_results$Expected_overlap,
  1
)


final_results$Fold_enrichment <- round(
  final_results$Fold_enrichment,
  3
)


# ============================================================
# 18. Print final results
# ============================================================

cat("\n============================================\n")
cat("FINAL HYPERGEOMETRIC RESULTS\n")
cat("============================================\n\n")


print(
  final_results
)


# ============================================================
# 19. Save final table
# ============================================================

write.csv(
  final_results,
  file.path(
    output_dir,
    "HyperTRIBE_STAMP_hypergeometric_results.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 20. Final summary
# ============================================================

cat("\n============================================\n")
cat("Hypergeometric analysis finished\n")
cat("============================================\n")

cat(
  "Root background:",
  N_root,
  "\n"
)

cat(
  "Cotyledon background:",
  N_cotyledon,
  "\n"
)

cat(
  "Output directory:",
  output_dir,
  "\n"
)

cat("\nDONE.\n")