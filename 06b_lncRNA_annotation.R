

#!/usr/bin/env Rscript

# ============================================================
# lncRNA annotation of shared novel AtGRP7 target genes
#
# This script identifies long non-coding RNAs (lncRNAs) among
# shared novel AtGRP7 candidate target genes detected by both
# HyperTRIBE and STAMP in roots and cotyledons.
#
# lncRNA annotation is based on the TAIR12 lncRNA reference.
#
# Analyses:
#   1. Identify lncRNAs among root shared novel targets
#   2. Identify lncRNAs among cotyledon shared novel targets
#   3. Identify lncRNA candidates detected in both tissues
#
# Inputs:
#   1. TAIR12 lncRNA gene annotation
#   2. Root HyperTRIBE + STAMP shared novel genes
#   3. Cotyledon HyperTRIBE + STAMP shared novel genes
#
# Outputs:
#   Root shared novel lncRNAs
#   Cotyledon shared novel lncRNAs
#   lncRNAs detected in both tissues
# ============================================================


# ============================================================
# 1. Package
# ============================================================

library(data.table)


# ============================================================
# 2. Input files
# ============================================================

# TAIR12 lncRNA reference
lncRNA_file <-
  "/PATH/TAIR12_lncRNA_gene_list.tsv"

# Root shared novel genes
root_file <-
  "/PATH/Root_HB_STAMP_shared_novel_genes.txt"

# Cotyledon shared novel genes
cotyledon_file <-
  "/PATH/Cotyledon_HB_STAMP_shared_novel_genes.txt"

# Output directory
output_dir <-
  "/PATH/OUTPUT"


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. Load TAIR12 lncRNA reference
# ============================================================

lncRNA_ref <- fread(
  lncRNA_file
)


lncRNA_ref[
  ,
  Gene_ID := toupper(
    trimws(Gene_ID)
  )
]


lncRNA_ref <- unique(
  lncRNA_ref,
  by = "Gene_ID"
)


cat(
  "TAIR12 annotated lncRNAs:",
  nrow(lncRNA_ref),
  "\n\n"
)


# ============================================================
# 4. Load shared novel gene lists
# ============================================================

root_shared <- fread(
  root_file,
  header = FALSE
)[[1]]


cotyledon_shared <- fread(
  cotyledon_file,
  header = FALSE
)[[1]]


root_shared <- unique(
  toupper(
    trimws(root_shared)
  )
)


cotyledon_shared <- unique(
  toupper(
    trimws(cotyledon_shared)
  )
)


cat(
  "Root shared novel genes:",
  length(root_shared),
  "\n"
)


cat(
  "Cotyledon shared novel genes:",
  length(cotyledon_shared),
  "\n\n"
)


# ============================================================
# 5. Overlap with TAIR12 lncRNA annotation
# ============================================================

root_lncRNA <- lncRNA_ref[
  Gene_ID %in% root_shared
]


cotyledon_lncRNA <- lncRNA_ref[
  Gene_ID %in% cotyledon_shared
]


# ============================================================
# 6. lncRNAs detected in both tissues
# ============================================================

shared_lncRNA <- intersect(
  root_lncRNA$Gene_ID,
  cotyledon_lncRNA$Gene_ID
)


# ============================================================
# 7. Summary
# ============================================================

cat(
  "========================================\n"
)

cat(
  "lncRNA ANNOTATION RESULTS\n"
)

cat(
  "========================================\n\n"
)


cat("ROOT\n")

cat(
  "Shared novel genes:",
  length(root_shared),
  "\n"
)

cat(
  "lncRNA candidates:",
  nrow(root_lncRNA),
  "\n"
)

cat(
  "Percentage:",
  round(
    100 *
      nrow(root_lncRNA) /
      length(root_shared),
    2
  ),
  "%\n\n"
)


cat("COTYLEDON\n")

cat(
  "Shared novel genes:",
  length(cotyledon_shared),
  "\n"
)

cat(
  "lncRNA candidates:",
  nrow(cotyledon_lncRNA),
  "\n"
)

cat(
  "Percentage:",
  round(
    100 *
      nrow(cotyledon_lncRNA) /
      length(cotyledon_shared),
    2
  ),
  "%\n\n"
)


cat(
  "lncRNA candidates detected in both tissues:",
  length(shared_lncRNA),
  "\n"
)


# ============================================================
# 8. Save results
# ============================================================

fwrite(
  root_lncRNA,
  file.path(
    output_dir,
    "Root_shared_novel_lncRNAs.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


fwrite(
  cotyledon_lncRNA,
  file.path(
    output_dir,
    "Cotyledon_shared_novel_lncRNAs.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


fwrite(
  data.table(
    Gene_ID = shared_lncRNA
  ),
  file.path(
    output_dir,
    "Root_Cotyledon_shared_lncRNAs.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 9. Final summary
# ============================================================

cat("\n========================================\n")
cat("lncRNA annotation finished\n")
cat("========================================\n")


cat(
  "Root lncRNA candidates:",
  nrow(root_lncRNA),
  "\n"
)


cat(
  "Cotyledon lncRNA candidates:",
  nrow(cotyledon_lncRNA),
  "\n"
)


cat(
  "lncRNA candidates detected in both tissues:",
  length(shared_lncRNA),
  "\n"
)


cat(
  "Output directory:",
  output_dir,
  "\n"
)


cat("\nDONE.\n")