

#!/usr/bin/env Rscript

# ============================================================
# GO enrichment of shared novel AtGRP7 target genes
#
# This script performs Gene Ontology (GO) Biological Process
# enrichment analysis for shared novel AtGRP7 candidate targets
# identified by both HyperTRIBE and STAMP.
#
# The tissue-specific expressed gene set is used as the
# background universe.
#
# Expressed gene criterion:
#   CPM >= 1 in at least 3 libraries
#
# GO enrichment:
#   Ontology: Biological Process (BP)
#   Multiple-testing correction: Benjamini-Hochberg (BH)
#   Significant terms: adjusted P < 0.05
#
# Redundant significant GO terms are simplified using:
#   cutoff = 0.7
#   by = "p.adjust"
#   select_fun = min
#
# Inputs:
#   1. featureCounts gene-count matrix
#   2. Shared novel HyperTRIBE + STAMP target gene list
#
# Outputs:
#   Tissue-specific expressed gene background
#   Shared novel genes used for GO analysis
#   All GO BP enrichment results
#   Significant GO BP terms (adjusted P < 0.05)
#   Simplified significant GO BP terms
# ============================================================


# ============================================================
# 1. Packages
# ============================================================

library(data.table)
library(edgeR)
library(clusterProfiler)
library(org.At.tair.db)
library(AnnotationDbi)


# ============================================================
# 2. Input files
# ============================================================

# featureCounts matrix
count_file <-
  "/PATH/gene_counts_24samples.txt"

# Shared novel HyperTRIBE + STAMP genes
shared_file <-
  "/PATH/HB_STAMP_shared_novel_genes.txt"

# Output directory
out_dir <-
  "/PATH/OUTPUT"


dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 3. Read featureCounts matrix
# ============================================================

fc <- fread(
  count_file,
  skip = 1
)

cat(
  "featureCounts dimensions:",
  nrow(fc),
  "genes x",
  ncol(fc) - 6,
  "samples\n"
)


# ============================================================
# 4. Extract cotyledon samples
#
# featureCounts columns:
#   1-6   = gene annotation
#   7-18  = root samples
#   19-30 = cotyledon samples
#
# Sample order was confirmed before analysis.
# ============================================================

cotyledon_counts <- as.matrix(
  fc[, 19:30]
)

rownames(cotyledon_counts) <- fc$Geneid

cat(
  "Cotyledon count matrix:",
  nrow(cotyledon_counts),
  "genes x",
  ncol(cotyledon_counts),
  "samples\n"
)


cat("\nCotyledon samples used:\n")

print(
  colnames(cotyledon_counts)
)


# ============================================================
# 5. Calculate CPM
# ============================================================

cotyledon_cpm <- edgeR::cpm(
  cotyledon_counts
)

cat(
  "\nCPM matrix:",
  nrow(cotyledon_cpm),
  "genes x",
  ncol(cotyledon_cpm),
  "samples\n"
)


# ============================================================
# 6. Define cotyledon-expressed gene universe
#
# CPM >= 1 in at least 3 of the 12 cotyledon libraries
# ============================================================

keep_cotyledon <- rowSums(
  cotyledon_cpm >= 1
) >= 3


cotyledon_expressed_genes <-
  rownames(cotyledon_cpm)[
    keep_cotyledon
  ]


cat(
  "Cotyledon-expressed genes:",
  length(cotyledon_expressed_genes),
  "\n"
)


# Save background gene list
write.table(
  cotyledon_expressed_genes,
  file.path(
    out_dir,
    "cotyledon_expressed_genes_CPM1_atleast3samples.txt"
  ),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


# ============================================================
# 7. Read shared novel gene list
# ============================================================

shared_genes <- fread(
  shared_file,
  header = FALSE
)[[1]]


shared_genes <- unique(
  trimws(shared_genes)
)


shared_genes <- shared_genes[
  shared_genes != ""
]


cat(
  "Shared novel genes:",
  length(shared_genes),
  "\n"
)


# ============================================================
# 8. Restrict foreground to genes present in
#    cotyledon-expressed universe
# ============================================================

shared_for_go <- intersect(
  shared_genes,
  cotyledon_expressed_genes
)


cat(
  "Shared novel genes present in cotyledon background:",
  length(shared_for_go),
  "\n"
)


cat(
  "Shared novel genes excluded because not cotyledon-expressed:",
  length(
    setdiff(
      shared_genes,
      cotyledon_expressed_genes
    )
  ),
  "\n"
)


# Save genes actually used for GO analysis
write.table(
  shared_for_go,
  file.path(
    out_dir,
    "shared_novel_genes_used_for_GO.txt"
  ),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


# ============================================================
# 9. GO annotation coverage QC
#
# This does not change the enrichment analysis.
# It reports how many foreground genes have BP GO annotation.
# ============================================================

foreground_GO_annotation <- AnnotationDbi::select(
  org.At.tair.db,
  keys = shared_for_go,
  keytype = "TAIR",
  columns = c(
    "GO",
    "ONTOLOGY"
  )
)


foreground_BP_genes <- unique(
  foreground_GO_annotation$TAIR[
    !is.na(
      foreground_GO_annotation$GO
    ) &
      foreground_GO_annotation$ONTOLOGY == "BP"
  ]
)


cat(
  "Shared novel genes with BP GO annotation:",
  length(foreground_BP_genes),
  "\n"
)


cat(
  "Percentage of GO-input genes with BP annotation:",
  round(
    100 *
      length(foreground_BP_genes) /
      length(shared_for_go),
    2
  ),
  "%\n"
)


# ============================================================
# 10. GO Biological Process enrichment
#
# universe = cotyledon-expressed genes
# ============================================================

go_shared_cotyledonBG <- enrichGO(
  gene          = shared_for_go,
  universe      = cotyledon_expressed_genes,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,
  qvalueCutoff  = 1,
  readable      = FALSE
)


# ============================================================
# 11. Convert GO result to data frame
# ============================================================

go_all <- as.data.frame(
  go_shared_cotyledonBG
)


go_all <- go_all[
  order(
    go_all$p.adjust
  ),
]


cat(
  "Total GO BP terms tested:",
  nrow(go_all),
  "\n"
)


cat(
  "Significant GO BP terms, padj < 0.05:",
  sum(
    go_all$p.adjust < 0.05,
    na.rm = TRUE
  ),
  "\n"
)


# Save all GO results
write.csv(
  go_all,
  file.path(
    out_dir,
    "shared_novel_GO_BP_cotyledonBG_ALL.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. Save significant GO terms
# ============================================================

go_sig_df <- go_all[
  !is.na(
    go_all$p.adjust
  ) &
    go_all$p.adjust < 0.05,
]


write.csv(
  go_sig_df,
  file.path(
    out_dir,
    "shared_novel_GO_BP_cotyledonBG_padj0.05.csv"
  ),
  row.names = FALSE
)


cat(
  "Significant GO terms saved:",
  nrow(go_sig_df),
  "\n"
)


# ============================================================
# 13. Create significant enrichResult object
# ============================================================

go_shared_cotyledonBG_sig <-
  go_shared_cotyledonBG


go_shared_cotyledonBG_sig@result <-
  go_shared_cotyledonBG_sig@result[
    !is.na(
      go_shared_cotyledonBG_sig@result$p.adjust
    ) &
      go_shared_cotyledonBG_sig@result$p.adjust < 0.05,
  ]


cat(
  "GO terms before simplify:",
  nrow(
    go_shared_cotyledonBG_sig@result
  ),
  "\n"
)


# ============================================================
# 14. Simplify redundant GO terms
# ============================================================

go_shared_cotyledonBG_simplified <- simplify(
  go_shared_cotyledonBG_sig,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)


simplified_df <- as.data.frame(
  go_shared_cotyledonBG_simplified
)


simplified_df <- simplified_df[
  order(
    simplified_df$p.adjust
  ),
]


cat(
  "GO terms after simplify:",
  nrow(simplified_df),
  "\n"
)


# Save simplified GO table
write.csv(
  simplified_df,
  file.path(
    out_dir,
    "shared_novel_GO_BP_cotyledonBG_SIMPLIFIED.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. Summary
# ============================================================

cat("\n")
cat("============================================\n")
cat("GO analysis finished\n")
cat("============================================\n")


cat(
  "Cotyledon-expressed background genes:",
  length(cotyledon_expressed_genes),
  "\n"
)


cat(
  "Shared novel genes:",
  length(shared_genes),
  "\n"
)


cat(
  "Shared novel genes used for GO:",
  length(shared_for_go),
  "\n"
)


cat(
  "Shared novel genes with BP annotation:",
  length(foreground_BP_genes),
  "\n"
)


cat(
  "Significant GO terms:",
  nrow(go_sig_df),
  "\n"
)


cat(
  "Simplified GO terms:",
  nrow(simplified_df),
  "\n"
)


cat(
  "Output directory:",
  out_dir,
  "\n"
)


cat("\nDONE.\n")