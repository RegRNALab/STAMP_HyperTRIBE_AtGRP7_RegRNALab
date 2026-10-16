

#!/usr/bin/env Rscript

# ============================================================
# GO enrichment of method-specific novel AtGRP7 target genes
#
# This script performs Gene Ontology (GO) Biological Process
# enrichment analysis separately for:
#
#   1. HyperTRIBE-only novel target genes
#   2. STAMP-only novel target genes
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
#   1. HyperTRIBE-only novel target genes
#   2. STAMP-only novel target genes
#   3. Tissue-specific expressed gene background
#
# Outputs:
#   GO input gene lists
#   All GO BP enrichment results
#   Significant GO BP terms
#   Simplified significant GO BP terms
# ============================================================


# ============================================================
# 1. Packages
# ============================================================

library(data.table)
library(clusterProfiler)
library(org.At.tair.db)


# ============================================================
# 2. Input files
# ============================================================

HB_only_file <-
  "/PATH/HB_only_novel_genes.txt"

STAMP_only_file <-
  "/PATH/STAMP_only_novel_genes.txt"

cotyledon_background_file <-
  "/PATH/cotyledon_expressed_genes_CPM1_atleast3samples.txt"

out_dir <-
  "/PATH/OUTPUT"


dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 3. Read cotyledon-expressed background
# ============================================================

cotyledon_expressed_genes <- fread(
  cotyledon_background_file,
  header = FALSE
)[[1]]


cotyledon_expressed_genes <- unique(
  trimws(cotyledon_expressed_genes)
)


cotyledon_expressed_genes <-
  cotyledon_expressed_genes[
    cotyledon_expressed_genes != ""
  ]


cat(
  "Cotyledon background genes:",
  length(cotyledon_expressed_genes),
  "\n"
)


# ============================================================
# 4. Read HyperTRIBE-only novel genes
# ============================================================

HB_only_genes <- fread(
  HB_only_file,
  header = FALSE
)[[1]]


HB_only_genes <- unique(
  trimws(HB_only_genes)
)


HB_only_genes <- HB_only_genes[
  HB_only_genes != ""
]


cat(
  "HyperTRIBE-only novel genes:",
  length(HB_only_genes),
  "\n"
)


# ============================================================
# 5. Read STAMP-only novel genes
# ============================================================

STAMP_only_genes <- fread(
  STAMP_only_file,
  header = FALSE
)[[1]]


STAMP_only_genes <- unique(
  trimws(STAMP_only_genes)
)


STAMP_only_genes <- STAMP_only_genes[
  STAMP_only_genes != ""
]


cat(
  "STAMP-only novel genes:",
  length(STAMP_only_genes),
  "\n"
)


# ============================================================
# 6. Restrict foregrounds to cotyledon-expressed background
# ============================================================

HB_only_for_go <- intersect(
  HB_only_genes,
  cotyledon_expressed_genes
)


STAMP_only_for_go <- intersect(
  STAMP_only_genes,
  cotyledon_expressed_genes
)


cat(
  "HyperTRIBE-only genes used for GO:",
  length(HB_only_for_go),
  "\n"
)


cat(
  "STAMP-only genes used for GO:",
  length(STAMP_only_for_go),
  "\n"
)


# ============================================================
# 7. Save actual GO input genes
# ============================================================

write.table(
  HB_only_for_go,
  file.path(
    out_dir,
    "HB_only_novel_genes_used_for_GO.txt"
  ),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


write.table(
  STAMP_only_for_go,
  file.path(
    out_dir,
    "STAMP_only_novel_genes_used_for_GO.txt"
  ),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


# ============================================================
# 8. HyperTRIBE-only GO enrichment
# ============================================================

go_HB_only <- enrichGO(
  gene          = HB_only_for_go,
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
# 9. STAMP-only GO enrichment
# ============================================================

go_STAMP_only <- enrichGO(
  gene          = STAMP_only_for_go,
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
# 10. Convert GO results to data frames
# ============================================================

HB_df <- as.data.frame(
  go_HB_only
)

STAMP_df <- as.data.frame(
  go_STAMP_only
)


HB_df <- HB_df[
  order(HB_df$p.adjust),
]

STAMP_df <- STAMP_df[
  order(STAMP_df$p.adjust),
]


cat(
  "HyperTRIBE-only significant GO terms:",
  sum(
    HB_df$p.adjust < 0.05,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "STAMP-only significant GO terms:",
  sum(
    STAMP_df$p.adjust < 0.05,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# 11. Save all GO results
# ============================================================

write.csv(
  HB_df,
  file.path(
    out_dir,
    "HB_only_novel_GO_BP_cotyledonBG_ALL.csv"
  ),
  row.names = FALSE
)


write.csv(
  STAMP_df,
  file.path(
    out_dir,
    "STAMP_only_novel_GO_BP_cotyledonBG_ALL.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. Significant GO terms: adjusted P < 0.05
# ============================================================

HB_sig_df <- HB_df[
  !is.na(HB_df$p.adjust) &
    HB_df$p.adjust < 0.05,
]


STAMP_sig_df <- STAMP_df[
  !is.na(STAMP_df$p.adjust) &
    STAMP_df$p.adjust < 0.05,
]


write.csv(
  HB_sig_df,
  file.path(
    out_dir,
    "HB_only_novel_GO_BP_cotyledonBG_padj0.05.csv"
  ),
  row.names = FALSE
)


write.csv(
  STAMP_sig_df,
  file.path(
    out_dir,
    "STAMP_only_novel_GO_BP_cotyledonBG_padj0.05.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. Print top 30 significant GO terms
# ============================================================

cat(
  "\n========== COTYLEDON HYPERTRIBE ONLY ==========\n"
)


print(
  head(
    HB_sig_df[
      ,
      c(
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "p.adjust",
        "Count"
      )
    ],
    30
  )
)


cat(
  "\n========== COTYLEDON STAMP ONLY ==========\n"
)


print(
  head(
    STAMP_sig_df[
      ,
      c(
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "p.adjust",
        "Count"
      )
    ],
    30
  )
)


# ============================================================
# 14. Make significant enrichResult objects
# ============================================================

go_HB_only_sig <- go_HB_only


go_HB_only_sig@result <-
  go_HB_only_sig@result[
    !is.na(
      go_HB_only_sig@result$p.adjust
    ) &
      go_HB_only_sig@result$p.adjust < 0.05,
  ]


go_STAMP_only_sig <- go_STAMP_only


go_STAMP_only_sig@result <-
  go_STAMP_only_sig@result[
    !is.na(
      go_STAMP_only_sig@result$p.adjust
    ) &
      go_STAMP_only_sig@result$p.adjust < 0.05,
  ]


# ============================================================
# 15. Simplify redundant GO terms
# ============================================================

go_HB_only_simplified <- simplify(
  go_HB_only_sig,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)


go_STAMP_only_simplified <- simplify(
  go_STAMP_only_sig,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)


HB_simplified_df <- as.data.frame(
  go_HB_only_simplified
)


STAMP_simplified_df <- as.data.frame(
  go_STAMP_only_simplified
)


HB_simplified_df <- HB_simplified_df[
  order(HB_simplified_df$p.adjust),
]


STAMP_simplified_df <- STAMP_simplified_df[
  order(STAMP_simplified_df$p.adjust),
]


# ============================================================
# 16. Save simplified GO results
# ============================================================

write.csv(
  HB_simplified_df,
  file.path(
    out_dir,
    "HB_only_novel_GO_BP_cotyledonBG_SIMPLIFIED.csv"
  ),
  row.names = FALSE
)


write.csv(
  STAMP_simplified_df,
  file.path(
    out_dir,
    "STAMP_only_novel_GO_BP_cotyledonBG_SIMPLIFIED.csv"
  ),
  row.names = FALSE
)


cat(
  "HyperTRIBE-only simplified terms:",
  nrow(HB_simplified_df),
  "\n"
)


cat(
  "STAMP-only simplified terms:",
  nrow(STAMP_simplified_df),
  "\n"
)


# ============================================================
# 17. Final summary
# ============================================================

cat("\n")
cat("============================================\n")
cat("COTYLEDON GO comparison finished\n")
cat("============================================\n")


cat(
  "Cotyledon background:",
  length(cotyledon_expressed_genes),
  "\n"
)


cat(
  "HyperTRIBE-only novel genes:",
  length(HB_only_genes),
  "\n"
)


cat(
  "HyperTRIBE-only used for GO:",
  length(HB_only_for_go),
  "\n"
)


cat(
  "HyperTRIBE-only significant GO:",
  nrow(HB_sig_df),
  "\n"
)


cat(
  "HyperTRIBE-only simplified GO:",
  nrow(HB_simplified_df),
  "\n\n"
)


cat(
  "STAMP-only novel genes:",
  length(STAMP_only_genes),
  "\n"
)


cat(
  "STAMP-only used for GO:",
  length(STAMP_only_for_go),
  "\n"
)


cat(
  "STAMP-only significant GO:",
  nrow(STAMP_sig_df),
  "\n"
)


cat(
  "STAMP-only simplified GO:",
  nrow(STAMP_simplified_df),
  "\n"
)


cat(
  "Output directory:",
  out_dir,
  "\n"
)


cat("\nDONE.\n")