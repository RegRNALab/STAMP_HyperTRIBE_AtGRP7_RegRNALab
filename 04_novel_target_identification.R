

#!/usr/bin/env Rscript

# ============================================================
# Identification of novel AtGRP7 target genes and clusters
#
# This script identifies novel HyperTRIBE and STAMP candidate
# targets by excluding previously published AtGRP7 targets
# identified by iCLIP or RIP.
#
# Novel target:
#   detected by HyperTRIBE and/or STAMP
#   AND not present in the published iCLIP target set
#   AND not present in the published RIP target set
#
# Novel genes are further classified as:
#   1. HyperTRIBE + STAMP shared novel genes
#   2. HyperTRIBE-only novel genes
#   3. STAMP-only novel genes
#
# For cluster-level analyses, a strict filtering rule is used:
# if a cluster overlaps any known iCLIP/RIP target gene,
# the entire cluster is excluded from the novel cluster set.
#
# Inputs:
#   1. HyperTRIBE cluster-to-gene annotation
#   2. STAMP cluster-to-gene annotation
#   3. Published AtGRP7 iCLIP target genes
#   4. Published AtGRP7 RIP target genes
#
# Outputs:
#   Novel target gene lists
#   Novel cluster annotation tables
#   BED files for downstream analyses
#   Novel target summary table
# ============================================================


# ============================================================
# 1. Package
# ============================================================

library(data.table)


# ============================================================
# 2. Input files
# ============================================================

hb_file <-
  "/PATH/HB_clusters_gene_annotation.tsv"

stamp_file <-
  "/PATH/STAMP_clusters_gene_annotation.tsv"

iclip_file <-
  "/PATH/AtGRP7_iCLIP.txt"

rip_file <-
  "/PATH/RIP_AtGRP7.txt"

out_dir <-
  "/PATH/OUTPUT"


dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 3. Load HyperTRIBE and STAMP cluster-gene annotation
# ============================================================

hb <- fread(
  hb_file
)

stamp <- fread(
  stamp_file
)


# ============================================================
# 4. Basic check
# ============================================================

cat("\n========================================\n")
cat("ORIGINAL DATA\n")
cat("========================================\n")

cat(
  "HyperTRIBE annotation rows:",
  nrow(hb),
  "\n"
)

cat(
  "HyperTRIBE unique clusters:",
  uniqueN(hb$region_id),
  "\n"
)

cat(
  "HyperTRIBE unique genes:",
  uniqueN(
    hb[
      !is.na(gene_id) &
        gene_id != "",
      gene_id
    ]
  ),
  "\n\n"
)

cat(
  "STAMP annotation rows:",
  nrow(stamp),
  "\n"
)

cat(
  "STAMP unique clusters:",
  uniqueN(stamp$region_id),
  "\n"
)

cat(
  "STAMP unique genes:",
  uniqueN(
    stamp[
      !is.na(gene_id) &
        gene_id != "",
      gene_id
    ]
  ),
  "\n"
)


# ============================================================
# 5. Read iCLIP / RIP genes
# ============================================================

iclip_raw <- fread(
  iclip_file,
  header = FALSE,
  fill = TRUE
)

rip_raw <- fread(
  rip_file,
  header = FALSE,
  fill = TRUE
)


# Use first column
iclip_genes <- trimws(
  as.character(iclip_raw[[1]])
)

rip_genes <- trimws(
  as.character(rip_raw[[1]])
)


# Remove blanks / NA
iclip_genes <- unique(
  iclip_genes[
    !is.na(iclip_genes) &
      iclip_genes != ""
  ]
)

rip_genes <- unique(
  rip_genes[
    !is.na(rip_genes) &
      rip_genes != ""
  ]
)


# ============================================================
# 6. Clean Arabidopsis gene IDs
#
# Protect against version suffixes such as:
# AT1G01010.1
# ============================================================

iclip_genes <- sub(
  "\\.[0-9]+$",
  "",
  iclip_genes
)

rip_genes <- sub(
  "\\.[0-9]+$",
  "",
  rip_genes
)

hb[, gene_id_clean := sub(
  "\\.[0-9]+$",
  "",
  gene_id
)]

stamp[, gene_id_clean := sub(
  "\\.[0-9]+$",
  "",
  gene_id
)]


# ============================================================
# 7. Previously known AtGRP7 target genes
#
# known = iCLIP UNION RIP
# ============================================================

known_genes <- union(
  iclip_genes,
  rip_genes
)


cat("\n========================================\n")
cat("KNOWN AtGRP7 TARGET SET\n")
cat("========================================\n")

cat(
  "iCLIP genes:",
  length(iclip_genes),
  "\n"
)

cat(
  "RIP genes:",
  length(rip_genes),
  "\n"
)

cat(
  "iCLIP + RIP union:",
  length(known_genes),
  "\n"
)

cat(
  "Overlap between iCLIP and RIP:",
  length(
    intersect(
      iclip_genes,
      rip_genes
    )
  ),
  "\n"
)


# ============================================================
# 8. Check recovery in HyperTRIBE / STAMP
# ============================================================

hb_genes_original <- unique(
  hb[
    !is.na(gene_id_clean) &
      gene_id_clean != "",
    gene_id_clean
  ]
)

stamp_genes_original <- unique(
  stamp[
    !is.na(gene_id_clean) &
      gene_id_clean != "",
    gene_id_clean
  ]
)


cat("\n========================================\n")
cat("KNOWN TARGET RECOVERY\n")
cat("========================================\n")

cat(
  "HB genes overlapping iCLIP:",
  length(
    intersect(
      hb_genes_original,
      iclip_genes
    )
  ),
  "\n"
)

cat(
  "HB genes overlapping RIP:",
  length(
    intersect(
      hb_genes_original,
      rip_genes
    )
  ),
  "\n"
)

cat(
  "HB genes overlapping iCLIP OR RIP:",
  length(
    intersect(
      hb_genes_original,
      known_genes
    )
  ),
  "\n\n"
)

cat(
  "STAMP genes overlapping iCLIP:",
  length(
    intersect(
      stamp_genes_original,
      iclip_genes
    )
  ),
  "\n"
)

cat(
  "STAMP genes overlapping RIP:",
  length(
    intersect(
      stamp_genes_original,
      rip_genes
    )
  ),
  "\n"
)

cat(
  "STAMP genes overlapping iCLIP OR RIP:",
  length(
    intersect(
      stamp_genes_original,
      known_genes
    )
  ),
  "\n"
)


# ============================================================
# 9. Identify novel genes at gene level
#
# novel target gene =
# detected by HyperTRIBE/STAMP
# AND not present in iCLIP
# AND not present in RIP
# ============================================================

hb_novel_genes <- setdiff(
  hb_genes_original,
  known_genes
)

stamp_novel_genes <- setdiff(
  stamp_genes_original,
  known_genes
)


cat("\n========================================\n")
cat("NOVEL TARGET GENES\n")
cat("========================================\n")

cat(
  "HB novel genes:",
  length(hb_novel_genes),
  "\n"
)

cat(
  "STAMP novel genes:",
  length(stamp_novel_genes),
  "\n"
)


# ============================================================
# 10. Determine whether each cluster touches a known gene
#
# Important for multi-gene clusters.
#
# STRICT RULE:
# If ONE cluster overlaps ANY known iCLIP/RIP gene,
# that cluster is excluded from the novel cluster set.
# ============================================================

hb_cluster_status <- hb[
  ,
  .(
    n_annotation_rows = .N,
    
    n_genes = uniqueN(
      gene_id_clean[
        !is.na(gene_id_clean) &
          gene_id_clean != ""
      ]
    ),
    
    genes = paste(
      unique(
        gene_id_clean[
          !is.na(gene_id_clean) &
            gene_id_clean != ""
        ]
      ),
      collapse = ","
    ),
    
    any_known_gene = any(
      gene_id_clean %in% known_genes,
      na.rm = TRUE
    ),
    
    any_annotated_gene = any(
      !is.na(gene_id_clean) &
        gene_id_clean != ""
    )
  ),
  by = region_id
]


stamp_cluster_status <- stamp[
  ,
  .(
    n_annotation_rows = .N,
    
    n_genes = uniqueN(
      gene_id_clean[
        !is.na(gene_id_clean) &
          gene_id_clean != ""
      ]
    ),
    
    genes = paste(
      unique(
        gene_id_clean[
          !is.na(gene_id_clean) &
            gene_id_clean != ""
        ]
      ),
      collapse = ","
    ),
    
    any_known_gene = any(
      gene_id_clean %in% known_genes,
      na.rm = TRUE
    ),
    
    any_annotated_gene = any(
      !is.na(gene_id_clean) &
        gene_id_clean != ""
    )
  ),
  by = region_id
]


# ============================================================
# 11. Strict novel cluster IDs
#
# Must:
# 1. map to at least one gene
# 2. map to ZERO known iCLIP/RIP genes
# ============================================================

hb_novel_region_ids <- hb_cluster_status[
  any_annotated_gene == TRUE &
    any_known_gene == FALSE,
  region_id
]

stamp_novel_region_ids <- stamp_cluster_status[
  any_annotated_gene == TRUE &
    any_known_gene == FALSE,
  region_id
]


# ============================================================
# 12. Extract full annotation rows for novel clusters
# ============================================================

hb_novel <- hb[
  region_id %in% hb_novel_region_ids &
    !is.na(gene_id_clean) &
    gene_id_clean != "" &
    !(gene_id_clean %in% known_genes)
]

stamp_novel <- stamp[
  region_id %in% stamp_novel_region_ids &
    !is.na(gene_id_clean) &
    gene_id_clean != "" &
    !(gene_id_clean %in% known_genes)
]


# ============================================================
# 13. Recalculate novel gene sets from retained clusters
#
# This ensures gene lists and cluster lists are fully consistent.
# ============================================================

hb_novel_genes_final <- unique(
  hb_novel$gene_id_clean
)

stamp_novel_genes_final <- unique(
  stamp_novel$gene_id_clean
)


# ============================================================
# 14. HyperTRIBE / STAMP overlap among novel target genes
# ============================================================

shared_novel_genes <- intersect(
  hb_novel_genes_final,
  stamp_novel_genes_final
)

hb_only_novel_genes <- setdiff(
  hb_novel_genes_final,
  stamp_novel_genes_final
)

stamp_only_novel_genes <- setdiff(
  stamp_novel_genes_final,
  hb_novel_genes_final
)


cat("\n========================================\n")
cat("NOVEL TARGET GENE OVERLAP\n")
cat("========================================\n")

cat(
  "HB novel genes:",
  length(hb_novel_genes_final),
  "\n"
)

cat(
  "STAMP novel genes:",
  length(stamp_novel_genes_final),
  "\n"
)

cat(
  "Shared HB + STAMP novel genes:",
  length(shared_novel_genes),
  "\n"
)

cat(
  "HB-only novel genes:",
  length(hb_only_novel_genes),
  "\n"
)

cat(
  "STAMP-only novel genes:",
  length(stamp_only_novel_genes),
  "\n"
)


# ============================================================
# 15. Add novel class
# ============================================================

hb_novel[
  ,
  novel_gene_class := fifelse(
    gene_id_clean %in% shared_novel_genes,
    "HB_STAMP_shared_novel_gene",
    "HB_only_novel_gene"
  )
]

stamp_novel[
  ,
  novel_gene_class := fifelse(
    gene_id_clean %in% shared_novel_genes,
    "HB_STAMP_shared_novel_gene",
    "STAMP_only_novel_gene"
  )
]


# ============================================================
# 16. Count final unique novel clusters
# ============================================================

hb_novel_cluster_n <- uniqueN(
  hb_novel$region_id
)

stamp_novel_cluster_n <- uniqueN(
  stamp_novel$region_id
)


cat("\n========================================\n")
cat("FINAL NOVEL CLUSTERS\n")
cat("========================================\n")

cat(
  "HB novel clusters:",
  hb_novel_cluster_n,
  "\n"
)

cat(
  "STAMP novel clusters:",
  stamp_novel_cluster_n,
  "\n"
)


# ============================================================
# 17. Make one row per cluster for BED output
#
# This removes duplicated annotation rows caused by
# multi-gene clusters.
# ============================================================

hb_novel_bed <- unique(
  hb_novel[
    ,
    .(
      chrom,
      start,
      end,
      region_id,
      score = ".",
      strand
    )
  ],
  by = "region_id"
)


stamp_novel_bed <- unique(
  stamp_novel[
    ,
    .(
      chrom,
      start,
      end,
      region_id,
      score = ".",
      strand
    )
  ],
  by = "region_id"
)


# ============================================================
# 18. BED files by novel-gene class
# ============================================================

hb_shared_bed <- unique(
  hb_novel[
    gene_id_clean %in% shared_novel_genes,
    .(
      chrom,
      start,
      end,
      region_id,
      score = ".",
      strand
    )
  ],
  by = "region_id"
)

stamp_shared_bed <- unique(
  stamp_novel[
    gene_id_clean %in% shared_novel_genes,
    .(
      chrom,
      start,
      end,
      region_id,
      score = ".",
      strand
    )
  ],
  by = "region_id"
)


hb_only_bed <- unique(
  hb_novel[
    gene_id_clean %in% hb_only_novel_genes,
    .(
      chrom,
      start,
      end,
      region_id,
      score = ".",
      strand
    )
  ],
  by = "region_id"
)


stamp_only_bed <- unique(
  stamp_novel[
    gene_id_clean %in% stamp_only_novel_genes,
    .(
      chrom,
      start,
      end,
      region_id,
      score = ".",
      strand
    )
  ],
  by = "region_id"
)


# ============================================================
# 19. Save full annotation tables
# ============================================================

fwrite(
  hb_novel,
  file.path(
    out_dir,
    "HB_novel_target_clusters_annotation.tsv"
  ),
  sep = "\t"
)

fwrite(
  stamp_novel,
  file.path(
    out_dir,
    "STAMP_novel_target_clusters_annotation.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 20. Save gene lists
# ============================================================

fwrite(
  data.table(
    gene_id = sort(hb_novel_genes_final)
  ),
  file.path(
    out_dir,
    "HB_novel_target_genes.txt"
  ),
  col.names = FALSE
)

fwrite(
  data.table(
    gene_id = sort(stamp_novel_genes_final)
  ),
  file.path(
    out_dir,
    "STAMP_novel_target_genes.txt"
  ),
  col.names = FALSE
)

fwrite(
  data.table(
    gene_id = sort(shared_novel_genes)
  ),
  file.path(
    out_dir,
    "HB_STAMP_shared_novel_genes.txt"
  ),
  col.names = FALSE
)

fwrite(
  data.table(
    gene_id = sort(hb_only_novel_genes)
  ),
  file.path(
    out_dir,
    "HB_only_novel_genes.txt"
  ),
  col.names = FALSE
)

fwrite(
  data.table(
    gene_id = sort(stamp_only_novel_genes)
  ),
  file.path(
    out_dir,
    "STAMP_only_novel_genes.txt"
  ),
  col.names = FALSE
)


# ============================================================
# 21. Save BED files for downstream analyses
# ============================================================

fwrite(
  hb_novel_bed,
  file.path(
    out_dir,
    "HB_novel_target_clusters.bed"
  ),
  sep = "\t",
  col.names = FALSE
)

fwrite(
  stamp_novel_bed,
  file.path(
    out_dir,
    "STAMP_novel_target_clusters.bed"
  ),
  sep = "\t",
  col.names = FALSE
)


# HyperTRIBE-only novel genes
fwrite(
  hb_only_bed,
  file.path(
    out_dir,
    "HB_only_novel_target_clusters.bed"
  ),
  sep = "\t",
  col.names = FALSE
)


# STAMP-only novel genes
fwrite(
  stamp_only_bed,
  file.path(
    out_dir,
    "STAMP_only_novel_target_clusters.bed"
  ),
  sep = "\t",
  col.names = FALSE
)


# Shared novel genes detected by HyperTRIBE
fwrite(
  hb_shared_bed,
  file.path(
    out_dir,
    "HB_clusters_in_shared_novel_genes.bed"
  ),
  sep = "\t",
  col.names = FALSE
)


# Shared novel genes detected by STAMP
fwrite(
  stamp_shared_bed,
  file.path(
    out_dir,
    "STAMP_clusters_in_shared_novel_genes.bed"
  ),
  sep = "\t",
  col.names = FALSE
)


# ============================================================
# 22. Save known gene lists for record
# ============================================================

fwrite(
  data.table(
    gene_id = sort(iclip_genes)
  ),
  file.path(
    out_dir,
    "iCLIP_known_genes.txt"
  ),
  col.names = FALSE
)

fwrite(
  data.table(
    gene_id = sort(rip_genes)
  ),
  file.path(
    out_dir,
    "RIP_known_genes.txt"
  ),
  col.names = FALSE
)

fwrite(
  data.table(
    gene_id = sort(known_genes)
  ),
  file.path(
    out_dir,
    "iCLIP_RIP_combined_known_genes.txt"
  ),
  col.names = FALSE
)


# ============================================================
# 23. Final summary table
# ============================================================

final_summary <- data.table(
  
  category = c(
    "Original HB",
    "Original STAMP",
    "HB novel",
    "STAMP novel",
    "HB+STAMP shared novel genes",
    "HB-only novel genes",
    "STAMP-only novel genes"
  ),
  
  clusters = c(
    uniqueN(hb$region_id),
    uniqueN(stamp$region_id),
    uniqueN(hb_novel$region_id),
    uniqueN(stamp_novel$region_id),
    NA_integer_,
    uniqueN(hb_only_bed$region_id),
    uniqueN(stamp_only_bed$region_id)
  ),
  
  genes = c(
    length(hb_genes_original),
    length(stamp_genes_original),
    length(hb_novel_genes_final),
    length(stamp_novel_genes_final),
    length(shared_novel_genes),
    length(hb_only_novel_genes),
    length(stamp_only_novel_genes)
  )
)


fwrite(
  final_summary,
  file.path(
    out_dir,
    "novel_target_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 24. Final summary
# ============================================================

cat("\n========================================\n")
cat("FINAL SUMMARY\n")
cat("========================================\n")

print(final_summary)

cat("\nAll files saved to:\n")
cat(out_dir, "\n")

cat("\nDONE.\n")