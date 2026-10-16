
#!/usr/bin/env Rscript

# ============================================================
# Transcription factor annotation of shared novel AtGRP7 targets
#
# This script identifies transcription factors (TFs) among
# shared novel AtGRP7 candidate target genes detected by both
# HyperTRIBE and STAMP in roots and cotyledons.
#
# TF annotations are obtained from PlantTFDB.
#
# Analyses:
#   1. Identify TFs among root shared novel targets
#   2. Identify TFs among cotyledon shared novel targets
#   3. Summarize TF family representation
#   4. Identify TF candidates detected in both tissues
#   5. Identify tissue-associated TF candidates
#
# Inputs:
#   1. Root HyperTRIBE + STAMP shared novel genes
#   2. Cotyledon HyperTRIBE + STAMP shared novel genes
#   3. Arabidopsis thaliana TF annotation from PlantTFDB
#
# Outputs:
#   Root shared novel TFs
#   Cotyledon shared novel TFs
#   TF family counts
#   TFs detected in both tissues
#   Root-associated TF candidates
#   Cotyledon-associated TF candidates
# ============================================================


# ============================================================
# 1. Package
# ============================================================

library(data.table)


# ============================================================
# 2. Input files
# ============================================================

root_file <-
  "/PATH/Root_HB_STAMP_shared_novel_genes.txt"

cotyledon_file <-
  "/PATH/Cotyledon_HB_STAMP_shared_novel_genes.txt"

tf_file <-
  "/PATH/Ath_TF_list_downloaded_PlantTFDB.txt"

out_dir <-
  "/PATH/OUTPUT"


dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 3. Read shared novel genes
# ============================================================

root_novel <- fread(
  root_file,
  header = FALSE
)[[1]]

root_novel <- unique(
  trimws(root_novel)
)

root_novel <- root_novel[
  root_novel != ""
]


cotyledon_novel <- fread(
  cotyledon_file,
  header = FALSE
)[[1]]

cotyledon_novel <- unique(
  trimws(cotyledon_novel)
)

cotyledon_novel <- cotyledon_novel[
  cotyledon_novel != ""
]


cat(
  "Root shared novel genes:",
  length(root_novel),
  "\n"
)

cat(
  "Cotyledon shared novel genes:",
  length(cotyledon_novel),
  "\n"
)


# ============================================================
# 4. Read PlantTFDB annotation
# ============================================================

tf <- fread(
  tf_file
)


cat(
  "PlantTFDB rows:",
  nrow(tf),
  "\n"
)

cat(
  "Unique TF genes:",
  uniqueN(tf$Gene_ID),
  "\n"
)


# ============================================================
# 5. Collapse transcript isoforms to unique genes
# ============================================================

tf_gene <- unique(
  tf[
    ,
    .(
      Gene_ID,
      Family
    )
  ]
)


cat(
  "Unique Gene_ID + Family combinations:",
  nrow(tf_gene),
  "\n"
)


# Check whether any gene is assigned to more than one TF family
multi_family <- tf_gene[
  ,
  .N,
  by = Gene_ID
][
  N > 1
]


cat(
  "Genes assigned to multiple TF families:",
  nrow(multi_family),
  "\n"
)


# ============================================================
# 6. Identify TFs among root shared novel targets
# ============================================================

root_TFs <- tf_gene[
  Gene_ID %in% root_novel
]

root_TFs <- root_TFs[
  order(
    Family,
    Gene_ID
  )
]


# ============================================================
# 7. Identify TFs among cotyledon shared novel targets
# ============================================================

cotyledon_TFs <- tf_gene[
  Gene_ID %in% cotyledon_novel
]

cotyledon_TFs <- cotyledon_TFs[
  order(
    Family,
    Gene_ID
  )
]


# ============================================================
# 8. Summary
# ============================================================

cat("\n====================================\n")
cat("TF ANNOTATION SUMMARY\n")
cat("====================================\n")


cat(
  "Root shared novel targets:",
  length(root_novel),
  "\n"
)

cat(
  "Root TF targets:",
  uniqueN(root_TFs$Gene_ID),
  "\n"
)

cat(
  "Root TF proportion:",
  round(
    100 *
      uniqueN(root_TFs$Gene_ID) /
      length(root_novel),
    2
  ),
  "%\n\n"
)


cat(
  "Cotyledon shared novel targets:",
  length(cotyledon_novel),
  "\n"
)

cat(
  "Cotyledon TF targets:",
  uniqueN(cotyledon_TFs$Gene_ID),
  "\n"
)

cat(
  "Cotyledon TF proportion:",
  round(
    100 *
      uniqueN(cotyledon_TFs$Gene_ID) /
      length(cotyledon_novel),
    2
  ),
  "%\n"
)


# ============================================================
# 9. TF family counts
# ============================================================

root_family_counts <- root_TFs[
  ,
  .(
    n_TF_genes = uniqueN(Gene_ID)
  ),
  by = Family
][
  order(-n_TF_genes)
]


cotyledon_family_counts <- cotyledon_TFs[
  ,
  .(
    n_TF_genes = uniqueN(Gene_ID)
  ),
  by = Family
][
  order(-n_TF_genes)
]


cat(
  "\n========== ROOT TF FAMILIES ==========\n"
)

print(
  root_family_counts
)


cat(
  "\n========== COTYLEDON TF FAMILIES ==========\n"
)

print(
  cotyledon_family_counts
)


# ============================================================
# 10. TFs shared between root and cotyledon
# ============================================================

shared_tissue_TFs <- intersect(
  unique(root_TFs$Gene_ID),
  unique(cotyledon_TFs$Gene_ID)
)


cat(
  "\nTFs found in shared-novel populations of BOTH tissues:",
  length(shared_tissue_TFs),
  "\n"
)


# ============================================================
# 11. Tissue-associated TF candidates
# ============================================================

root_only_TFs <- setdiff(
  unique(root_TFs$Gene_ID),
  unique(cotyledon_TFs$Gene_ID)
)


cotyledon_only_TFs <- setdiff(
  unique(cotyledon_TFs$Gene_ID),
  unique(root_TFs$Gene_ID)
)


both_tissue_TFs <- intersect(
  unique(root_TFs$Gene_ID),
  unique(cotyledon_TFs$Gene_ID)
)


cat(
  "Root-associated TF candidates:",
  length(root_only_TFs),
  "\n"
)

cat(
  "Cotyledon-associated TF candidates:",
  length(cotyledon_only_TFs),
  "\n"
)

cat(
  "TF candidates detected in both tissues:",
  length(both_tissue_TFs),
  "\n"
)


# ============================================================
# 12. Save results
# ============================================================

fwrite(
  root_TFs,
  file.path(
    out_dir,
    "root_shared_novel_TFs.csv"
  )
)


fwrite(
  cotyledon_TFs,
  file.path(
    out_dir,
    "cotyledon_shared_novel_TFs.csv"
  )
)


fwrite(
  root_family_counts,
  file.path(
    out_dir,
    "root_shared_novel_TF_family_counts.csv"
  )
)


fwrite(
  cotyledon_family_counts,
  file.path(
    out_dir,
    "cotyledon_shared_novel_TF_family_counts.csv"
  )
)


fwrite(
  data.table(
    Gene_ID = shared_tissue_TFs
  ),
  file.path(
    out_dir,
    "shared_between_tissues_novel_TFs.csv"
  )
)


fwrite(
  data.table(
    Gene_ID = root_only_TFs
  ),
  file.path(
    out_dir,
    "root_only_shared_novel_TFs.csv"
  )
)


fwrite(
  data.table(
    Gene_ID = cotyledon_only_TFs
  ),
  file.path(
    out_dir,
    "cotyledon_only_shared_novel_TFs.csv"
  )
)


fwrite(
  data.table(
    Gene_ID = both_tissue_TFs
  ),
  file.path(
    out_dir,
    "both_tissues_shared_novel_TFs.csv"
  )
)


# ============================================================
# 13. Final summary
# ============================================================

cat("\n====================================\n")
cat("TF annotation finished\n")
cat("====================================\n")

cat(
  "Root TF targets:",
  uniqueN(root_TFs$Gene_ID),
  "\n"
)

cat(
  "Cotyledon TF targets:",
  uniqueN(cotyledon_TFs$Gene_ID),
  "\n"
)

cat(
  "TF candidates detected in both tissues:",
  length(both_tissue_TFs),
  "\n"
)

cat(
  "Output directory:",
  out_dir,
  "\n"
)

cat("\nDONE.\n")