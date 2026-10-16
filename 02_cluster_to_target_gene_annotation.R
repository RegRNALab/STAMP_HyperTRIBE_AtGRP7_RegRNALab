


#!/usr/bin/env Rscript

# ============================================================
# Gene annotation of high-confidence editing clusters
# Cotyledon samples
#
# Inputs:
#   1. HyperTRIBE high-confidence clusters
#   2. STAMP high-confidence clusters
#   3. TAIR12 genome annotation (GTF/GFF3)
#
# Outputs:
#   HyperTRIBE:
#     HB_clusters_gene_annotation.tsv
#     HB_unique_genes.tsv
#
#   STAMP:
#     STAMP_clusters_gene_annotation.tsv
#     STAMP_unique_genes.tsv
#
#   Summary:
#     HB_STAMP_gene_annotation_summary.tsv
#     HB_STAMP_gene_overlap_summary.tsv
#
#   Gene lists:
#     HB_STAMP_shared_genes.tsv
#     HB_only_genes.tsv
#     STAMP_only_genes.tsv
#
# IMPORTANT:
# Genome annotation must correspond to the same genome assembly
# used for read alignment and FLARE analysis.
# ============================================================


# ============================================================
# 1. Packages
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(GenomicRanges)
  library(rtracklayer)
})


# ============================================================
# 2. Input files
# ============================================================

# User-defined input paths

hb_file <-
  "/PATH/HB_cotyledon_high_confidence_clusters/HB_cotyledon_high_confidence_clusters.tsv"

stamp_file <-
  "/PATH/STAMP_cotyledon_high_confidence_clusters/STAMP_cotyledon_high_confidence_clusters.tsv"


# Genome annotation corresponding to the same genome assembly
# used for read alignment and FLARE analysis

annotation_file <-
  "/PATH/TAIR12_cleaned.gtf"


# Output directory

output_dir <-
  "/PATH/OUTPUT"


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 3. Check files
# ============================================================

for (f in c(
  hb_file,
  stamp_file,
  annotation_file
)) {
  
  if (!file.exists(f)) {
    stop("File does not exist: ", f)
  }
}


# ============================================================
# 4. Read cluster files
# ============================================================

hb <- read_tsv(
  hb_file,
  show_col_types = FALSE
)

stamp <- read_tsv(
  stamp_file,
  show_col_types = FALSE
)


required_cluster_columns <- c(
  "region_id",
  "chrom",
  "start",
  "end",
  "strand"
)


if (!all(required_cluster_columns %in% colnames(hb))) {
  stop("Required columns missing from HyperTRIBE cluster file.")
}


if (!all(required_cluster_columns %in% colnames(stamp))) {
  stop("Required columns missing from STAMP cluster file.")
}


cat("\nInput clusters:\n")
cat("HyperTRIBE:", nrow(hb), "\n")
cat("STAMP:", nrow(stamp), "\n")


# ============================================================
# 5. Read genome annotation
# ============================================================

message("Reading genome annotation...")

annotation <- import(
  annotation_file
)


cat(
  "\nAnnotation feature types:\n"
)

print(
  sort(
    table(annotation$type),
    decreasing = TRUE
  )
)


# ============================================================
# 6. Keep gene features
# ============================================================

genes <- annotation[
  annotation$type == "gene"
]


if (length(genes) == 0) {
  stop(
    "No features with type='gene' found in annotation file."
  )
}


cat(
  "\nNumber of gene features:",
  length(genes),
  "\n"
)


# ============================================================
# 7. Inspect available metadata columns
# ============================================================

cat(
  "\nGene annotation metadata columns:\n"
)

print(
  colnames(
    mcols(genes)
  )
)


# ============================================================
# 8. Standardise gene information
#
# Different GFF files may use:
#
#   gene_id
#   ID
#   Name
#   gene
#   locus_tag
#
# Therefore this section tries several possibilities.
# ============================================================

gene_meta <- as.data.frame(
  mcols(genes)
)


get_first_existing <- function(
    df,
    candidates
) {
  
  available <- candidates[
    candidates %in% colnames(df)
  ]
  
  if (length(available) == 0) {
    return(
      rep(
        NA_character_,
        nrow(df)
      )
    )
  }
  
  as.character(
    df[[available[1]]]
  )
}


genes$gene_id_standard <- get_first_existing(
  gene_meta,
  c(
    "gene_id",
    "ID",
    "gene",
    "locus_tag",
    "Name"
  )
)


genes$gene_name_standard <- get_first_existing(
  gene_meta,
  c(
    "gene_name",
    "Name",
    "gene",
    "locus_tag",
    "gene_id",
    "ID"
  )
)


# ============================================================
# 9. Function: annotate one cluster dataset
# ============================================================

annotate_clusters <- function(
    cluster_df,
    method_name
) {
  
  message(
    "Annotating ",
    method_name,
    " clusters..."
  )
  
  
  # ----------------------------------------------------------
  # Convert clusters to genomic ranges
  # ----------------------------------------------------------
  
  cluster_gr <- GRanges(
    
    seqnames =
      cluster_df$chrom,
    
    ranges =
      IRanges(
        start = cluster_df$start,
        end = cluster_df$end
      ),
    
    strand =
      cluster_df$strand
  )
  
  
  cluster_gr$cluster_row <-
    seq_len(
      nrow(cluster_df)
    )
  
  
  # ----------------------------------------------------------
  # Strand-specific overlap
  #
  # ignore.strand = FALSE:
  # cluster must overlap gene on matching strand.
  # ----------------------------------------------------------
  
  hits <- findOverlaps(
    cluster_gr,
    genes,
    ignore.strand = FALSE
  )
  
  
  # ----------------------------------------------------------
  # No overlap case
  # ----------------------------------------------------------
  
  if (length(hits) == 0) {
    
    warning(
      "No gene overlaps detected for ",
      method_name,
      ". Check chromosome naming / assembly."
    )
    
    annotated <- cluster_df %>%
      mutate(
        gene_id = NA_character_,
        gene_name = NA_character_,
        gene_chrom = NA_character_,
        gene_start = NA_integer_,
        gene_end = NA_integer_,
        gene_strand = NA_character_
      )
    
    return(annotated)
  }
  
  
  # ----------------------------------------------------------
  # Build cluster -> gene table
  # ----------------------------------------------------------
  
  hit_df <- tibble(
    
    cluster_row =
      queryHits(hits),
    
    gene_row =
      subjectHits(hits)
  ) %>%
    
    mutate(
      
      gene_id =
        as.character(
          genes$gene_id_standard[
            gene_row
          ]
        ),
      
      gene_name =
        as.character(
          genes$gene_name_standard[
            gene_row
          ]
        ),
      
      gene_chrom =
        as.character(
          seqnames(genes)[
            gene_row
          ]
        ),
      
      gene_start =
        start(genes)[
          gene_row
        ],
      
      gene_end =
        end(genes)[
          gene_row
        ],
      
      gene_strand =
        as.character(
          strand(genes)[
            gene_row
          ]
        )
    )
  
  
  # ----------------------------------------------------------
  # Join back to original cluster information
  #
  # One cluster overlapping multiple genes will produce
  # multiple rows.
  # ----------------------------------------------------------
  
  cluster_with_index <- cluster_df %>%
    mutate(
      cluster_row =
        row_number()
    )
  
  
  annotated <- cluster_with_index %>%
    
    left_join(
      hit_df,
      by = "cluster_row"
    ) %>%
    
    select(
      -cluster_row,
      -gene_row
    )
  
  
  annotated
}


# ============================================================
# 10. Annotate HyperTRIBE
# ============================================================

hb_annotated <- annotate_clusters(
  hb,
  "HyperTRIBE"
)


write_tsv(
  hb_annotated,
  file.path(
    output_dir,
    "HB_clusters_gene_annotation.tsv"
  )
)


# ============================================================
# 11. Annotate STAMP
# ============================================================

stamp_annotated <- annotate_clusters(
  stamp,
  "STAMP"
)


write_tsv(
  stamp_annotated,
  file.path(
    output_dir,
    "STAMP_clusters_gene_annotation.tsv"
  )
)


# ============================================================
# 12. Unique gene lists
# ============================================================

hb_genes <- hb_annotated %>%
  
  filter(
    !is.na(gene_id)
  ) %>%
  
  distinct(
    gene_id,
    gene_name
  ) %>%
  
  arrange(
    gene_id
  )


stamp_genes <- stamp_annotated %>%
  
  filter(
    !is.na(gene_id)
  ) %>%
  
  distinct(
    gene_id,
    gene_name
  ) %>%
  
  arrange(
    gene_id
  )


write_tsv(
  hb_genes,
  file.path(
    output_dir,
    "HB_unique_genes.tsv"
  )
)


write_tsv(
  stamp_genes,
  file.path(
    output_dir,
    "STAMP_unique_genes.tsv"
  )
)


# ============================================================
# 13. Annotation statistics
# ============================================================

hb_annotated_regions <- hb_annotated %>%
  
  filter(
    !is.na(gene_id)
  ) %>%
  
  distinct(
    region_id
  ) %>%
  
  nrow()


stamp_annotated_regions <- stamp_annotated %>%
  
  filter(
    !is.na(gene_id)
  ) %>%
  
  distinct(
    region_id
  ) %>%
  
  nrow()


summary_table <- tibble(
  
  method = c(
    "HyperTRIBE",
    "STAMP"
  ),
  
  input_clusters = c(
    nrow(hb),
    nrow(stamp)
  ),
  
  clusters_overlapping_gene = c(
    hb_annotated_regions,
    stamp_annotated_regions
  ),
  
  proportion_clusters_annotated = c(
    hb_annotated_regions / nrow(hb),
    stamp_annotated_regions / nrow(stamp)
  ),
  
  unique_genes = c(
    nrow(hb_genes),
    nrow(stamp_genes)
  )
)


write_tsv(
  summary_table,
  file.path(
    output_dir,
    "HB_STAMP_gene_annotation_summary.tsv"
  )
)


cat("\n")
cat("========================================\n")
cat("Gene annotation summary\n")
cat("========================================\n\n")

print(summary_table)


# ============================================================
# 14. HyperTRIBE vs STAMP gene-level overlap
# ============================================================

hb_gene_ids <- unique(
  hb_genes$gene_id
)

stamp_gene_ids <- unique(
  stamp_genes$gene_id
)


shared_genes <- intersect(
  hb_gene_ids,
  stamp_gene_ids
)

hb_only_genes <- setdiff(
  hb_gene_ids,
  stamp_gene_ids
)

stamp_only_genes <- setdiff(
  stamp_gene_ids,
  hb_gene_ids
)


# ============================================================
# Summary
# ============================================================

overlap_summary <- tibble(
  
  category = c(
    "HyperTRIBE total",
    "STAMP total",
    "Shared",
    "HyperTRIBE only",
    "STAMP only"
  ),
  
  n_genes = c(
    length(hb_gene_ids),
    length(stamp_gene_ids),
    length(shared_genes),
    length(hb_only_genes),
    length(stamp_only_genes)
  )
)


cat("\n========================================\n")
cat("HyperTRIBE vs STAMP gene overlap\n")
cat("========================================\n\n")

print(overlap_summary)


cat(
  "\n% of HyperTRIBE genes recovered by STAMP:",
  round(
    100 *
      length(shared_genes) /
      length(hb_gene_ids),
    2
  ),
  "%\n"
)


cat(
  "% of STAMP genes recovered by HyperTRIBE:",
  round(
    100 *
      length(shared_genes) /
      length(stamp_gene_ids),
    2
  ),
  "%\n"
)


# ============================================================
# Save
# ============================================================

write_tsv(
  overlap_summary,
  file.path(
    output_dir,
    "HB_STAMP_gene_overlap_summary.tsv"
  )
)


write_tsv(
  tibble(
    gene_id = shared_genes
  ),
  file.path(
    output_dir,
    "HB_STAMP_shared_genes.tsv"
  )
)


write_tsv(
  tibble(
    gene_id = hb_only_genes
  ),
  file.path(
    output_dir,
    "HB_only_genes.tsv"
  )
)


write_tsv(
  tibble(
    gene_id = stamp_only_genes
  ),
  file.path(
    output_dir,
    "STAMP_only_genes.tsv"
  )
)