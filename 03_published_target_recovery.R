

#!/usr/bin/env Rscript

# ============================================================
# Recovery of published AtGRP7 target genes
#
# This script compares HyperTRIBE and STAMP candidate target
# genes with previously published AtGRP7 iCLIP and RIP target
# datasets.
#
# Inputs:
#   1. HyperTRIBE candidate target genes
#   2. STAMP candidate target genes
#   3. Published AtGRP7 iCLIP target genes
#   4. Published AtGRP7 RIP target genes
#
# Output:
#   iCLIP_RIP_recovery_summary.tsv
# ============================================================


# ============================================================
# 1. Packages
# ============================================================

library(readr)
library(dplyr)


# ============================================================
# 2. Input files
# ============================================================

hb_candidate_file <-
  "/PATH/HB_unique_genes.tsv"

stamp_candidate_file <-
  "/PATH/STAMP_unique_genes.tsv"

iclip_file <-
  "/PATH/AtGRP7_iCLIP.txt"

rip_file <-
  "/PATH/RIP_AtGRP7.txt"

output_dir <-
  "/PATH/OUTPUT"


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. Read HyperTRIBE / STAMP candidate genes
# ============================================================

hb_df <- read_tsv(
  hb_candidate_file,
  show_col_types = FALSE
)

stamp_df <- read_tsv(
  stamp_candidate_file,
  show_col_types = FALSE
)


HB_all <- hb_df %>%
  filter(!is.na(gene_id)) %>%
  pull(gene_id) %>%
  toupper() %>%
  unique()


STAMP_all <- stamp_df %>%
  filter(!is.na(gene_id)) %>%
  pull(gene_id) %>%
  toupper() %>%
  unique()


# ============================================================
# 4. Read published iCLIP / RIP target genes
# ============================================================

iclip_genes <- read_lines(
  iclip_file
) %>%
  trimws() %>%
  toupper()

iclip_genes <- unique(
  iclip_genes[iclip_genes != ""]
)


rip_genes <- read_lines(
  rip_file
) %>%
  trimws() %>%
  toupper()

rip_genes <- unique(
  rip_genes[rip_genes != ""]
)


# ============================================================
# 5. Check input sizes
# ============================================================

cat("\nInput gene numbers:\n\n")

cat(
  "HyperTRIBE candidates:",
  length(HB_all),
  "\n"
)

cat(
  "STAMP candidates:",
  length(STAMP_all),
  "\n"
)

cat(
  "iCLIP reference genes:",
  length(iclip_genes),
  "\n"
)

cat(
  "RIP reference genes:",
  length(rip_genes),
  "\n\n"
)


# ============================================================
# 6. Calculate reference target recovery
# ============================================================

calculate_recovery <- function(
    candidate_genes,
    reference_genes,
    method_name,
    reference_name
) {
  
  overlap_ids <- intersect(
    candidate_genes,
    reference_genes
  )
  
  n_recovered <- length(overlap_ids)
  n_reference <- length(reference_genes)
  
  tibble(
    method = method_name,
    reference = reference_name,
    candidate_genes = length(candidate_genes),
    reference_genes = n_reference,
    recovered_genes = n_recovered,
    recovery_percent =
      100 * n_recovered / n_reference
  )
}


recovery_summary <- bind_rows(
  
  calculate_recovery(
    HB_all,
    iclip_genes,
    "HyperTRIBE",
    "iCLIP"
  ),
  
  calculate_recovery(
    STAMP_all,
    iclip_genes,
    "STAMP",
    "iCLIP"
  ),
  
  calculate_recovery(
    HB_all,
    rip_genes,
    "HyperTRIBE",
    "RIP"
  ),
  
  calculate_recovery(
    STAMP_all,
    rip_genes,
    "STAMP",
    "RIP"
  )
)


# ============================================================
# 7. Print recovery summary
# ============================================================

cat("\n========================================\n")
cat("Reference target recovery\n")
cat("========================================\n\n")

print(recovery_summary)


# ============================================================
# 8. Save recovery summary
# ============================================================

write_tsv(
  recovery_summary,
  file.path(
    output_dir,
    "iCLIP_RIP_recovery_summary.tsv"
  )
)


# ============================================================
# 9. Finish
# ============================================================

cat("\n")
cat("Recovery analysis complete.\n")

cat(
  "Output:",
  file.path(
    output_dir,
    "iCLIP_RIP_recovery_summary.tsv"
  ),
  "\n"
)
