
#!/usr/bin/env Rscript

# ============================================================
# HyperTRIBE or STAMP high-confidence cluster filtering
#
# This script processes FLARE Mode 2 outputs and identifies
# high-confidence HyperTRIBE or STAMP editing clusters in root or
# cotyledon samples.
#
# High-confidence definition:
#
#   1. HyperTRIBE or STAMP editing is detected in >= 2/3 biological
#      replicates
#
#   AND
#
#   2. The editing fraction in >= 2/3 HyperTRIBE or STAMP replicates
#      is greater than the maximum editing fraction observed
#      across all three CTRL replicates.
#
# Pooled editing fraction:
#
#   sum(edited bases across three replicates) /
#   sum(target bases across three replicates)
#
# Delta pooled editing fraction:
#
#   HyperTRIBE or STAMP pooled EF - CTRL pooled EF
#
# The same analysis is applied independently to root and
# cotyledon samples.
# ============================================================


# ============================================================
# 1. Load packages
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(tidyr)
})


# ============================================================
# 2. User settings
# ============================================================

# Select tissue:
#   "root"
#   "cotyledon"

tissue <- "cotyledon"


# ============================================================
# 3. Input files
# ============================================================

if (tissue == "root") {
  
  files <- c(
    HB1   = "data/root/HB1_edit_fraction_in_regions.tsv",
    HB2   = "data/root/HB2_edit_fraction_in_regions.tsv",
    HB3   = "data/root/HB3_edit_fraction_in_regions.tsv",
    CTRL1 = "data/root/CTRL1_edit_fraction_in_regions.tsv",
    CTRL2 = "data/root/CTRL2_edit_fraction_in_regions.tsv",
    CTRL3 = "data/root/CTRL3_edit_fraction_in_regions.tsv"
  )
  
} else if (tissue == "cotyledon") {
  
  files <- c(
    HB1   = "data/cotyledon/HB1_edit_fraction_in_regions.tsv",
    HB2   = "data/cotyledon/HB2_edit_fraction_in_regions.tsv",
    HB3   = "data/cotyledon/HB3_edit_fraction_in_regions.tsv",
    CTRL1 = "data/cotyledon/CTRL1_edit_fraction_in_regions.tsv",
    CTRL2 = "data/cotyledon/CTRL2_edit_fraction_in_regions.tsv",
    CTRL3 = "data/cotyledon/CTRL3_edit_fraction_in_regions.tsv"
  )
  
} else {
  
  stop("tissue must be either 'root' or 'cotyledon'.")
  
}


output_dir <- file.path(
  "results",
  tissue,
  "HyperTRIBE_high_confidence_clusters"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


expected_samples <- c(
  "HB1", "HB2", "HB3",
  "CTRL1", "CTRL2", "CTRL3"
)

if (!identical(names(files), expected_samples)) {
  stop(
    "Input files must be named exactly: ",
    paste(expected_samples, collapse = ", ")
  )
}


# ============================================================
# 4. Read one FLARE Mode 2 file
# ============================================================

read_mode2 <- function(file, sample_name) {
  
  if (!file.exists(file)) {
    stop("File does not exist: ", file)
  }
  
  message("Reading ", sample_name, ": ", file)
  
  df <- read_tsv(
    file,
    show_col_types = FALSE,
    progress = FALSE
  )
  
  
  required_columns <- c(
    "chrom",
    "start",
    "end",
    "edit_fraction",
    "strand",
    "target_bases",
    "edited_bases"
  )
  
  
  missing_columns <- setdiff(
    required_columns,
    colnames(df)
  )
  
  
  if (length(missing_columns) > 0) {
    
    stop(
      "Missing required columns in ",
      file,
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  
  # ----------------------------------------------------------
  # Identify region ID column
  # ----------------------------------------------------------
  
  if ("region_id" %in% colnames(df)) {
    
    # Already correctly named.
    
  } else if ("region" %in% colnames(df)) {
    
    df <- df %>%
      rename(region_id = region)
    
  } else {
    
    extra_columns <- setdiff(
      colnames(df),
      required_columns
    )
    
    if (length(extra_columns) == 1) {
      
      df <- df %>%
        rename(
          region_id = all_of(extra_columns)
        )
      
    } else {
      
      stop(
        "Could not identify region ID column in ",
        file
      )
    }
  }
  
  
  df %>%
    transmute(
      
      region_id = as.character(region_id),
      
      chrom = as.character(chrom),
      
      start = as.integer(start),
      
      end = as.integer(end),
      
      strand = as.character(strand),
      
      sample = sample_name,
      
      condition = if_else(
        grepl("^HB", sample_name),
        "HB",
        "CTRL"
      ),
      
      replicate = case_when(
        grepl("1$", sample_name) ~ 1L,
        grepl("2$", sample_name) ~ 2L,
        grepl("3$", sample_name) ~ 3L,
        TRUE ~ NA_integer_
      ),
      
      edit_fraction = as.numeric(edit_fraction),
      
      target_bases = as.numeric(target_bases),
      
      edited_bases = as.numeric(edited_bases)
    )
}


# ============================================================
# 5. Read and combine all six samples
# ============================================================

all_data <- imap_dfr(
  files,
  ~ read_mode2(
    file = .x,
    sample_name = .y
  )
)


message(
  "Combined dimensions: ",
  nrow(all_data),
  " rows x ",
  ncol(all_data),
  " columns"
)


# ============================================================
# 6. Basic data checks
# ============================================================

if (anyNA(all_data$region_id)) {
  stop("Missing region IDs detected.")
}


if (anyNA(all_data$target_bases)) {
  stop("NA values detected in target_bases.")
}


if (anyNA(all_data$edited_bases)) {
  stop("NA values detected in edited_bases.")
}


if (any(all_data$target_bases < 0, na.rm = TRUE)) {
  stop("Negative target_bases values detected.")
}


if (any(all_data$edited_bases < 0, na.rm = TRUE)) {
  stop("Negative edited_bases values detected.")
}


duplicate_check <- all_data %>%
  count(
    sample,
    region_id,
    name = "n"
  ) %>%
  filter(n != 1)


if (nrow(duplicate_check) > 0) {
  
  stop(
    "Some sample-region combinations do not occur exactly once."
  )
}


region_sample_check <- all_data %>%
  count(
    region_id,
    name = "n_samples"
  ) %>%
  filter(n_samples != 6)


if (nrow(region_sample_check) > 0) {
  
  stop(
    nrow(region_sample_check),
    " regions are not represented in all six samples."
  )
}


# ============================================================
# 7. Convert replicate-level data to wide format
# ============================================================

ef_wide <- all_data %>%
  
  select(
    region_id,
    chrom,
    start,
    end,
    strand,
    sample,
    edit_fraction,
    target_bases,
    edited_bases
  ) %>%
  
  pivot_wider(
    
    names_from = sample,
    
    values_from = c(
      edit_fraction,
      target_bases,
      edited_bases
    ),
    
    names_glue = "{sample}_{.value}"
  )


# ============================================================
# 8. Calculate replicate consistency and editing fractions
# ============================================================

region_summary <- ef_wide %>%
  
  mutate(
    
    # --------------------------------------------------------
    # Editing detection in individual replicates
    # --------------------------------------------------------
    
    HB1_detected =
      HB1_edited_bases > 0,
    
    HB2_detected =
      HB2_edited_bases > 0,
    
    HB3_detected =
      HB3_edited_bases > 0,
    
    
    CTRL1_detected =
      CTRL1_edited_bases > 0,
    
    CTRL2_detected =
      CTRL2_edited_bases > 0,
    
    CTRL3_detected =
      CTRL3_edited_bases > 0,
    
    
    # --------------------------------------------------------
    # Number of replicates with detected editing
    # --------------------------------------------------------
    
    HB_reps =
      as.integer(HB1_detected) +
      as.integer(HB2_detected) +
      as.integer(HB3_detected),
    
    
    CTRL_reps =
      as.integer(CTRL1_detected) +
      as.integer(CTRL2_detected) +
      as.integer(CTRL3_detected),
    
    
    # --------------------------------------------------------
    # Maximum CTRL editing fraction
    # --------------------------------------------------------
    
    CTRL_max_EF = pmax(
      CTRL1_edit_fraction,
      CTRL2_edit_fraction,
      CTRL3_edit_fraction,
      na.rm = TRUE
    ),
    
    
    # --------------------------------------------------------
    # Compare each HB replicate with maximum CTRL EF
    # --------------------------------------------------------
    
    HB1_above_CTRLmax =
      HB1_edit_fraction > CTRL_max_EF,
    
    HB2_above_CTRLmax =
      HB2_edit_fraction > CTRL_max_EF,
    
    HB3_above_CTRLmax =
      HB3_edit_fraction > CTRL_max_EF,
    
    
    n_HB_above_CTRLmax =
      as.integer(HB1_above_CTRLmax) +
      as.integer(HB2_above_CTRLmax) +
      as.integer(HB3_above_CTRLmax),
    
    
    # --------------------------------------------------------
    # Pooled editing fractions
    # --------------------------------------------------------
    
    HB_total_target =
      HB1_target_bases +
      HB2_target_bases +
      HB3_target_bases,
    
    
    HB_total_edited =
      HB1_edited_bases +
      HB2_edited_bases +
      HB3_edited_bases,
    
    
    CTRL_total_target =
      CTRL1_target_bases +
      CTRL2_target_bases +
      CTRL3_target_bases,
    
    
    CTRL_total_edited =
      CTRL1_edited_bases +
      CTRL2_edited_bases +
      CTRL3_edited_bases,
    
    
    HB_pooled_EF = if_else(
      HB_total_target > 0,
      HB_total_edited / HB_total_target,
      NA_real_
    ),
    
    
    CTRL_pooled_EF = if_else(
      CTRL_total_target > 0,
      CTRL_total_edited / CTRL_total_target,
      NA_real_
    ),
    
    
    delta_pooled_EF =
      HB_pooled_EF -
      CTRL_pooled_EF,
    
    
    # --------------------------------------------------------
    # Reproducibility class
    # --------------------------------------------------------
    
    reproducibility_class =
      paste0(
        "HB",
        HB_reps,
        "_CTRL",
        CTRL_reps
      ),
    
    
    # --------------------------------------------------------
    # High-confidence criteria
    # --------------------------------------------------------
    
    pass_HB_reproducibility =
      HB_reps >= 2,
    
    
    pass_strong_consistency =
      n_HB_above_CTRLmax >= 2,
    
    
    high_confidence =
      pass_HB_reproducibility &
      pass_strong_consistency
  )


# ============================================================
# 9. Save complete region-level table
# ============================================================

write_tsv(
  region_summary,
  file.path(
    output_dir,
    paste0(
      tissue,
      "_HyperTRIBE_all_regions_with_consistency.tsv"
    )
  )
)


# ============================================================
# 10. Extract high-confidence clusters
# ============================================================

high_confidence_clusters <- region_summary %>%
  
  filter(
    high_confidence
  )


write_tsv(
  high_confidence_clusters,
  file.path(
    output_dir,
    paste0(
      tissue,
      "_HyperTRIBE_high_confidence_clusters.tsv"
    )
  )
)


# ============================================================
# 11. Summary statistics
# ============================================================

summary_table <- tibble(
  
  metric = c(
    
    "Total treatment-specific regions",
    
    "HB editing detected >=2/3",
    
    ">=2 HB replicates > max CTRL EF",
    
    "Final high-confidence clusters"
  ),
  
  
  n_regions = c(
    
    nrow(region_summary),
    
    sum(
      region_summary$pass_HB_reproducibility,
      na.rm = TRUE
    ),
    
    sum(
      region_summary$pass_strong_consistency,
      na.rm = TRUE
    ),
    
    sum(
      region_summary$high_confidence,
      na.rm = TRUE
    )
  )
) %>%
  
  mutate(
    
    proportion_of_total =
      n_regions /
      first(n_regions)
  )


write_tsv(
  summary_table,
  file.path(
    output_dir,
    paste0(
      tissue,
      "_HyperTRIBE_high_confidence_summary.tsv"
    )
  )
)


# ============================================================
# 12. Save reproducibility-class summary
# ============================================================

reproducibility_summary <- region_summary %>%
  
  count(
    HB_reps,
    CTRL_reps,
    n_HB_above_CTRLmax,
    high_confidence,
    name = "n_regions"
  ) %>%
  
  arrange(
    desc(HB_reps),
    CTRL_reps,
    desc(n_HB_above_CTRLmax)
  )


write_tsv(
  reproducibility_summary,
  file.path(
    output_dir,
    paste0(
      tissue,
      "_HyperTRIBE_reproducibility_summary.tsv"
    )
  )
)


# ============================================================
# 13. Finish
# ============================================================

cat("\n")
cat("============================================\n")
cat("HyperTRIBE high-confidence analysis complete\n")
cat("============================================\n\n")

cat("Tissue:", tissue, "\n")

cat(
  "Total treatment-specific regions:",
  nrow(region_summary),
  "\n"
)

cat(
  "High-confidence clusters:",
  nrow(high_confidence_clusters),
  "\n\n"
)

cat(
  "Output directory:",
  output_dir,
  "\n"
)