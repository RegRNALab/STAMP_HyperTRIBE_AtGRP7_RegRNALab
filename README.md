# AtGRP7 HyperTRIBE and STAMP Analysis Pipeline

R scripts for the identification and functional annotation of transcriptome-wide RNA targets of AtGRP7 in *Arabidopsis thaliana* roots and cotyledons using HyperTRIBE and STAMP.

Associated publication:

> Yan X & Fernandez-Valverde SL (2025). Expanding the RNA-binding landscape of AtGRP7 in *Arabidopsis thaliana* seedlings using rapid transient HyperTRIBE and STAMP. *[Journal]*.

---

## Overview

This pipeline takes per-replicate FLARE Mode 2 editing-fraction outputs as input and proceeds through the following steps:

1. High-confidence editing cluster identification
2. Cluster-to-gene annotation
3. Recovery of previously published AtGRP7 targets
4. Novel candidate target identification
5. Gene Ontology enrichment analysis
6. Transcription factor and lncRNA annotation
7. Statistical assessment of HyperTRIBE–STAMP target overlap

All analyses are performed independently for roots and cotyledons unless otherwise noted.

---

## Requirements

### R version
R >= 4.1.0

### R packages
| Package | Source | Use |
|---|---|---|
| `readr` | CRAN | File I/O |
| `dplyr` | CRAN | Data manipulation |
| `tidyr` | CRAN | Data reshaping |
| `purrr` | CRAN | Functional programming |
| `data.table` | CRAN | Large file handling |
| `edgeR` | Bioconductor | CPM normalisation |
| `GenomicRanges` | Bioconductor | Genomic interval operations |
| `rtracklayer` | Bioconductor | GTF/GFF3 import |
| `clusterProfiler` | Bioconductor | GO enrichment analysis |
| `org.At.tair.db` | Bioconductor | Arabidopsis gene annotations |

Install CRAN packages with `install.packages()` and Bioconductor packages with `BiocManager::install()`.

---

## Input Data

| File | Description |
|---|---|
| FLARE Mode 2 `.tsv` files | Per-replicate editing fraction outputs from [FLARE](https://github.com/YeoLab/FLARE); 3 treatment and 3 control replicates per tissue |
| `TAIR12_cleaned.gtf` | TAIR12 genome annotation; must match the assembly used for read alignment and FLARE analysis |
| `gene_counts_24samples.txt` | featureCounts gene-level count matrix across all 24 libraries |
| `iCLIP_targets.txt` | Previously published AtGRP7 iCLIP target gene list (Meyer et al., 2017) |
| `RIPseq_targets.txt` | Previously published AtGRP7 RIP-seq target gene list (Meyer et al., 2017) |
| PlantTFDB TF annotation | Arabidopsis TF list from [PlantTFDB v5.0](https://planttfdb.gao-lab.org/) |
| TAIR12 lncRNA annotation | Non-coding RNA annotation from [Arabase](https://arabidopsis.ac.cn/) |

---

## Scripts

### `01_high_confidence_cluster_filtering.R`

Identifies high-confidence editing clusters from FLARE Mode 2 outputs.

**High-confidence criteria (both must be satisfied):**
- Editing detected in ≥ 2 of 3 treatment replicates
- Editing fraction in ≥ 2 treatment replicates exceeds the maximum editing fraction observed across all 3 control replicates

**User settings:** Set `tissue` to `"root"` or `"cotyledon"` and `method` to `"HyperTRIBE"` or `"STAMP"` at the top of the script.

**Outputs** (written to `results/<tissue>/HyperTRIBE_high_confidence_clusters/`):
- `<tissue>_HyperTRIBE_all_regions_with_consistency.tsv` — all candidate regions with replicate-level consistency metrics
- `<tissue>_HyperTRIBE_high_confidence_clusters.tsv` — filtered high-confidence clusters
- `<tissue>_HyperTRIBE_high_confidence_summary.tsv` — filtering step summary statistics
- `<tissue>_HyperTRIBE_reproducibility_summary.tsv` — breakdown by replicate reproducibility class

---

### `02_cluster_to_target_gene_annotation.R`

Assigns high-confidence HyperTRIBE and STAMP editing clusters to genes using strand-specific genomic overlap with the TAIR12 annotation. Identifies shared and method-specific candidate target gene sets.

**User settings:** Set input paths for HyperTRIBE clusters, STAMP clusters, genome annotation file, and output directory.

**Outputs:**
- `HB_clusters_gene_annotation.tsv` — HyperTRIBE clusters with gene assignments
- `STAMP_clusters_gene_annotation.tsv` — STAMP clusters with gene assignments
- `HB_unique_genes.tsv` — unique HyperTRIBE target genes
- `STAMP_unique_genes.tsv` — unique STAMP target genes
- `HB_STAMP_shared_genes.tsv` — genes identified by both methods
- `HB_only_genes.tsv` — HyperTRIBE-specific target genes
- `STAMP_only_genes.tsv` — STAMP-specific target genes
- `HB_STAMP_gene_annotation_summary.tsv` — cluster annotation statistics
- `HB_STAMP_gene_overlap_summary.tsv` — gene-level overlap summary

---

### `03_published_target_recovery.R`

Calculates the proportion of previously published AtGRP7 iCLIP and RIP-seq target genes (Meyer et al., 2017) recovered by HyperTRIBE and STAMP in each tissue.

**User settings:** Set `tissue` to `"root"` or `"cotyledon"`.

**Outputs:** Recovery statistics for iCLIP and RIP-seq targets by method and tissue.

---

### `04_novel_target_identification.R`

Identifies novel candidate AtGRP7 targets by excluding genes previously reported by iCLIP or RIP-seq. Produces shared and method-specific novel candidate gene sets for downstream analyses.

**Outputs:** Novel candidate gene lists stratified by method (shared, HyperTRIBE-specific, STAMP-specific) and tissue.

---

### `05a_shared_novel_GO_enrichment.R`

Performs Gene Ontology (GO) Biological Process enrichment analysis on shared novel AtGRP7 candidate targets identified by both HyperTRIBE and STAMP. Tissue-specific expressed genes (CPM ≥ 1 in ≥ 3 libraries) are used as the background gene universe. P values are adjusted using the Benjamini–Hochberg method; terms with adjusted P < 0.05 are retained. Redundant GO terms are reduced using semantic similarity (cutoff = 0.7).

**Outputs:** Enriched GO terms for shared novel targets in roots and cotyledons.

---

### `05b_method_specific_novel_GO_enrichment.R`

Performs the same GO enrichment analysis as `05a` for HyperTRIBE-specific and STAMP-specific novel candidate target genes.

**Outputs:** Enriched GO terms for method-specific novel targets in roots and cotyledons.

---

### `06a_TF_annotation.R`

Identifies transcription factors (TFs) among shared novel AtGRP7 candidate targets using the PlantTFDB v5.0 Arabidopsis TF annotation. Reports TF family distributions across roots and cotyledons.

**Outputs:** TF candidates with family annotations in roots and cotyledons.

---

### `06b_lncRNA_annotation.R`

Identifies long non-coding RNAs (lncRNAs) among shared novel AtGRP7 candidate targets using the TAIR12 non-coding RNA annotation from Arabase.

**Outputs:** lncRNA candidates in roots and cotyledons.

---

### `07_HyperTRIBE_STAMP_overlap_hypergeometric_test.R`

Tests whether the observed gene-level overlap between HyperTRIBE and STAMP candidate target sets is greater than expected by chance using a right-tailed hypergeometric test. Tissue-specific expressed genes (CPM ≥ 1 in ≥ 3 libraries) serve as the background. Analyses are performed separately for roots and cotyledons.

**Outputs:** Expressed-gene background lists and hypergeometric test results per tissue.

---

## Recommended Execution Order

```
01_high_confidence_cluster_filtering.R     # run for each method × tissue combination
02_cluster_to_target_gene_annotation.R    # run per tissue
03_published_target_recovery.R             # run per tissue
04_novel_target_identification.R           # run once
05a_shared_novel_GO_enrichment.R          # run once
05b_method_specific_novel_GO_enrichment.R # run once
06a_TF_annotation.R                        # run once
06b_lncRNA_annotation.R                    # run once
07_HyperTRIBE_STAMP_overlap_hypergeometric_test.R  # run once
```

---

## Notes

- Script `01` must be run for all four method × tissue combinations (HyperTRIBE root, HyperTRIBE cotyledon, STAMP root, STAMP cotyledon) before running script `02`.
- All file paths marked `/PATH/` must be updated to reflect your local directory structure before running.
- Genome annotation and read alignment must use the same assembly version (TAIR12 used here).
- Computational analyses were performed on the Katana HPC cluster at UNSW Sydney (https://doi.org/10.26190/669x-a286).

---

## Contact

Selene L. Fernandez-Valverde
School of Biotechnology and Biomolecular Sciences / UNSW RNA Institute
The University of New South Wales, Sydney, Australia
