# This script will serialize a FacileTcgaDataSet to `fds.dir` from data that has
# been downloaded from xena and available locally in a parent directory
# `xena.dir`.
#
# We assume that the names of the local files are unchanged from their original
# downloaded names.
#
# The description when entired into the FacileDataSet will be taken from the
# `label` and `dataSubype` fields of its corresponding `*.json` file.
#
# In time, we will include increasingly more assays to the dataset, but for now
# we we have:
#
# 1. Estimated RSEM gene-level counts (tcga_gene_expected_count.gz)
# 2. Estimated Kallisto transcript-level counts (tcga_Kallisto_est_counts.gz)
# 3. real valued, gene-level copy number (Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz)
# 4. thresholed, gene-level copy number (Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz)
# 5. batch corrected miRNA (pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena.gz)
devtools::load_all(".")

library(FacileData)
library(dplyr)

xena.dir <- "~/workspace/data/FacileData/consortia/tcga/xena"
fds.dir <- "~/workspace/data/FacileData/consortia/tcga/FacileTcgaDataSet"

# 0. Load Description per indication/dataset
ind.description <- local({
  fn <- system.file("extdata", "tcga-indications.csv",
                    package = "FacileTcgaDataSet")
  info <- read.csv(fn, stringsAsFactors = FALSE)
  out <- lapply(split(info, info$code), function(x) {
    as.list(x[-1])
  })
})

# 1. Parse RSEM gene expression estimates into a a list of DGELists, one DGEList
#    per indication. The first round of sample-level covariates will be
#    extracted from the $samples data.frame from each DGEList.
gene.expr <- prep_main_gene_expression(xena.dir)
gc(verbose = TRUE, full = TRUE, reset = TRUE)

tcga.fds <- as.FacileDataSet(
  gene.expr,
  fds.dir,
  dataset_name = "FacileTcgaDataSet",
  dataset_meta = ind.description,
  assay_name = "rsem_gene",
  assay_description = "TOIL RSEM expected_count (PANCAN dataset)",
  assay_type = "rnaseq",
  organism = "Homo sapiens")
