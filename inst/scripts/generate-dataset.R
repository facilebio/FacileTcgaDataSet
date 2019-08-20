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
# 2. real valued, gene-level copy number (Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz)
# 3. thresholed, gene-level copy number (Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz)
# 4. batch corrected miRNA (pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena.gz)
# 5. Estimated Kallisto transcript-level counts (tcga_Kallisto_est_counts.gz)
devtools::load_all(".")

library(FacileData)
library(dplyr)

xena.dir <- "~/workspace/data/FacileData/consortia/tcga/xena"
fds.dir <- "~/workspace/data/FacileData/consortia/tcga/FacileTcgaDataSet"
stopifnot(
  dir.exists(xena.dir),
  dir.exists(dirname(fds.dir)))

# 0. Load Description per indication/dataset ===================================
# Generates a named (by indication) list of metadata per internal dataset
# (indication). This is not required for a FacileDataSet in general, but the
# more metadata we can store about the interanl data, the better / more handy
# it will be for downstream consumers.
ind.description <- local({
  fn <- system.file("extdata", "tcga-indications.csv",
                    package = "FacileTcgaDataSet")
  info <- read.csv(fn, stringsAsFactors = FALSE)
  out <- lapply(split(info, info$code), function(x) {
    as.list(x[-1])
  })
})

# 1. RSEM gene expected gene-level counts ======================================
# Parse RSEM gene expression estimates into a a list of DGELists, one DGEList
# per indication. The first round of sample-level covariates will be
# extracted from the $samples data.frame from each DGEList.
gene.expr <- prep_main_gene_expression(xena.dir)
gc(verbose = TRUE, full = TRUE, reset = TRUE)

# At max load, FacileDataSet creation eats up ~30GB of RAM.
#
# I think we can optimize this by serializing individual DGELists separately
# into the hdf5 file, then calculating the librarysize and normfactors in
# batches as is currently done in
# https://github.com/denalitherapeutics/archs4
tcga.fds <- as.FacileDataSet(
  gene.expr,
  fds.dir,
  dataset_name = "FacileTcgaDataSet",
  dataset_meta = ind.description,
  assay_name = "rsem_gene",
  assay_description = "TOIL RSEM expected_count (PANCAN dataset)",
  assay_type = "rnaseq",
  organism = "Homo sapiens")
# Time taken: 7.8 minutes

# 2. real valued gene-level copy number ========================================
# Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz
cnv.dat <- prep_gene_cnv(xena.dir, samples = samples(tcga.fds), discrete = FALSE)
cnv.f <- addFacileAssaySet(
  tcga.fds,
  cnv.dat$data,
  facile_feature_info = cnv.dat$features,
  facile_assay_name = "cnv_score",
  facile_assay_type = "real",
  facile_feature_type = "ensgid",
  facile_assay_description = "TCGA PANCAN gene-level copy number (gistic2)",
  storage_mode = "numeric")

# 3. discrete gene-level copy number ===========================================
# Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz
cnvd.dat <- prep_gene_cnv(xena.dir, samples = samples(tcga.fds), discrete = TRUE)
addFacileAssaySet(
  tcga.fds,
  cnvd.dat$data,
  facile_feature_info = cnvd.dat$features,
  facile_assay_name = "cnv_discrete",
  facile_assay_type = "discrete",
  facile_feature_type = "ensgid",
  facile_assay_description = paste(
    "TCGA PANCAN gene-level copy number (thresholded gistic2).",
    "-2: homozygous deletion; -1: single copy deletion; 0: diploid normal;",
    "1: low-level copy number amplification; 2: high-level copy number amplification"),
  storage_mode = "integer")

# 4. batch corrected miRNA (pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena.gz)
# 5. Estimated Kallisto transcript-level counts (tcga_Kallisto_est_counts.gz)
