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
library(dplyR)

fds.dir <- "~/workspace/data/FacileData/consortia/tcga/FacileTcgaDataSet"

