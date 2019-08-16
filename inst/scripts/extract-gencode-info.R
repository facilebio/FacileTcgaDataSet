# As of 2019-08-15, the expression data in xena was processed against the
# GENCODE v23 annotations.
#
# The metadata provided by xena is *almost* all that we need
# (ie. `gencode.v23.annotation.transcript.probemap`), however I want to put the
# gene biotype information in the `meta` columne of the feature_info
# table.
#
# This script works over the gencode.v23.annotation.gtf.gz to extract this
# information. I wanted to initially save these data in the `inst/extdata` of
# this package, but the gene and tx info when gzipped is still 1.4M and 4.6M,
# respectively, so ... let's not.
#
# Put the `gencode.v23.annotation.gtf.gz` file in your local
# `XENA_DATA_DIR` directory, and this script will parse that and serialize
# gene-info.rds and transcript-info.rds files.
#
# Why not *.csv.gz files? I want to store the order of the biotype levels,
# and this is easy and space efficient enough.
xena.data.dir <- "~/workspace/data/FacileData/consortia/tcga/xena"
gtf.path <- file.path(xena.data.dir, "gencode.v23.annotation.gtf.gz")

stopifnot(
  dir.exists(xena.data.dir),
  file.exists(gtf.path))

out.fn <- c(
  "gene" = file.path(xena.data.dir, "gene-info.rds"),
  "tx"   = file.path(xena.data.dir, "transcript-info.rds"))

library(FacileData)
library(dplyr)

info <- FacileData::extract_transcribed_info_from_ensembl_gtf(gtf.path)
saveRDS(info$gene_info, out.fn['gene'])
saveRDS(info$transcript_info, out.fn['tx'])

# Now that that's done, let's see how it compares to the data that xena provides
# in their "gencode.v23.annotation.gene.probemap" and
# "gencode.v23.annotation.transcript.probemap" files
library(readr)
xena.dir <- dirname(gtf.path)
xena.ginfo <- file.path(xena.data.dir, "gencode.v23.annotation.gene.probemap") %>%
  read_tsv()
xena.txinfo <- file.path(xena.dir, "gencode.v23.annotation.transcript.probemap") %>%
  read_tsv()

# All gene information matches up bueno
cmp.gene <- full_join(info$gene_info, xena.ginfo, by = c("gene_id" = "id"))
setequal(cmp.gene$gene_id, info$gene_info$gene_id)
all.equal(cmp.gene$symbol, cmp.gene$gene)
all.equal(cmp.gene$start, cmp.gene$chromStart)

# All transcript information matches up bueno
cmp.tx <- full_join(info$transcript_info, xena.txinfo,
                    by = c("transcript_id" = "id"))
setequal(cmp.tx$transcript_id, info$transcript_info$transcript_id)
all.equal(cmp.tx$gene_name, cmp.tx$gene)
all.equal(cmp.tx$start, cmp.tx$chromStart)
