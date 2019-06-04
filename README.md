
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Overview

This package leverages the MultiAssayExperiment and curatedTCGAData
packages to download the relevant data, and assemble a singular,
pan-cancer, FacileTcgaDataSet using a subset of the assays made
availalbe from those resources.

Gene-level quantitaion (expression, CNV, etc) will be mapped to the
ensembl gene-level universe. Entries that cannot be mapped that way will
be dropped.

## Datasets

We will create a dataset across all indications using the following
assays:

### mRNA Abundance

The following assays can be pulled from the curatedTCGAData package:

  - `RNASeq2GeneNorm`: fpkm values (sad)
  - `RNASeqGene`: v2 fpkm values?
  - `miRNASeqGene`: miRNA-seq data

We might consider getting gene and transcript level counts from the
recount2 project, though:

<https://jhubiostatistics.shinyapps.io/recount/>

They used Gencode v25 GFF3 annotations, which we can parse with the
utility funcionts in GemomicsStudyDb

### Copy Number

  - `GISTIC_AllByGene`: real valued copy number per gene
  - `GISTIC_ThresholdedByGene`: threshold (duplication / deletion /
    normal) gene-level CNV scores

### GenomicVariants

  - `Mutation`: mutation status? What does this look like?

### Other

  - `RPPAArray`: Reverse Phase Protein Array (NOISY)
  - `Methylation`, `Methylation_methyl27`, or `Methylation_methyl450`

Individual indication will be downloaded separately and assembled into
one final FacileTcgaDataSet.

## Data Assembly

Let’s follow along with the [curatedTCGAData
vignette](http://bioconductor.org/packages/release/data/experiment/vignettes/curatedTCGAData/inst/doc/curatedTCGAData.html)
to figure out how to do this.
