#' Prepare gene expression as main datasource for dataset construction.
#'
#' The estimated RSEM gene-level counts from `tcga_gene_expected_count.gz` will
#' be used as the main assay for the dataset. This means a number of things, but
#' one thing it means is that if we don't make an attempt to add the basel-level
#' covariates for samples not found in the `tcga_gene_expected_count.gz`, they
#' won't appear anywhere in the final FacileTcgaDataSet.
#'
#' This function will parse this count file into a list of DGEList objects
#' (one per indication). The DGELists will have basic annotation tied to them,
#' as generated from the [parse_base_sample_covariates()] function.
#'
#' To achieve this, we wrangle / harmonize three sources of data:
#'
#' 1. The counts data (`tcga_gene_expected_count.gz`)
#' 2. The gene-level data (generated in `scripts/extract-gencode-info.R`)
#' 3. A base set of sample-level covariate data, generated in
#'    [parse_base_sample_covariates()].
#'
#' @export
#' @importFrom jsonlite fromJSON
#' @importFrom readr read_tsv
#' @importFrom edgeR DGEList
#'
#' @param xena.dir The path to the directory where you've downloaded the xena
#'   data files. This function will load the `tcga_gene_expected_count.gz` and
#'   `tcga_gene_expected_count.json` files
#' @return A list of DGELists with some $samples info, names are TCGA
#'   indications
prep_main_gene_expression <- function(xena.dir, ..., debug = FALSE) {
  if (FALSE) {
    xena.dir <- "~/workspace/data/FacileData/consortia/tcga/xena"
  }
  assert_directory(xena.dir, "r")
  dat.fn <- file.path(xena.dir, "tcga_gene_expected_count.gz") %>%
    assert_file_exists()
  meta.fn <- file.path(xena.dir, "tcga_gene_expected_count.json") %>%
    assert_file_exists()

  # counts .....................................................................
  # counts <- read_tsv(dat.fn)
  counts <- data.table::fread(dat.fn)
  colnames(counts)[1L] <- "feature_id"
  counts[["feature_id"]] <- strip_ensembl_version(counts[["feature_id"]])

  # gene-level meta data .......................................................
  ginfo <- load_gene_info(xena.dir, strict_cast = TRUE)
  ginfo[["meta"]] <- as.character(ginfo[["meta"]])
  if (!setequal(ginfo$feature_id, counts$feature_id)) {
    warning("Mismatch of features in annotation and those in gene expression")
    if (debug) {
      browser()
    } else {
      stop("You probably downloaded the wrong GENCODE annotation file, ",
           "v23 is required")
    }
  }
  ginfo <- ginfo[match(counts$feature_id, ginfo$feature_id),]
  stopifnot(isTRUE(all.equal(ginfo$feature_id, counts$feature_id)))

  # sample-level meta data .....................................................
  scovs <- parse_base_sample_covariates(xena.dir)

  keep.ids <- intersect(colnames(counts), scovs[["sample_id"]])
  no.sample.meta <- setdiff(colnames(counts)[-1L], keep.ids)
  if (length(no.sample.meta)) {
    warning("Missing sample covariate information for ", length(no.sample.meta),
            " expression samples. These expression samples will be dropped")
  }
  no.expr.data <- setdiff(scovs[["sample_id"]], keep.ids)
  if (length(no.expr.data)) {
    warning("There are ", length(no.expr.data), " samples that do not have ",
            "expression data")
  }

  scovs <- filter(scovs, sample_id %in% keep.ids)
  E <- as.matrix(counts[, keep.ids, with = FALSE])
  rownames(E) <- counts$feature_id
  rm(counts)
  gc()

  inds <- sort(unique(scovs$indication))
  dats <- sapply(inds, function(ind) {
    message("==== Constructing ", ind, " DGEList ==============================")
    sinfo <- filter(scovs, indication == ind) %>% as.data.frame()
    rownames(sinfo) <- sinfo[["sample_id"]]

    # These data are `log2(expected_count+1)`, so let's transform these back to
    # to the count scale, and round to integer since this saves mucho space
    # and doesn't lose any real appreciable information.
    cnts <- round(2**E[, sinfo$sample_id] - 1)
    storage.mode(cnts) <- "integer"
    edgeR::DGEList(cnts, samples = sinfo, genes = ginfo)
  }, simplify = FALSE)
  dats
}

#' Load the gene-level information
#'
#' You need to have run the `inst/scripts/extract-gencode-info.R` script over
#' the gencode v23 files and save the outputs to `xena.dir`.
#'
#' The gene level information will be returned in the order expected for the
#' FacileDataSet
#'
#' @param xena.dir the top-level data directory for the xena data
#' @param strict_cast When TRUE (default), the columns are cast to the type
#'   expected in the FacileDataSet constructors. Currently we have strict checks
#'   for character, numeric, etc. for the given column names.
load_gene_info <- function(xena.dir, strict_cast = TRUE, ...) {
  assert_directory(xena.dir, "r")
  gene.fn <- assert_file_exists(file.path(xena.dir, "gene-info.rds"))
  gene.info <- readRDS(gene.fn)

  finfo <- gene.info %>%
    transmute(
      feature_type = "ensgid",
      feature_id = strip_ensembl_version(gene_id),
      name = symbol,
      meta = gene_type,
      seqnames, start, end, strand = ifelse(strand == "+", 1L, -1L),
      effective_length = as.integer(length),
      source = "GENCODEv23")
  rownames(finfo) <- finfo$feature_id

  if (strict_cast) {
    finfo[["meta"]] <- as.character(finfo[["meta"]])
    finfo[["seqnames"]] <- as.character(finfo[["seqnames"]])
  }

  finfo
}

#' @noRd
strip_ensembl_version <- function(x, ...) {
  sub("\\.\\d+$", "", x)
}
