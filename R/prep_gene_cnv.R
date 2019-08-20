#' @noRd
#' @export
#' @param xena.dir path to local xena data directory
#' @param samples a dataset,sample_id data.frame of samples to restrict
#'   data extraction to. When `NULL` (default), all samples are used
#' @param discrete when `FALSE` (default), the real valued scores are procssed,
#'   otherwise we parse the tresholded scores
prep_gene_cnv <- function(xena.dir, samples = NULL, discrete = FALSE, ...) {
  assert_directory(xena.dir, "r")
  if (is.null(samples)) {
    samples <- parse_base_sample_covariates(xena.dir)
  } else {
    samples <- collect(samples, n = Inf)
  }
  assert_multi_class(samples, c("data.frame", "tbl"))
  assert_subset(c("dataset", "sample_id"), colnames(samples))

  assert_flag(discrete)
  fn <- paste0(
    "Gistic2_CopyNumber_Gistic2_all_",
    if (discrete) "thresholded.by_genes.gz" else "data_by_genes.gz")
  fn <- assert_file_exists(file.path(xena.dir, fn))

  message("=== Loading ", fn, " ====================================================")
  xdat <- data.table::fread(fn)
  colnames(xdat)[1] <- "symbol"

  features <- match_gene_symbol_to_feature(xena.dir, xdat$symbol)
  keep <- !is.na(features$feature_id)

  xdat <- xdat[keep,]
  features <- features[keep,]
  stopifnot(isTRUE(all.equal(xdat$symbol, features$name)))

  cnv.dat <- as.matrix(xdat[, -1, with = FALSE])
  rownames(cnv.dat) <- features$feature_id
  if (discrete) {
    storage.mode(cnv.dat) <- "integer"
  }
  keep.dat <- intersect(colnames(cnv.dat), samples[["sample_id"]])
  # no.sample.meta <- setdiff(colnames(cnv.dat), keep.dat)
  # if (length(no.sample.meta)) {
  #   warning("Missing sample covariate information for ", length(no.sample.meta),
  #           " cnv samples. These data will be dropped")
  # }
  # no.cnv.data <- setdiff(samples[["sample_id"]], keep.dat)
  # if (length(no.cnv.data)) {
  #   warning("There are ", length(no.cnv.data), " samples that do not have ",
  #           "cnv data")
  # }

  message("Inserting CNV data for ", length(keep.dat), " / ", ncol(cnv.dat), " samples")
  out <- sapply(unique(samples[["dataset"]]), function(ind) {
    message("==== Extracting ", ind, " CNV data  ==============================")
    sinfo <- filter(samples, dataset == ind)
    cnv.dat[, intersect(sinfo[["sample_id"]], colnames(cnv.dat))]
  }, simplify = FALSE)

  list(features = features, data = out)
}

match_gene_symbol_to_feature <- function(xena.dir, symbol, ...) {
  ginfo <- load_gene_info(xena.dir, strict_cast = TRUE)
  xref <- match(symbol, ginfo[["name"]])
  ginfo[xref,]
}
