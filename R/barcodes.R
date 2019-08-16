#' Extracts the sample type (tumor/normal) from the barcode
#'
#' The sample type is encoded in the TCGA barcodes at the `[]` location,
#' `TCGA-02-0001-[01]C-01D-0183-01`. Codes between 01-09 are tumors, 10-19
#' are normal.
#'
#' More detailed code information can be found in the following file:
#' `inst/extdata/tcga-barcode-codebook.xlsx,sample_type`
#'
#' @export
#' @param x a vector of TCGA barcodes
barcode_sample_type <- function(x, ...) {
  x <- assert_tcga_barcode(x)

}

#' Converts a barcode into a data.frame of semantic pieces.
#'
#' Each bit of the TCGA barcode means something, as outlined here:
#' https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
#'
#' A full barcode follows this format: `TCGA-02-0001-[01]C-01D-0183-01`, however
#' the data downloaded from xena mostly looks like: `TCGA-02-0002-01`
#' ``
#'
#' @export
#' @importFrom stringr str_match
#' @param x a character vector of TCGA barcodes
#' @examples
#' barcodes <- c("TCGA-CV-5976-11", "TCGA-EJ-5532-01", "TCGA-D3-A1Q7-06",
#'               "TCGA-06-0211-02", "TCGA-AB-2989-03R")
#' parsed <- parse_tcga_barcode(barcodes)
parse_tcga_barcode <- function(x, ...) {
  assert_character(x, min.chars = 15)
  parsed <- str_match(x, .barcode_regex)
  bad.barcode <- which(is.na(parsed[,1]))
  if (length(bad.barcode)) {
    baddies <- head(bad.barcode, 10)
    msg <- paste(
      length(bad.barcode), "illegal barcodes found at position(s):",
      paste(baddies, collapse = ","))
    if (length(baddies) < length(bad.barcode)) {
      msg <- paste0(msg, "...")
    }
    stop(msg)
  }

  stype <- as.integer(parsed[,4])

  out <- tibble(
    barcode = x,
    patient_id = substr(x, 1, 12),
    sample_type = case_when(
      stype < 10 ~ "tumor",
      stype < 20 ~ "normal",
      TRUE       ~ "control"),
    sample_type_code = parsed[, 4],
    vial = parsed[, 5],
    analyte = .analyte_codes[parsed[, 7]])
  out[["sample_type"]] <- factor(out[["sample_type"]],
                                 c("normal", "tumor", "control"))
  out
}

.barcode_regex <- paste(
  "^TCGA-",
  TSS = "([A-Z0-9]{2})-",
  participant = "([A-Z0-9]{4})-",
  sample = "([0-9]{2})",
  vial = "([A-Z])?-?",
  portion = "(..)?",
  analyte = "([DGHRTWX])?-?",
  rest = "(-.*?$)?", # plate-center
  collapse = "", sep = "")

# this is in inst/extdata/tcga-barcode-codebook.xlsx, but ...
.analyte_codes <- c(
  D = "dna",
  G = "wga_rubicon",
  H = "mirVana",
  R = "rna",
  T = "rna_total",
  W = "wga_qiagen",
  X = "wga_qiagen2")


