context("TCGA Barcodes")

test_that("we find meaing in parse_tcga_barcode", {
  dat <- tribble(
    ~barcode,                       ~sample_type, ~analyte,
    "TCGA-CV-5976-11",              "normal",     NA_character_,
    "TCGA-EJ-5532-01",              "tumor",      NA_character_,
    "TCGA-D3-A1Q7-06",              "tumor",      NA_character_,
    "TCGA-06-0211-02A",             "tumor",      NA_character_,
    "TCGA-AB-2989-03B",             "tumor",      NA_character_,
    "TCGA-02-0001-01C-01D-0183-01", "tumor",      "dna")
  dat$sample_type <- factor(dat$sample_type, c("normal", "tumor", "control"))

  parsed <- parse_tcga_barcode(dat$barcode)
  expect_equal(parsed$barcode, dat$barcode)
  expect_equal(parsed$sample_type, dat$sample_type)
  expect_equal(parsed$analyte, dat$analyte)
})
