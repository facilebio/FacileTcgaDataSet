#' Base sample covariates
#'
#' Split sample IDs indications, and provide sex, sample_type, and
#' stage annotation
#'
#' @param xena.dir
#' @return a tibble of id,indication,sex,sample_type
parse_base_sample_covariates <- function(xena.dir, ...) {
  if (FALSE) {
    xena.dir <- "~/workspace/data/FacileData/consortia/tcga/xena"
  }
  assert_directory(xena.dir, "r")
  dat.fn <- xena.dir %>%
    file.path("Survival_SupplementalTable_S1_20171025_xena_sp.gz") %>%
    assert_file_exists()
  sdat <- read_tsv(dat.fn)
  # According to xena, we should have covariate info from 12,591 samples
  stopifnot(nrow(sdat) == 12591)

  out <- sdat %>%
    transmute(dataset = `cancer type abbreviation`,
              sample_id = sample,
              patient_id = `_PATIENT`,
              indication = dataset,
              sex = tolower(gender),
              sample_type = parse_tcga_barcode(sample_id)$sample_type)
  droplevels(out)
}
