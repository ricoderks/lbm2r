library(testthat)
library(lbm2r)

# Path to the test LBM2 file (relative to the package root when using
# testthat::test_package(), or set via env var for standalone runs)
lbm2_path <- test_path("library_test.lbm2")

test_that("read_lbm2 returns a non-empty list", {
  skip_if_not(file.exists(lbm2_path),
              paste("LBM2 test file not found:", lbm2_path))

  records <- read_lbm2(lbm2_path, verbose = FALSE)

  expect_type(records, "list")
  expect_gt(length(records), 0L)
})

test_that("each record has the expected field names", {
  skip_if_not(file.exists(lbm2_path),
              paste("LBM2 test file not found:", lbm2_path))

  records <- read_lbm2(lbm2_path, verbose = FALSE)
  rec     <- records[[1L]]

  expected_fields <- c(
    "scan_id", "precursor_mz", "chrom_xs", "ion_mode", "peaks",
    "name", "formula_bean", "ontology", "smiles", "inchi_key",
    "adduct", "collision_cross_section", "isotopic_peaks", "quant_mass",
    "compound_class", "comment", "instrument", "instrument_type",
    "links", "collision_energy", "database_id", "charge", "ms_level"
  )

  expect_true(all(expected_fields %in% names(rec)),
              info = paste("Missing fields:",
                           paste(setdiff(expected_fields, names(rec)), collapse = ", ")))
})

test_that("ion_mode is 'Positive', 'Negative', or 'Both'", {
  skip_if_not(file.exists(lbm2_path),
              paste("LBM2 test file not found:", lbm2_path))

  records   <- read_lbm2(lbm2_path, verbose = FALSE)
  ion_modes <- vapply(records, `[[`, character(1L), "ion_mode")

  expect_true(all(ion_modes %in% c("Positive", "Negative", "Both")))
})

test_that("peaks field is a data.frame with correct columns", {
  skip_if_not(file.exists(lbm2_path),
              paste("LBM2 test file not found:", lbm2_path))

  records <- read_lbm2(lbm2_path, verbose = FALSE)

  # Check first 100 records for speed
  for (i in seq_len(length(records))) {
    pk <- records[[i]]$peaks
    expect_true(inherits(pk, "data.frame"),
                info = paste("Record", i, "peaks is not a data.frame"))
    expect_true(all(c("mz", "intensity", "comment") %in% names(pk)),
                info = paste("Record", i, "peaks missing columns"))
  }
})

test_that("peak count matches nrow(peaks)", {
  skip_if_not(file.exists(lbm2_path),
              paste("LBM2 test file not found:", lbm2_path))

  records <- read_lbm2(lbm2_path, verbose = FALSE)

  for (i in seq_len(length(records))) {
    rec <- records[[i]]
    expect_equal(length(rec$peaks$mz), nrow(rec$peaks),
                 info = paste("Record", i, "peaks row count inconsistent"))
  }
})

test_that("precursor_mz values are positive", {
  skip_if_not(file.exists(lbm2_path),
              paste("LBM2 test file not found:", lbm2_path))

  records <- read_lbm2(lbm2_path, verbose = FALSE)
  mzs     <- vapply(records, `[[`, double(1L), "precursor_mz")

  expect_true(all(mzs > 0),
              info = paste("Non-positive precursor_mz found in",
                           sum(mzs <= 0), "records"))
})

test_that("adduct field has expected sub-fields", {
  skip_if_not(file.exists(lbm2_path),
              paste("LBM2 test file not found:", lbm2_path))

  records <- read_lbm2(lbm2_path, verbose = FALSE)
  adduct  <- records[[1L]]$adduct

  expect_true(all(c("accurate_mass", "xmer", "name", "charge",
                    "ion_mode", "format_check") %in% names(adduct)))
})
