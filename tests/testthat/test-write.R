library(testthat)
library(lbm2r)

lbm2_path <- test_path("library_test.lbm2")

# ---------------------------------------------------------------------------
# Helper: create a small list of synthetic records for write-only tests
# ---------------------------------------------------------------------------
make_test_records <- function(n = 5L) {
  lapply(seq_len(n), function(i) {
    rec                    <- new_lbm2_record()
    rec$scan_id            <- i - 1L
    rec$name               <- paste0("TestLipid_", i)
    rec$compound_class     <- "PC"
    rec$precursor_mz       <- 700.0 + i
    rec$ion_mode           <- if (i %% 2 == 0) "Negative" else "Positive"
    rec$smiles             <- "C"
    rec$inchi_key          <- "AAAAAAAAAAAAAA-AAAAAAAAAA-A"
    rec$adduct$name        <- "[M+H]+"
    rec$adduct$accurate_mass <- 1.0073
    rec$adduct$charge      <- 1L
    rec$adduct$xmer        <- 1L
    rec$adduct$ion_mode    <- "Positive"
    rec$adduct$format_check <- TRUE
    rec$formula_bean$formula_string <- "C42H82NO8P"
    rec$formula_bean$c  <- 42L
    rec$formula_bean$h  <- 82L
    rec$formula_bean$n  <- 1L
    rec$formula_bean$o  <- 8L
    rec$formula_bean$p  <- 1L
    rec$peaks <- data.frame(
      mz               = c(184.07, 496.34, 700.0 + i),
      intensity        = c(999.0, 500.0, 200.0),
      comment          = c("head", "NL", "precursor"),
      peak_quality     = c(0L, 0L, 0L),
      peak_id          = c(0L, 1L, 2L),
      spectrum_comment = c(0L, 0L, 0L),
      is_required      = c(FALSE, FALSE, FALSE),
      stringsAsFactors = FALSE
    )
    rec
  })
}


# ---------------------------------------------------------------------------
# new_lbm2_record tests
# ---------------------------------------------------------------------------
test_that("new_lbm2_record returns a list with all required fields", {
  rec <- new_lbm2_record()

  expected_fields <- c(
    "scan_id", "precursor_mz", "chrom_xs", "ion_mode", "peaks",
    "name", "formula_bean", "ontology", "smiles", "inchi_key",
    "adduct", "collision_cross_section", "isotopic_peaks", "quant_mass",
    "compound_class", "comment", "instrument", "instrument_type",
    "links", "collision_energy", "database_id", "charge", "ms_level",
    "retention_time_tolerance", "mass_tolerance", "minimum_peak_height",
    "is_target_molecule"
  )

  expect_true(all(expected_fields %in% names(rec)),
              info = paste("Missing fields:",
                           paste(setdiff(expected_fields, names(rec)), collapse = ", ")))
  expect_s3_class(rec$peaks, "data.frame")
  expect_equal(rec$ion_mode, "Positive")
})


# ---------------------------------------------------------------------------
# write_lbm2 + read_lbm2 round-trip tests (synthetic records, no real file)
# ---------------------------------------------------------------------------
test_that("write_lbm2 creates a file", {
  tmp  <- tempfile(fileext = ".lbm2")
  recs <- make_test_records(3L)
  write_lbm2(recs, tmp, verbose = FALSE)
  expect_true(file.exists(tmp))
  expect_gt(file.size(tmp), 0L)
  unlink(tmp)
})

test_that("round-trip preserves record count", {
  tmp  <- tempfile(fileext = ".lbm2")
  recs <- make_test_records(5L)
  write_lbm2(recs, tmp, verbose = FALSE)
  back <- read_lbm2(tmp, verbose = FALSE)
  expect_equal(length(back), length(recs))
  unlink(tmp)
})

test_that("round-trip preserves record names and scalar fields", {
  tmp  <- tempfile(fileext = ".lbm2")
  recs <- make_test_records(5L)
  write_lbm2(recs, tmp, verbose = FALSE)
  back <- read_lbm2(tmp, verbose = FALSE)

  for (i in seq_along(recs)) {
    expect_equal(back[[i]]$name,           recs[[i]]$name,           info = paste("name, record", i))
    expect_equal(back[[i]]$compound_class, recs[[i]]$compound_class, info = paste("compound_class, record", i))
    expect_equal(back[[i]]$ion_mode,       recs[[i]]$ion_mode,       info = paste("ion_mode, record", i))
    expect_equal(back[[i]]$precursor_mz,   recs[[i]]$precursor_mz,   tolerance = 1e-4,
                 info = paste("precursor_mz, record", i))
  }
  unlink(tmp)
})

test_that("round-trip preserves peaks data.frame", {
  tmp  <- tempfile(fileext = ".lbm2")
  recs <- make_test_records(3L)
  write_lbm2(recs, tmp, verbose = FALSE)
  back <- read_lbm2(tmp, verbose = FALSE)

  for (i in seq_along(recs)) {
    orig_pk <- recs[[i]]$peaks
    back_pk <- back[[i]]$peaks
    expect_equal(nrow(back_pk), nrow(orig_pk), info = paste("nrow(peaks), record", i))
    expect_equal(back_pk$mz,        orig_pk$mz,        tolerance = 1e-4,
                 info = paste("mz, record", i))
    expect_equal(back_pk$intensity, orig_pk$intensity,  tolerance = 1e-2,
                 info = paste("intensity, record", i))
    expect_equal(back_pk$comment,   orig_pk$comment,   info = paste("comment, record", i))
  }
  unlink(tmp)
})

test_that("round-trip preserves adduct fields", {
  tmp  <- tempfile(fileext = ".lbm2")
  recs <- make_test_records(2L)
  write_lbm2(recs, tmp, verbose = FALSE)
  back <- read_lbm2(tmp, verbose = FALSE)

  for (i in seq_along(recs)) {
    expect_equal(back[[i]]$adduct$name,   recs[[i]]$adduct$name,   info = paste("adduct name, record", i))
    expect_equal(back[[i]]$adduct$charge, recs[[i]]$adduct$charge, info = paste("adduct charge, record", i))
  }
  unlink(tmp)
})

test_that("multi-chunk round-trip works correctly", {
  tmp  <- tempfile(fileext = ".lbm2")
  recs <- make_test_records(7L)
  # Force 3 chunks of size 3, 3, 1
  write_lbm2(recs, tmp, chunk_size = 3L, verbose = FALSE)
  back <- read_lbm2(tmp, verbose = FALSE)
  expect_equal(length(back), length(recs))
  expect_equal(vapply(back, `[[`, character(1L), "name"),
               vapply(recs, `[[`, character(1L), "name"))
  unlink(tmp)
})


# ---------------------------------------------------------------------------
# Round-trip against the real LBM2 file (skipped if not present)
# ---------------------------------------------------------------------------
test_that("round-trip on real LBM2 file preserves first 50 records", {
  skip_if_not(file.exists(lbm2_path),
              paste("LBM2 test file not found:", lbm2_path))

  records <- read_lbm2(lbm2_path, verbose = FALSE)
  subset  <- records[seq_len(min(50L, length(records)))]

  tmp <- tempfile(fileext = ".lbm2")
  write_lbm2(subset, tmp, verbose = FALSE)
  back <- read_lbm2(tmp, verbose = FALSE)
  unlink(tmp)

  expect_equal(length(back), length(subset))

  for (i in seq_along(subset)) {
    expect_equal(back[[i]]$name,         subset[[i]]$name,         info = paste("name, record", i))
    expect_equal(back[[i]]$precursor_mz, subset[[i]]$precursor_mz, tolerance = 1e-4,
                 info = paste("precursor_mz, record", i))
    expect_equal(nrow(back[[i]]$peaks),  nrow(subset[[i]]$peaks),  info = paste("nrow(peaks), record", i))
  }
})
