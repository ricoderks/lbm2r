# write_lbm2.R
# Public API: write_lbm2()

#' Write lipid records to an MSDIAL LBM2 file
#'
#' Serializes a list of lipid records (as returned by \code{\link{read_lbm2}})
#' into a binary LBM2 file that can be loaded by MS-DIAL.  Records are
#' serialized with MessagePack and compressed with LZ4 in the same chunked
#' format used by MS-DIAL itself.
#'
#' @param records A list of lipid records, each as returned by
#'   \code{\link{read_lbm2}}.  All fields must be present; use
#'   \code{\link{new_lbm2_record}} to create a blank template.
#' @param path Character string. Path of the output \code{.lbm2} file.
#'   Existing files are overwritten.
#' @param chunk_size Integer. Maximum number of records per LZ4-compressed
#'   chunk (default \code{50000L}).  Larger values use more memory but produce
#'   fewer chunks.
#' @param verbose Logical. If \code{TRUE} (default) print progress messages.
#'
#' @return Invisibly returns \code{path}.
#'
#' @examples
#' \dontrun{
#' records <- read_lbm2("original.lbm2")
#' # modify records as needed ...
#' write_lbm2(records, "modified.lbm2")
#' }
#'
#' @seealso \code{\link{read_lbm2}}, \code{\link{new_lbm2_record}}
#' @export
write_lbm2 <- function(records, path, chunk_size = 50000L, verbose = TRUE) {
  if (!is.list(records)) stop("'records' must be a list of lipid record lists.")
  chunk_size <- as.integer(chunk_size)
  if (is.na(chunk_size) || chunk_size < 1L) stop("'chunk_size' must be a positive integer.")

  n      <- length(records)
  starts <- seq(1L, max(1L, n), by = chunk_size)
  n_chunks <- length(starts)

  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)

  for (ci in seq_along(starts)) {
    idx_from <- starts[ci]
    idx_to   <- min(idx_from + chunk_size - 1L, n)
    batch    <- records[idx_from:idx_to]

    # Convert each record to a positional list
    raw_list <- lapply(batch, .record_to_raw)

    # Serialize the batch as a MessagePack array
    msgpack_bytes <- RcppMsgPack::msgpack_pack(raw_list)

    # LZ4 compress
    compressed <- lz4_compress(msgpack_bytes)

    # Write 11-byte chunk header + compressed payload
    .write_chunk_header(con,
                        lz4_length        = length(compressed),
                        uncompressed_size = length(msgpack_bytes))
    writeBin(compressed, con)

    if (verbose) {
      message(sprintf("  Chunk %d/%d: records %d-%d (%d bytes -> %d bytes compressed)",
                      ci, n_chunks, idx_from, idx_to,
                      length(msgpack_bytes), length(compressed)))
    }
  }

  if (verbose) {
    message(sprintf("Done. Wrote %d records in %d chunk(s) to: %s", n, n_chunks, path))
  }

  invisible(path)
}


#' Create a blank LBM2 lipid record
#'
#' Returns a named list pre-filled with default (empty/NA) values for every
#' field expected by \code{\link{write_lbm2}}.  Use this as a template when
#' adding new entries to an LBM2 database.
#'
#' @return A named list with all record fields set to sensible defaults,
#'   matching the MoleculeMsReference schema used by MS-DIAL 5.
#'
#' @examples
#' rec <- new_lbm2_record()
#' rec$name           <- "PC 16:0/18:1"
#' rec$compound_class <- "PC"
#' rec$precursor_mz   <- 760.5849
#' rec$ion_mode       <- "Positive"
#' rec$adduct$name    <- "[M+H]+"
#' rec$peaks <- data.frame(
#'   mz               = c(184.0733, 760.5849),
#'   intensity        = c(999, 500),
#'   comment          = c("", ""),
#'   peak_quality     = c(0L, 0L),
#'   peak_id          = c(0L, 1L),
#'   spectrum_comment = c(0L, 0L),
#'   is_required      = c(FALSE, FALSE),
#'   stringsAsFactors = FALSE
#' )
#'
#' @export
new_lbm2_record <- function() {
  list(
    scan_id                    = 0L,
    precursor_mz               = 0.0,
    chrom_xs                   = list(
      rt        = list(),
      ri        = list(),
      drift     = list(),
      mz_chrom  = list(),
      main_type = 0L
    ),
    ion_mode                   = "Positive",
    peaks                      = data.frame(
      mz               = double(0),
      intensity        = double(0),
      comment          = character(0),
      peak_quality     = integer(0),
      peak_id          = integer(0),
      spectrum_comment = integer(0),
      is_required      = logical(0),
      stringsAsFactors = FALSE
    ),
    name                       = "",
    formula_bean               = list(
      formula_string        = "",
      mass                  = 0.0,
      m1_isotope            = 0.0,
      m2_isotope            = 0.0,
      c  = 0L, n  = 0L, h  = 0L, o  = 0L,
      s  = 0L, p  = 0L, f  = 0L, cl = 0L,
      br = 0L, i  = 0L, si = 0L,
      tms_count             = 0L,
      meox_count            = 0L,
      c13 = 0L, n15 = 0L, h2 = 0L, o18 = 0L,
      s34 = 0L, cl37 = 0L, br81 = 0L,
      is_correctly_imported = FALSE,
      se                    = 0L,
      element2count         = list()
    ),
    ontology                   = "",
    smiles                     = "",
    inchi_key                  = "",
    adduct                     = list(
      accurate_mass = 0.0,
      xmer          = 1L,
      name          = "",
      charge        = 1L,
      ion_mode      = "Positive",
      format_check  = FALSE,
      m1_intensity  = 0.0,
      m2_intensity  = 0.0,
      is_radical    = FALSE,
      is_included   = TRUE
    ),
    collision_cross_section    = 0.0,
    isotopic_peaks             = list(),
    quant_mass                 = 0.0,
    compound_class             = "",
    comment                    = "",
    instrument                 = "",
    instrument_type            = "",
    links                      = "",
    collision_energy           = 0.0,
    database_id                = 0L,
    charge                     = 1L,
    ms_level                   = 2L,
    retention_time_tolerance   = 0.05,
    mass_tolerance             = 0.05,
    minimum_peak_height        = 1000.0,
    is_target_molecule         = TRUE,
    database_unique_identifier = "",
    fragmentation_condition    = ""
  )
}
