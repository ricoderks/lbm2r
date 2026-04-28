# read_lbm2.R
# Public API: read_lbm2()

#' Read an MSDIAL LBM2 lipid database file
#'
#' Parses a binary LBM2 file produced by MS-DIAL and returns its contents as a
#' list of lipid records.  Each record corresponds to one
#' \code{MspFormatCompoundInformationBean} entry in the underlying
#' MessagePack/LZ4 stream.
#'
#' @param path Character string. Path to the \code{.lbm2} file.
#' @param n_max Integer or \code{Inf}. Maximum number of records to return.
#'   Defaults to \code{Inf} (read everything).  Set to e.g. \code{1000L} to
#'   read only the first 1000 entries — useful when the full file exceeds
#'   available memory.
#' @param verbose Logical. If \code{TRUE} (default) print progress messages
#'   (chunk count, total records loaded).
#'
#' @return A named list of lipid records.  Every element is itself a named list
#'   with the following fields:
#'   \describe{
#'     \item{id}{Integer record index (as stored in the file).}
#'     \item{bin_id}{Integer bin identifier.}
#'     \item{name}{Lipid name, e.g. \code{"PC 16:0/18:1"}.}
#'     \item{compound_class}{Lipid class string, e.g. \code{"PC"}.}
#'     \item{formula}{Molecular formula string, e.g. \code{"C42H82NO8P"}.}
#'     \item{formula_bean}{Named list with individual atom counts and exact
#'       mass; see Details.}
#'     \item{precursor_mz}{Precursor m/z value (numeric).}
#'     \item{adduct}{Named list describing the adduct ion; see Details.}
#'     \item{ion_mode}{Character: \code{"Positive"} or \code{"Negative"}.}
#'     \item{retention_time}{Predicted retention time in minutes (numeric).}
#'     \item{retention_index}{Retention index (numeric, \code{-1} if absent).}
#'     \item{drift_time}{Ion mobility drift time (numeric, \code{-1} if absent).}
#'     \item{collision_cross_section}{CCS value in Å² (numeric, \code{-1} if absent).}
#'     \item{smiles}{SMILES string.}
#'     \item{inchi_key}{InChIKey string.}
#'     \item{ontology}{Ontology annotation string.}
#'     \item{comment}{Free-text comment.}
#'     \item{links}{External database links.}
#'     \item{intensity}{Reference intensity (numeric, \code{-1} if absent).}
#'     \item{instrument}{Instrument name string.}
#'     \item{instrument_type}{Instrument type string.}
#'     \item{collision_energy}{Collision energy string.}
#'     \item{quant_mass}{Quantitation mass (numeric, \code{-1} if absent).}
#'     \item{isotope_ratio_list}{Numeric vector of isotope ratios.}
#'     \item{peak_number}{Integer number of MS/MS peaks.}
#'     \item{peaks}{A \code{data.frame} with columns \code{mz} (numeric),
#'       \code{intensity} (numeric), \code{comment} (character), and
#'       \code{frag} (character) — one row per MS/MS fragment ion.}
#'   }
#'
#' @details
#' **\code{formula_bean} fields:**
#' \code{formula_string} (character), \code{mass} (double),
#' \code{m1_isotope}, \code{m2_isotope} (double),
#' \code{c}, \code{n}, \code{h}, \code{o}, \code{s}, \code{p},
#' \code{f}, \code{cl}, \code{br}, \code{i}, \code{si},
#' \code{tms_count}, \code{meox_count} (all integer).
#'
#' **\code{adduct} fields:**
#' \code{name} (character, e.g. \code{"[M+H]+"}),
#' \code{accurate_mass} (double), \code{xmer} (integer),
#' \code{charge} (integer), \code{ion_type} (character:
#' \code{"Positive"}, \code{"Negative"}, or \code{"Unknown"}),
#' \code{format_check} (logical).
#'
#' @examples
#' \dontrun{
#' # Read entire file
#' records <- read_lbm2("path/to/library.lbm2")
#'
#' # Read only the first 1000 entries (low memory)
#' records <- read_lbm2("path/to/library.lbm2", n_max = 1000)
#'
#' length(records)            # number of lipid entries
#' records[[1]]$name          # name of the first entry
#' records[[1]]$peaks         # MS/MS peak table of the first entry
#' }
#'
#' @export
read_lbm2 <- function(path, n_max = Inf, verbose = TRUE) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
  if (!is.numeric(n_max) || length(n_max) != 1L || n_max < 1) {
    stop("'n_max' must be a positive number or Inf")
  }
  n_max <- if (is.infinite(n_max)) .Machine$integer.max else as.integer(n_max)

  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)

  all_records <- list()
  chunk_index <- 0L
  done        <- FALSE

  while (!done) {
    hdr <- .read_chunk_header(con)
    if (is.null(hdr)) break   # EOF

    chunk_index <- chunk_index + 1L

    # Read LZ4-compressed bytes
    compressed <- readBin(con, what = "raw", n = hdr$lz4_length, endian = "little")
    if (length(compressed) < hdr$lz4_length) {
      warning(sprintf("Chunk %d: expected %d compressed bytes, got %d; file may be truncated",
                      chunk_index, hdr$lz4_length, length(compressed)))
      break
    }
    
    # Decompress
    raw_msgpack <- lz4_decompress(compressed, hdr$uncompressed_size)

    # Unpack MessagePack -> R list of positional lists
    chunk_list <- RcppMsgPack::msgpack_unpack(raw_msgpack)

    # Trim chunk if it would exceed n_max
    remaining <- n_max - length(all_records)
    if (length(chunk_list) >= remaining) {
      chunk_list <- chunk_list[seq_len(remaining)]
      done <- TRUE
    }

    # Convert each positional list to a named record
    chunk_records <- lapply(chunk_list, .record_from_raw)
    all_records   <- c(all_records, chunk_records)

    if (verbose) {
      message(sprintf("  Chunk %d: %d records (total so far: %d)",
                      chunk_index, length(chunk_records), length(all_records)))
    }
  }

  if (verbose) {
    message(sprintf("Done. Read %d chunks, %d records total.", chunk_index, length(all_records)))
  }

  all_records
}
