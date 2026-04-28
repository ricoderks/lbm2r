# hdf5.R
# Convert between LBM2 and HDF5 formats.
#
# Design: single-pass, extensible (resizable) HDF5 datasets.
# The LBM2 file is read one LBM2-chunk at a time; each MessagePack chunk is
# unpacked, its records are converted to flat column vectors, and those vectors
# are appended to the HDF5 datasets.  Peak data uses a CSR-style layout:
# a flat mz/intensity/etc. array plus an `offsets` array of length n_records+1.
#
# HDF5 layout
# -----------
# /records/
#   scan_id                 int[n]
#   precursor_mz            double[n]
#   ion_mode                int[n]        0=Positive 1=Negative 2=Both
#   name                    string[n]
#   ontology                string[n]
#   smiles                  string[n]
#   inchi_key               string[n]
#   compound_class          string[n]
#   comment                 string[n]
#   collision_cross_section double[n]
#   quant_mass              double[n]
#   collision_energy        double[n]
#   retention_time          double[n]     value from chrom_xs$rt
#   database_id             int[n]
#   charge                  int[n]
#   ms_level                int[n]
#   adduct_name             string[n]
#   adduct_accurate_mass    double[n]
#   adduct_charge           int[n]
#   adduct_xmer             int[n]
#   adduct_ion_mode         int[n]
#   formula_string          string[n]
#   formula_mass            double[n]
# /peaks/
#   mz               double[total_peaks]
#   intensity        double[total_peaks]
#   comment          string[total_peaks]
#   spectrum_comment int[total_peaks]
#   offsets          int[n+1]    # peaks for record i = offsets[i]:(offsets[i+1]-1) (0-based)


# ---------------------------------------------------------------------------
# Internal helper: extract the numeric RT value from a parsed ChromXs list
# ---------------------------------------------------------------------------
.rt_from_chrom_xs <- function(cx) {
  rt_raw <- cx$rt
  if (is.list(rt_raw) && length(rt_raw) >= 2L) {
    inner <- rt_raw[[2L]]
    if (is.list(inner) && length(inner) >= 1L) {
      v <- inner[[1L]]
      if (is.numeric(v)) return(as.double(v))
    }
  }
  NA_real_
}


# ---------------------------------------------------------------------------
# Internal: create an extensible 1-D HDF5 dataset starting at size 0
# ---------------------------------------------------------------------------
.h5_create_ext <- function(grp, name, type, chunk = 10000L) {
  dtype <- switch(type,
    int    = hdf5r::h5types$H5T_NATIVE_INT32,
    int64  = hdf5r::h5types$H5T_NATIVE_INT64,
    double = hdf5r::h5types$H5T_NATIVE_DOUBLE,
    string = hdf5r::H5T_STRING$new(size = Inf)
  )
  sp <- hdf5r::H5S$new("simple", dims = 0L, maxdims = Inf)
  grp$create_dataset(name, dtype = dtype, space = sp, chunk_dims = chunk)
}


# ---------------------------------------------------------------------------
# Internal: append a vector to an extensible 1-D HDF5 dataset
# ---------------------------------------------------------------------------
.h5_append <- function(ds, values) {
  if (length(values) == 0L) return(invisible(NULL))
  cur      <- ds$dims
  new_size <- cur + length(values)
  ds$set_extent(new_size)
  idx <- seq.int(cur + 1L, new_size)
  ds[idx] <- values
}


#' Convert an LBM2 file to HDF5
#'
#' Reads an LBM2 file one LZ4-compressed chunk at a time and writes a flat
#' columnar HDF5 file.  Memory use is bounded to one decompressed+unpacked
#' chunk at a time (typically ~1.7 GB for 100 000 records).
#'
#' Once the HDF5 file exists, retention times (and any other scalar fields) can
#' be modified in place with a few lines of `hdf5r` code, and the file can be
#' converted back to LBM2 with \code{\link{hdf5_to_lbm2}}.
#'
#' @param lbm2_path Character. Path to the source \code{.lbm2} file.
#' @param h5_path   Character. Path for the output \code{.h5} file.
#'   Existing files are overwritten.
#' @param chunk_size Integer. Number of records (MessagePack items) to process
#'   per HDF5 write batch.  Does not limit memory beyond what one LBM2 chunk
#'   occupies; lower values just flush to HDF5 more frequently (default
#'   \code{50000L}).
#' @param verbose Logical. Print progress messages (default \code{TRUE}).
#'
#' @return Invisibly returns \code{h5_path}.
#'
#' @examples
#' \dontrun{
#' lbm2_to_hdf5("library.lbm2", "library.h5")
#'
#' # Adjust retention times (5 % correction)
#' library(hdf5r)
#' h5 <- hdf5r::H5File$new("library.h5", mode = "r+")
#' rt <- h5[["records/retention_time"]][]
#' h5[["records/retention_time"]][] <- rt * 1.05
#' h5$close_all()
#'
#' # Write corrected LBM2
#' hdf5_to_lbm2("library.h5", "library_corrected.lbm2")
#' }
#'
#' @seealso \code{\link{hdf5_to_lbm2}}
#' @export
lbm2_to_hdf5 <- function(lbm2_path, h5_path, chunk_size = 50000L, verbose = TRUE) {
  if (!requireNamespace("hdf5r", quietly = TRUE))
    stop("Package 'hdf5r' is required. Install it with install.packages('hdf5r').")
  if (!file.exists(lbm2_path)) stop("File not found: ", lbm2_path)
  chunk_size <- as.integer(chunk_size)

  # ---- Create HDF5 file with extensible datasets --------------------------
  if (file.exists(h5_path)) file.remove(h5_path)
  h5 <- hdf5r::H5File$new(h5_path, mode = "w")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)

  rg <- h5$create_group("records")
  pg <- h5$create_group("peaks")

  for (nm in c("scan_id", "ion_mode", "database_id", "charge", "ms_level",
               "adduct_charge", "adduct_xmer", "adduct_ion_mode"))
    .h5_create_ext(rg, nm, "int")
  for (nm in c("precursor_mz", "collision_cross_section", "quant_mass",
               "collision_energy", "retention_time", "adduct_accurate_mass",
               "formula_mass"))
    .h5_create_ext(rg, nm, "double")
  for (nm in c("name", "ontology", "smiles", "inchi_key", "compound_class",
               "comment", "adduct_name", "formula_string"))
    .h5_create_ext(rg, nm, "string")

  .h5_create_ext(pg, "mz",               "double")
  .h5_create_ext(pg, "intensity",        "double")
  .h5_create_ext(pg, "spectrum_comment", "int")
  .h5_create_ext(pg, "comment",          "string")
  # offsets: CSR-style; first element is always 0
  .h5_create_ext(pg, "offsets", "int64")
  .h5_append(pg[["offsets"]], 0L)   # sentinel: offset[0] = 0

  # ---- Single pass: read LBM2 chunks, write HDF5 --------------------------
  con <- file(lbm2_path, open = "rb")
  on.exit(try(close(con), silent = TRUE), add = TRUE)

  total_records <- 0L
  total_peaks   <- 0L
  lbm2_chunk_idx <- 0L

  repeat {
    hdr <- .read_chunk_header(con)
    if (is.null(hdr)) break
    lbm2_chunk_idx <- lbm2_chunk_idx + 1L

    if (verbose)
      message(sprintf("Decompressing LZ4 chunk %d  (%.0f MB compressed -> %.0f MB)",
                      lbm2_chunk_idx,
                      hdr$lz4_length / 1e6,
                      hdr$uncompressed_size / 1e6))

    compressed  <- readBin(con, "raw", n = hdr$lz4_length)
    raw_msgpack <- lz4_decompress(compressed, hdr$uncompressed_size)
    rm(compressed); gc(verbose = FALSE)

    if (verbose)
      message(sprintf("  Unpacking MessagePack for LZ4 chunk %d...", lbm2_chunk_idx))
    chunk_list <- RcppMsgPack::msgpack_unpack(raw_msgpack)
    rm(raw_msgpack); gc(verbose = FALSE)

    n_all <- length(chunk_list)

    # Process in sub-batches of chunk_size to limit intermediate R objects
    batch_starts <- seq(1L, n_all, by = chunk_size)
    for (bs in batch_starts) {
      be  <- min(bs + chunk_size - 1L, n_all)
      m   <- be - bs + 1L
      sub <- chunk_list[bs:be]

      # Pre-allocate batch vectors
      v_scan_id <- integer(m);    v_pmz   <- numeric(m)
      v_imode   <- integer(m);    v_rt    <- numeric(m)
      v_ccs     <- numeric(m);    v_qm    <- numeric(m)
      v_ce      <- numeric(m);    v_dbid  <- integer(m)
      v_chg     <- integer(m);    v_mslvl <- integer(m)
      v_name    <- character(m);  v_ont   <- character(m)
      v_smi     <- character(m);  v_inchi <- character(m)
      v_cls     <- character(m);  v_cmt   <- character(m)
      v_aname   <- character(m);  v_amass <- numeric(m)
      v_achg    <- integer(m);    v_axmer <- integer(m)
      v_aimode  <- integer(m)
      v_fstr    <- character(m);  v_fmass <- numeric(m)

      pk_mz  <- vector("list", m)
      pk_int <- vector("list", m)
      pk_spc <- vector("list", m)
      pk_cmt <- vector("list", m)

      for (j in seq_len(m)) {
        r <- sub[[j]]

        v_scan_id[j] <- .get_int(r, 0L)
        v_pmz[j]     <- .get_dbl(r, 1L)
        v_imode[j]   <- .get_int(r, 3L)

        cx           <- .get_list(r, 2L)
        v_rt[j]      <- .rt_from_chrom_xs(.chrom_xs_from_raw(cx))

        pk_raw       <- .get_list(r, 4L)
        np_j         <- length(pk_raw)

        if (np_j > 0L) {
          pk_mz[[j]]  <- vapply(pk_raw, function(p) .get_dbl(p, 0L),  numeric(1L))
          pk_int[[j]] <- vapply(pk_raw, function(p) .get_dbl(p, 1L),  numeric(1L))
          pk_spc[[j]] <- vapply(pk_raw, function(p) .get_int(p, 11L), integer(1L))
          pk_cmt[[j]] <- vapply(pk_raw, function(p) .get_str(p, 2L),  character(1L))
        } else {
          pk_mz[[j]]  <- numeric(0);   pk_int[[j]] <- numeric(0)
          pk_spc[[j]] <- integer(0);   pk_cmt[[j]] <- character(0)
        }

        v_name[j]  <- .get_str(r, 5L)
        fb         <- .get_list(r, 6L)
        v_fstr[j]  <- .get_str(fb, 0L)
        v_fmass[j] <- .get_dbl(fb, 1L)
        v_ont[j]   <- .get_str(r, 7L)
        v_smi[j]   <- .get_str(r, 8L)
        v_inchi[j] <- .get_str(r, 9L)
        ad         <- .get_list(r, 10L)
        if (is.list(ad) && length(ad) == 1L && is.list(ad[[1L]])) ad <- ad[[1L]]
        v_amass[j]  <- .get_dbl(ad, 0L)
        v_axmer[j]  <- .get_int(ad, 1L)
        v_aname[j]  <- .get_str(ad, 2L)
        v_achg[j]   <- .get_int(ad, 3L)
        v_aimode[j] <- .get_int(ad, 4L)
        v_ccs[j]    <- .get_dbl(r, 11L)
        v_qm[j]     <- .get_dbl(r, 13L)
        v_cls[j]    <- .get_str(r, 14L)
        v_cmt[j]    <- .get_str(r, 15L)
        v_ce[j]     <- .get_dbl(r, 19L)
        v_dbid[j]   <- .get_int(r, 20L)
        v_chg[j]    <- .get_int(r, 21L)
        v_mslvl[j]  <- .get_int(r, 22L)
      }
      rm(sub); gc(verbose = FALSE)

      # Append record columns
      .h5_append(rg[["scan_id"]],                 v_scan_id)
      .h5_append(rg[["precursor_mz"]],            v_pmz)
      .h5_append(rg[["ion_mode"]],                v_imode)
      .h5_append(rg[["retention_time"]],          v_rt)
      .h5_append(rg[["collision_cross_section"]], v_ccs)
      .h5_append(rg[["quant_mass"]],              v_qm)
      .h5_append(rg[["collision_energy"]],        v_ce)
      .h5_append(rg[["database_id"]],             v_dbid)
      .h5_append(rg[["charge"]],                  v_chg)
      .h5_append(rg[["ms_level"]],                v_mslvl)
      .h5_append(rg[["name"]],                    v_name)
      .h5_append(rg[["ontology"]],                v_ont)
      .h5_append(rg[["smiles"]],                  v_smi)
      .h5_append(rg[["inchi_key"]],               v_inchi)
      .h5_append(rg[["compound_class"]],          v_cls)
      .h5_append(rg[["comment"]],                 v_cmt)
      .h5_append(rg[["adduct_name"]],             v_aname)
      .h5_append(rg[["adduct_accurate_mass"]],    v_amass)
      .h5_append(rg[["adduct_charge"]],           v_achg)
      .h5_append(rg[["adduct_xmer"]],             v_axmer)
      .h5_append(rg[["adduct_ion_mode"]],         v_aimode)
      .h5_append(rg[["formula_string"]],          v_fstr)
      .h5_append(rg[["formula_mass"]],            v_fmass)

      # Append peak data + cumulative offsets
      pk_counts  <- vapply(pk_mz, length, integer(1L))
      cum_offset <- total_peaks + cumsum(pk_counts)
      .h5_append(pg[["offsets"]],          as.integer(cum_offset))
      .h5_append(pg[["mz"]],              unlist(pk_mz,  use.names = FALSE))
      .h5_append(pg[["intensity"]],       unlist(pk_int, use.names = FALSE))
      .h5_append(pg[["spectrum_comment"]],unlist(pk_spc, use.names = FALSE))
      .h5_append(pg[["comment"]],         unlist(pk_cmt, use.names = FALSE))

      total_peaks   <- total_peaks   + sum(pk_counts)
      total_records <- total_records + m

      rm(v_scan_id, v_pmz, v_imode, v_rt, v_ccs, v_qm, v_ce,
         v_dbid, v_chg, v_mslvl, v_name, v_ont, v_smi, v_inchi,
         v_cls, v_cmt, v_aname, v_amass, v_achg, v_axmer, v_aimode,
         v_fstr, v_fmass, pk_mz, pk_int, pk_spc, pk_cmt, pk_counts)
      gc(verbose = FALSE)

      if (verbose)
        message(sprintf("  [LZ4 chunk %d]  records written so far: %d  (peaks: %d)",
                        lbm2_chunk_idx, total_records, total_peaks))
    }

    rm(chunk_list); gc(verbose = FALSE)
  }
  close(con)

  # Write metadata attributes
  h5$create_attr("n_records",
    space = hdf5r::H5S$new("scalar"),
    dtype = hdf5r::h5types$H5T_NATIVE_INT32)$write(as.integer(total_records))
  h5$create_attr("n_peaks",
    space = hdf5r::H5S$new("scalar"),
    dtype = hdf5r::h5types$H5T_NATIVE_INT32)$write(as.integer(total_peaks))
  h5$create_attr("lbm2r_version",
    space = hdf5r::H5S$new("scalar"),
    dtype = hdf5r::H5T_STRING$new(size = Inf))$write("0.1.0")

  h5$close_all()
  if (verbose)
    message(sprintf("Done. %d records, %d peaks -> %s", total_records, total_peaks, h5_path))
  invisible(h5_path)
}


#' Convert an HDF5 file (created by lbm2_to_hdf5) back to LBM2
#'
#' Reads the HDF5 file produced by \code{\link{lbm2_to_hdf5}} chunk-by-chunk
#' and writes a new LBM2 file.  All scalar fields stored in HDF5 are restored;
#' the ChromXs, IsotopicPeaks, and other nested objects that are not stored in
#' HDF5 are reconstructed with sensible defaults.
#'
#' @param h5_path   Character. Path to the source \code{.h5} file.
#' @param lbm2_path Character. Path for the output \code{.lbm2} file.
#'   Existing files are overwritten.
#' @param chunk_size Integer. Number of records to write per LBM2 chunk
#'   (default \code{50000L}).
#' @param verbose Logical. Print progress messages (default \code{TRUE}).
#'
#' @return Invisibly returns \code{lbm2_path}.
#'
#' @seealso \code{\link{lbm2_to_hdf5}}
#' @export
hdf5_to_lbm2 <- function(h5_path, lbm2_path, chunk_size = 50000L, verbose = TRUE) {
  if (!requireNamespace("hdf5r", quietly = TRUE))
    stop("Package 'hdf5r' is required. Install it with install.packages('hdf5r').")
  if (!file.exists(h5_path)) stop("File not found: ", h5_path)
  chunk_size <- as.integer(chunk_size)

  h5 <- hdf5r::H5File$new(h5_path, mode = "r")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)

  rg      <- h5[["records"]]
  pg      <- h5[["peaks"]]
  offsets <- pg[["offsets"]][]   # integer vector length n_records+1 (0-based)
  n_records <- rg[["scan_id"]]$dims

  ion_mode_levels <- c("Positive", "Negative", "Both")

  con_out <- file(lbm2_path, open = "wb")
  on.exit(try(close(con_out), silent = TRUE), add = TRUE)

  starts <- seq(1L, n_records, by = chunk_size)

  for (ci in seq_along(starts)) {
    from <- starts[ci]
    to   <- min(from + chunk_size - 1L, n_records)
    m    <- to - from + 1L
    idx  <- seq(from, to)

    # Read record columns for this batch
    scan_id  <- rg[["scan_id"]][idx]
    pmz      <- rg[["precursor_mz"]][idx]
    ion_mode <- rg[["ion_mode"]][idx]
    rt       <- rg[["retention_time"]][idx]
    ccs      <- rg[["collision_cross_section"]][idx]
    qm       <- rg[["quant_mass"]][idx]
    ce       <- rg[["collision_energy"]][idx]
    dbid     <- rg[["database_id"]][idx]
    chg      <- rg[["charge"]][idx]
    mslvl    <- rg[["ms_level"]][idx]
    name     <- rg[["name"]][idx]
    ont      <- rg[["ontology"]][idx]
    smi      <- rg[["smiles"]][idx]
    inchi    <- rg[["inchi_key"]][idx]
    cls      <- rg[["compound_class"]][idx]
    cmt      <- rg[["comment"]][idx]
    aname    <- rg[["adduct_name"]][idx]
    amass    <- rg[["adduct_accurate_mass"]][idx]
    achg     <- rg[["adduct_charge"]][idx]
    axmer    <- rg[["adduct_xmer"]][idx]
    aimode   <- rg[["adduct_ion_mode"]][idx]
    fstr     <- rg[["formula_string"]][idx]
    fmass    <- rg[["formula_mass"]][idx]

    # Read all peaks for this batch in one slice (offsets are 0-based)
    pk_from <- offsets[from] + 1L     # convert to 1-based R index
    pk_to   <- offsets[to + 1L]       # still 0-based end = last 1-based index
    has_peaks <- pk_to >= pk_from && pk_from >= 1L

    if (has_peaks) {
      pk_idx   <- seq(pk_from, pk_to)
      all_mz   <- pg[["mz"]][pk_idx]
      all_int  <- pg[["intensity"]][pk_idx]
      all_spc  <- pg[["spectrum_comment"]][pk_idx]
      all_pcmt <- pg[["comment"]][pk_idx]
    }

    # Build record list
    records <- vector("list", m)
    for (j in seq_len(m)) {
      rec <- new_lbm2_record()

      rec$scan_id      <- as.integer(scan_id[j])
      rec$precursor_mz <- pmz[j]
      im_int           <- as.integer(ion_mode[j])
      rec$ion_mode     <- if (!is.na(im_int) && im_int >= 0L && im_int <= 2L)
                            ion_mode_levels[im_int + 1L] else "Positive"

      # Rebuild ChromXs with the (possibly modified) RT value
      rt_val       <- rt[j]
      rec$chrom_xs <- list(
        rt        = list(1L, list(rt_val, 0L, 0L)),
        ri        = list(2L, list(-1.0,   1L, 4L)),
        drift     = list(4L, list(-1.0,   2L, 2L)),
        mz_chrom  = list(3L, list(-1.0,   3L, 3L)),
        main_type = 0L
      )

      rec$name                    <- name[j]
      rec$ontology                <- ont[j]
      rec$smiles                  <- smi[j]
      rec$inchi_key               <- inchi[j]
      rec$compound_class          <- cls[j]
      rec$comment                 <- cmt[j]
      rec$collision_cross_section <- ccs[j]
      rec$quant_mass              <- qm[j]
      rec$collision_energy        <- ce[j]
      rec$database_id             <- as.integer(dbid[j])
      rec$charge                  <- as.integer(chg[j])
      rec$ms_level                <- as.integer(mslvl[j])

      rec$adduct$name          <- aname[j]
      rec$adduct$accurate_mass <- amass[j]
      rec$adduct$charge        <- as.integer(achg[j])
      rec$adduct$xmer          <- as.integer(axmer[j])
      aim <- as.integer(aimode[j])
      rec$adduct$ion_mode      <- if (!is.na(aim) && aim >= 0L && aim <= 2L)
                                    ion_mode_levels[aim + 1L] else "Positive"

      rec$formula_bean$formula_string <- fstr[j]
      rec$formula_bean$mass           <- fmass[j]

      # Peaks: offsets are 0-based; convert to 1-based positions within pk_idx
      p_abs_start <- offsets[from + j - 1L] + 1L   # 1-based absolute
      p_abs_end   <- offsets[from + j]              # 1-based absolute (end inclusive)
      if (has_peaks && p_abs_end >= p_abs_start) {
        # Translate to positions within the already-sliced peak vectors
        p_start <- p_abs_start - pk_from + 1L
        p_end   <- p_abs_end   - pk_from + 1L
        p_idx   <- seq(p_start, p_end)
        rec$peaks <- data.frame(
          mz               = all_mz[p_idx],
          intensity        = all_int[p_idx],
          comment          = all_pcmt[p_idx],
          peak_quality     = 0L,
          peak_id          = seq_along(p_idx) - 1L,
          spectrum_comment = as.integer(all_spc[p_idx]),
          is_required      = FALSE,
          stringsAsFactors = FALSE
        )
      }

      records[[j]] <- rec
    }

    # Serialize chunk to LBM2
    raw_list      <- lapply(records, .record_to_raw)
    msgpack_bytes <- RcppMsgPack::msgpack_pack(raw_list)
    compressed    <- lz4_compress(msgpack_bytes)
    .write_chunk_header(con_out,
                        lz4_length        = length(compressed),
                        uncompressed_size = length(msgpack_bytes))
    writeBin(compressed, con_out)

    if (verbose)
      message(sprintf("  Chunk %d/%d: records %d-%d written",
                      ci, length(starts), from, to))
    rm(records, raw_list, msgpack_bytes, compressed,
       scan_id, pmz, ion_mode, rt, ccs, qm, ce, dbid, chg, mslvl,
       name, ont, smi, inchi, cls, cmt, aname, amass, achg, axmer,
       aimode, fstr, fmass)
    if (has_peaks) rm(all_mz, all_int, all_spc, all_pcmt)
    gc(verbose = FALSE)
  }

  close(con_out)
  h5$close_all()
  if (verbose) message(sprintf("Done. LBM2 written to: %s", lbm2_path))
  invisible(lbm2_path)
}
