# convert.R
# Conversion between raw MessagePack positional lists (as returned by
# RcppMsgPack::msgpack_unpack) and human-readable named R lists.
#
# The real LBM2 file uses the MoleculeMsReference class (newer MSDIAL5 format),
# NOT the old MspFormatCompoundInformationBean.
#
# MoleculeMsReference key layout (0-28):
#   0  scan_id                     int
#   1  precursor_mz                double
#   2  chrom_xs                    ChromXs (list)
#   3  ion_mode                    int  (0=Positive, 1=Negative, 2=Both)
#   4  peaks                       List<SpectrumPeak>
#   5  name                        string
#   6  formula_bean                Formula (list of 27)
#   7  ontology                    string
#   8  smiles                      string
#   9  inchi_key                   string
#  10  adduct                      AdductIon (list of 10)
#  11  collision_cross_section     double
#  12  isotopic_peaks              List<IsotopicPeak>
#  13  quant_mass                  double
#  14  compound_class              string
#  15  comment                     string
#  16  instrument                  string
#  17  instrument_type             string
#  18  links                       string
#  19  collision_energy            float
#  20  database_id                 int
#  21  charge                      int
#  22  ms_level                    int
#  23  retention_time_tolerance    float
#  24  mass_tolerance              float
#  25  minimum_peak_height         float
#  26  is_target_molecule          bool
#  27  database_unique_identifier  string
#  28  fragmentation_condition     string
#
# SpectrumPeak key layout (sparse, IgnoreMember slots left as NULL):
#   0  mz             double
#   1  intensity      double
#   2  comment        string
#   3-5 IgnoreMember (NULL)
#   6  peak_quality   int
#   7  peak_id        int
#   8-10 IgnoreMember (NULL)
#  11  spectrum_comment  int (flags)
#  12  is_required    bool
#
# AdductIon key layout (0-9):
#   0  accurate_mass  double
#   1  xmer           int
#   2  name           string
#   3  charge         int
#   4  ion_mode       int  (0=Positive, 1=Negative, 2=Both)
#   5  format_check   bool
#   6  m1_intensity   double
#   7  m2_intensity   double
#   8  is_radical     bool
#   9  is_included    bool
#
# Formula key layout (0-26):
#   0  formula_string     string
#   1  mass               double
#   2  m1_isotope         double
#   3  m2_isotope         double
#   4  c                  int
#   5  n                  int
#   6  h                  int
#   7  o                  int
#   8  s                  int
#   9  p                  int
#  10  f                  int
#  11  cl                 int
#  12  br                 int
#  13  i                  int
#  14  si                 int
#  15  tms_count          int
#  16  meox_count         int
#  17  c13                int
#  18  n15                int
#  19  h2                 int
#  20  o18                int
#  21  s34                int
#  22  cl37               int
#  23  br81               int
#  24  is_correctly_imported  bool
#  25  se                 int
#  26  element2count      map (named int vector)
#
# ChromXs key layout (0-4):
#   0  rt    list(2): [type_int, list(value, type_int, unit_int)]
#   1  ri    list(2)
#   2  drift list(2)
#   3  mz    list(2)
#   4  main_type  int (ChromXType enum)
#
# IsotopicPeak key layout (0-4):
#   0  relative_abundance                 double
#   1  mass                               double
#   2  mass_diff_from_monoisotopic_ion    double
#   3  comment                            string
#   4  ???                                (5th field observed)

.ION_MODE_LEVELS   <- c("Positive", "Negative", "Both")
.MSP_FIELD_COUNT   <- 29L   # keys 0-28 (29 total, note key 27 out-of-order in C# but stored at idx 28)


# ---------------------------------------------------------------------------
# Helper: safely extract element k (0-based) from a positional list
# ---------------------------------------------------------------------------
.get <- function(lst, k, default = NULL) {
  idx <- k + 1L
  if (idx > length(lst)) return(default)
  val <- lst[[idx]]
  if (is.null(val)) default else val
}

# Scalar unwrapper: MessagePack sometimes stores enum/scalar values as a
# 1-element fixarray (a list with one element). Unwrap it if needed.
.unwrap_scalar <- function(v) {
  if (is.list(v) && length(v) == 1L) v[[1L]] else v
}

.get_str  <- function(lst, k) { v <- .unwrap_scalar(.get(lst, k, "")); if (is.null(v)) "" else as.character(v) }
.get_int  <- function(lst, k) { v <- .unwrap_scalar(.get(lst, k, NA_integer_)); if (is.null(v)) NA_integer_ else as.integer(v) }
.get_dbl  <- function(lst, k) { v <- .unwrap_scalar(.get(lst, k, NA_real_));    if (is.null(v)) NA_real_    else as.double(v) }
.get_bool <- function(lst, k) { v <- .unwrap_scalar(.get(lst, k, FALSE));       if (is.null(v)) FALSE       else as.logical(v) }
.get_list <- function(lst, k) { v <- .get(lst, k, list());      if (is.null(v)) list()      else v }


# ---------------------------------------------------------------------------
# SpectrumPeak: raw list -> data.frame row values
# ---------------------------------------------------------------------------
.peaks_from_raw <- function(peak_list) {
  n <- length(peak_list)
  if (n == 0L) {
    return(data.frame(
      mz               = double(0),
      intensity        = double(0),
      comment          = character(0),
      peak_quality     = integer(0),
      peak_id          = integer(0),
      spectrum_comment = integer(0),
      is_required      = logical(0),
      stringsAsFactors = FALSE
    ))
  }
  mz               <- numeric(n)
  intensity        <- numeric(n)
  comment          <- character(n)
  peak_quality     <- integer(n)
  peak_id          <- integer(n)
  spectrum_comment <- integer(n)
  is_required      <- logical(n)
  for (i in seq_len(n)) {
    p                   <- peak_list[[i]]
    mz[i]               <- .get_dbl(p, 0L)
    intensity[i]        <- .get_dbl(p, 1L)
    comment[i]          <- .get_str(p, 2L)
    peak_quality[i]     <- .get_int(p, 6L)
    peak_id[i]          <- .get_int(p, 7L)
    spectrum_comment[i] <- .get_int(p, 11L)
    is_required[i]      <- .get_bool(p, 12L)
  }
  data.frame(
    mz = mz, intensity = intensity, comment = comment,
    peak_quality = peak_quality, peak_id = peak_id,
    spectrum_comment = spectrum_comment, is_required = is_required,
    stringsAsFactors = FALSE
  )
}

# SpectrumPeak: data.frame -> positional list (13-element, sparse)
.peaks_to_raw <- function(df) {
  lapply(seq_len(nrow(df)), function(i) {
    p <- vector("list", 13L)
    p[[1]]  <- df$mz[i]               # key 0
    p[[2]]  <- df$intensity[i]        # key 1
    p[[3]]  <- df$comment[i]          # key 2
    # keys 3-5: IgnoreMember -> NULL (already NULL)
    p[[7]]  <- as.integer(df$peak_quality[i])     # key 6
    p[[8]]  <- as.integer(df$peak_id[i])           # key 7
    # keys 8-10: IgnoreMember -> NULL
    p[[12]] <- as.integer(df$spectrum_comment[i])  # key 11
    p[[13]] <- as.logical(df$is_required[i])       # key 12
    p
  })
}


# ---------------------------------------------------------------------------
# AdductIon: raw list -> named list
# ---------------------------------------------------------------------------
.adduct_from_raw <- function(lst) {
  # lst may be wrapped in a 1-element list (list of adducts with one entry)
  if (is.list(lst) && length(lst) == 1L && is.list(lst[[1L]])) {
    lst <- lst[[1L]]
  }
  ion_mode_int <- .get_int(lst, 4L)
  list(
    accurate_mass = .get_dbl(lst,  0L),
    xmer          = .get_int(lst,  1L),
    name          = .get_str(lst,  2L),
    charge        = .get_int(lst,  3L),
    ion_mode      = if (!is.na(ion_mode_int) && ion_mode_int >= 0L && ion_mode_int <= 2L)
                      .ION_MODE_LEVELS[ion_mode_int + 1L] else "Positive",
    format_check  = .get_bool(lst, 5L),
    m1_intensity  = .get_dbl(lst,  6L),
    m2_intensity  = .get_dbl(lst,  7L),
    is_radical    = .get_bool(lst, 8L),
    is_included   = .get_bool(lst, 9L)
  )
}

# AdductIon: named list -> positional list (10-element)
.adduct_to_raw <- function(a) {
  ion_mode_int <- match(a$ion_mode, .ION_MODE_LEVELS) - 1L
  if (is.na(ion_mode_int)) ion_mode_int <- 0L
  list(
    a$accurate_mass,          # key 0
    a$xmer,                   # key 1
    a$name,                   # key 2
    a$charge,                 # key 3
    ion_mode_int,             # key 4
    a$format_check,           # key 5
    a$m1_intensity,           # key 6
    a$m2_intensity,           # key 7
    a$is_radical,             # key 8
    a$is_included             # key 9
  )
}


# ---------------------------------------------------------------------------
# Formula: raw list -> named list
# ---------------------------------------------------------------------------
.formula_from_raw <- function(lst) {
  # element2count at key 26 is a msgpack map returned as data.frame with
  # columns 'key' and 'value' by RcppMsgPack, or as a named list
  e2c_raw <- .get(lst, 26L, list())
  if (inherits(e2c_raw, "data.frame") && all(c("key", "value") %in% names(e2c_raw))) {
    keys   <- sapply(e2c_raw$key,   as.character)
    values <- sapply(e2c_raw$value, as.integer)
    e2c <- stats::setNames(as.list(values), keys)
  } else if (is.list(e2c_raw) && !is.null(names(e2c_raw))) {
    e2c <- as.list(e2c_raw)
  } else {
    e2c <- list()
  }
  list(
    formula_string          = .get_str(lst,   0L),
    mass                    = .get_dbl(lst,   1L),
    m1_isotope              = .get_dbl(lst,   2L),
    m2_isotope              = .get_dbl(lst,   3L),
    c                       = .get_int(lst,   4L),
    n                       = .get_int(lst,   5L),
    h                       = .get_int(lst,   6L),
    o                       = .get_int(lst,   7L),
    s                       = .get_int(lst,   8L),
    p                       = .get_int(lst,   9L),
    f                       = .get_int(lst,  10L),
    cl                      = .get_int(lst,  11L),
    br                      = .get_int(lst,  12L),
    i                       = .get_int(lst,  13L),
    si                      = .get_int(lst,  14L),
    tms_count               = .get_int(lst,  15L),
    meox_count              = .get_int(lst,  16L),
    c13                     = .get_int(lst,  17L),
    n15                     = .get_int(lst,  18L),
    h2                      = .get_int(lst,  19L),
    o18                     = .get_int(lst,  20L),
    s34                     = .get_int(lst,  21L),
    cl37                    = .get_int(lst,  22L),
    br81                    = .get_int(lst,  23L),
    is_correctly_imported   = .get_bool(lst, 24L),
    se                      = .get_int(lst,  25L),
    element2count           = e2c
  )
}

# Formula: named list -> positional list (27-element)
.formula_to_raw <- function(f) {
  # element2count as msgpack map: use a named vector -> RcppMsgPack will
  # serialise it as a map
  e2c <- f$element2count
  if (is.null(e2c) || length(e2c) == 0L) {
    e2c_raw <- stats::setNames(integer(0), character(0))
  } else {
    e2c_raw <- stats::setNames(as.integer(unlist(e2c)), names(e2c))
  }
  list(
    f$formula_string,          # key 0
    f$mass,                    # key 1
    f$m1_isotope,              # key 2
    f$m2_isotope,              # key 3
    as.integer(f$c),           # key 4
    as.integer(f$n),           # key 5
    as.integer(f$h),           # key 6
    as.integer(f$o),           # key 7
    as.integer(f$s),           # key 8
    as.integer(f$p),           # key 9
    as.integer(f$f),           # key 10
    as.integer(f$cl),          # key 11
    as.integer(f$br),          # key 12
    as.integer(f$i),           # key 13
    as.integer(f$si),          # key 14
    as.integer(f$tms_count),   # key 15
    as.integer(f$meox_count),  # key 16
    as.integer(f$c13),         # key 17
    as.integer(f$n15),         # key 18
    as.integer(f$h2),          # key 19
    as.integer(f$o18),         # key 20
    as.integer(f$s34),         # key 21
    as.integer(f$cl37),        # key 22
    as.integer(f$br81),        # key 23
    as.logical(f$is_correctly_imported),  # key 24
    as.integer(f$se),          # key 25
    e2c_raw                    # key 26
  )
}


# ---------------------------------------------------------------------------
# ChromXs: raw list -> named list (retain as opaque list for round-trip)
# ---------------------------------------------------------------------------
.chrom_xs_from_raw <- function(lst) {
  # ChromXs has keys 0-4.  Keys 0-3 each hold a [type_int, [value, type, unit]]
  # pair serialised by the interface properties. Key 4 is main_type (int).
  # We store these opaquely for round-trip fidelity.
  list(
    rt        = .get(lst, 0L, list()),
    ri        = .get(lst, 1L, list()),
    drift     = .get(lst, 2L, list()),
    mz_chrom  = .get(lst, 3L, list()),
    main_type = .get_int(lst, 4L)
  )
}

.chrom_xs_to_raw <- function(cx) {
  list(
    cx$rt,                       # key 0
    cx$ri,                       # key 1
    cx$drift,                    # key 2
    cx$mz_chrom,                 # key 3
    as.integer(cx$main_type)     # key 4
  )
}


# ---------------------------------------------------------------------------
# IsotopicPeak: raw list -> named list
# ---------------------------------------------------------------------------
.isotopic_peaks_from_raw <- function(ip_list) {
  lapply(ip_list, function(p) {
    list(
      relative_abundance               = .get_dbl(p, 0L),
      mass                             = .get_dbl(p, 1L),
      mass_diff_from_monoisotopic_ion  = .get_dbl(p, 2L),
      comment                          = .get_str(p, 3L),
      extra                            = .get(p, 4L, NULL)
    )
  })
}

.isotopic_peaks_to_raw <- function(ip_list) {
  lapply(ip_list, function(ip) {
    list(
      ip$relative_abundance,
      ip$mass,
      ip$mass_diff_from_monoisotopic_ion,
      ip$comment,
      ip$extra
    )
  })
}


# ---------------------------------------------------------------------------
# MoleculeMsReference: raw positional list -> named R list
# ---------------------------------------------------------------------------
#' @keywords internal
.record_from_raw <- function(lst) {
  ion_mode_int <- .get_int(lst, 3L)
  ion_mode_str <- if (!is.na(ion_mode_int) && ion_mode_int >= 0L && ion_mode_int <= 2L)
                    .ION_MODE_LEVELS[ion_mode_int + 1L] else "Positive"

  list(
    scan_id                    = .get_int(lst,  0L),
    precursor_mz               = .get_dbl(lst,  1L),
    chrom_xs                   = .chrom_xs_from_raw(.get_list(lst, 2L)),
    ion_mode                   = ion_mode_str,
    peaks                      = .peaks_from_raw(.get_list(lst, 4L)),
    name                       = .get_str(lst,  5L),
    formula_bean               = .formula_from_raw(.get_list(lst, 6L)),
    ontology                   = .get_str(lst,  7L),
    smiles                     = .get_str(lst,  8L),
    inchi_key                  = .get_str(lst,  9L),
    adduct                     = .adduct_from_raw(.get_list(lst, 10L)),
    collision_cross_section    = .get_dbl(lst, 11L),
    isotopic_peaks             = .isotopic_peaks_from_raw(.get_list(lst, 12L)),
    quant_mass                 = .get_dbl(lst, 13L),
    compound_class             = .get_str(lst, 14L),
    comment                    = .get_str(lst, 15L),
    instrument                 = .get_str(lst, 16L),
    instrument_type            = .get_str(lst, 17L),
    links                      = .get_str(lst, 18L),
    collision_energy           = .get_dbl(lst, 19L),
    database_id                = .get_int(lst, 20L),
    charge                     = .get_int(lst, 21L),
    ms_level                   = .get_int(lst, 22L),
    retention_time_tolerance   = .get_dbl(lst, 23L),
    mass_tolerance             = .get_dbl(lst, 24L),
    minimum_peak_height        = .get_dbl(lst, 25L),
    is_target_molecule         = .get_bool(lst, 26L),
    database_unique_identifier = .get_str(lst, 27L),
    fragmentation_condition    = .get_str(lst, 28L)
  )
}


# ---------------------------------------------------------------------------
# MoleculeMsReference: named R list -> raw positional list
# ---------------------------------------------------------------------------
#' @keywords internal
.record_to_raw <- function(rec) {
  ion_mode_int <- match(rec$ion_mode, .ION_MODE_LEVELS) - 1L
  if (is.na(ion_mode_int)) ion_mode_int <- 0L

  # Build 29-element list (keys 0-28)
  r <- vector("list", 29L)
  r[[1]]  <- as.integer(rec$scan_id)                  # key 0
  r[[2]]  <- as.double(rec$precursor_mz)              # key 1
  r[[3]]  <- .chrom_xs_to_raw(rec$chrom_xs)           # key 2
  r[[4]]  <- as.integer(ion_mode_int)                 # key 3
  r[[5]]  <- .peaks_to_raw(rec$peaks)                 # key 4
  r[[6]]  <- rec$name                                 # key 5
  r[[7]]  <- .formula_to_raw(rec$formula_bean)        # key 6
  r[[8]]  <- rec$ontology                             # key 7
  r[[9]]  <- rec$smiles                               # key 8
  r[[10]] <- rec$inchi_key                            # key 9
  r[[11]] <- .adduct_to_raw(rec$adduct)               # key 10
  r[[12]] <- as.double(rec$collision_cross_section)   # key 11
  r[[13]] <- .isotopic_peaks_to_raw(rec$isotopic_peaks)  # key 12
  r[[14]] <- as.double(rec$quant_mass)                # key 13
  r[[15]] <- rec$compound_class                       # key 14
  r[[16]] <- rec$comment                              # key 15
  r[[17]] <- rec$instrument                           # key 16
  r[[18]] <- rec$instrument_type                      # key 17
  r[[19]] <- rec$links                                # key 18
  r[[20]] <- as.double(rec$collision_energy)          # key 19
  r[[21]] <- as.integer(rec$database_id)              # key 20
  r[[22]] <- as.integer(rec$charge)                   # key 21
  r[[23]] <- as.integer(rec$ms_level)                 # key 22
  r[[24]] <- as.double(rec$retention_time_tolerance)  # key 23
  r[[25]] <- as.double(rec$mass_tolerance)            # key 24
  r[[26]] <- as.double(rec$minimum_peak_height)       # key 25
  r[[27]] <- as.logical(rec$is_target_molecule)       # key 26
  r[[28]] <- rec$database_unique_identifier           # key 27
  r[[29]] <- rec$fragmentation_condition              # key 28
  r
}
