# lbm2r

An R package for reading and writing **MSDIAL LBM2** lipid database files.

The LBM2 format is the binary in-silico MS/MS spectral library used by
[MS-DIAL](https://systemsomicslab.github.io/compms/msdial/main.html).
Internally it is a sequence of LZ4-compressed MessagePack chunks, each
containing a batch of serialized lipid records.

## Installation

### System requirement

The package links against the system LZ4 library. Install it before
installing the R package:

```bash
# Debian / Ubuntu
sudo apt-get install liblz4-dev

# Fedora / RHEL
sudo dnf install lz4-devel

# macOS (Homebrew)
brew install lz4
```

### Install from source

```r
# Install dependencies first
install.packages(c("Rcpp", "RcppMsgPack"))

# Install lbm2r (replace path with the actual location)
install.packages("path/to/lbm2r", repos = NULL, type = "source")
```

Or, if the package is on GitHub:

```r
# install.packages("remotes")
remotes::install_github("your-org/lbm2r")
```

## Usage

### Read an LBM2 file

```r
library(lbm2r)

records <- read_lbm2("path/to/library.lbm2")

length(records)          # number of lipid entries
records[[1]]$name        # "PC 34:1"
records[[1]]$ion_mode    # "Positive"
records[[1]]$precursor_mz
records[[1]]$peaks       # data.frame: mz / intensity / comment / frag
```

### Modify records and write back

```r
# Change the retention time of all PC entries
for (i in seq_along(records)) {
  if (records[[i]]$compound_class == "PC") {
    records[[i]]$retention_time <- records[[i]]$retention_time * 1.05
  }
}

write_lbm2(records, "modified.lbm2")
```

### Add a new lipid entry

```r
rec <- new_lbm2_record()
rec$name           <- "PC 16:0/18:1"
rec$compound_class <- "PC"
rec$formula        <- "C42H82NO8P"
rec$precursor_mz   <- 760.5849
rec$ion_mode       <- "Positive"
rec$adduct$name    <- "[M+H]+"
rec$adduct$charge  <- 1L
rec$retention_time <- 12.3
rec$smiles         <- "CCCCCCCCCCCCCCCC(=O)OC..."
rec$inchi_key      <- "XXXXXXXXXXXXXX-YYYYYYYYYY-Z"
rec$peaks <- data.frame(
  mz        = c(184.0733, 496.3402, 760.5849),
  intensity = c(999, 500, 200),
  comment   = c("head group", "NL183", "precursor"),
  frag      = c("", "", ""),
  stringsAsFactors = FALSE
)
rec$peak_number <- nrow(rec$peaks)

records <- c(records, list(rec))
write_lbm2(records, "extended.lbm2")
```

## Record structure

Each element of the list returned by `read_lbm2()` is a named list:

| Field | Type | Description |
|---|---|---|
| `id` | integer | Record index |
| `bin_id` | integer | Bin identifier |
| `name` | character | Lipid name |
| `compound_class` | character | Lipid class (e.g. `"PC"`) |
| `formula` | character | Molecular formula |
| `formula_bean` | named list | Atom counts + exact mass (17 fields) |
| `precursor_mz` | numeric | Precursor m/z |
| `adduct` | named list | Adduct ion info (6 fields) |
| `ion_mode` | character | `"Positive"` or `"Negative"` |
| `retention_time` | numeric | Predicted RT (min) |
| `retention_index` | numeric | Retention index (`-1` if absent) |
| `drift_time` | numeric | Ion mobility drift time (`-1` if absent) |
| `collision_cross_section` | numeric | CCS in Å² (`-1` if absent) |
| `smiles` | character | SMILES string |
| `inchi_key` | character | InChIKey |
| `ontology` | character | Ontology string |
| `comment` | character | Free-text comment |
| `links` | character | External DB links |
| `intensity` | numeric | Reference intensity |
| `instrument` | character | Instrument name |
| `instrument_type` | character | Instrument type |
| `collision_energy` | character | Collision energy |
| `quant_mass` | numeric | Quantitation mass |
| `isotope_ratio_list` | numeric vector | Isotope ratios |
| `peak_number` | integer | Number of MS/MS peaks |
| `peaks` | data.frame | MS/MS peaks: `mz`, `intensity`, `comment`, `frag` |

## Running the tests

```r
# From within the package directory
devtools::test()

# To also run tests against your own LBM2 file, set the env var:
Sys.setenv(LBM2_TEST_FILE = "path/to/your/library.lbm2")
devtools::test()
```

## LBM2 file format reference

The LBM2 file is a concatenation of chunks:

```
[Ext32 header: 6 bytes] [Int32 uncompressed size: 5 bytes] [LZ4 payload: N bytes]
[Ext32 header: 6 bytes] [Int32 uncompressed size: 5 bytes] [LZ4 payload: N bytes]
...
```

Each LZ4 payload decompresses to a MessagePack array of
`MspFormatCompoundInformationBean` objects serialized as positional arrays
(keys 0-25). The source of truth is the MS-DIAL C# class at
`src/Common/CommonStandard/Database/MspFormatCompoundInformationBean.cs`
in the [MsdialWorkbench](https://github.com/systemsomicslab/MsdialWorkbench)
repository.

## License

MIT
