# chunks.R
# Internal helpers for reading and writing the 11-byte LBM2 chunk headers.
#
# The LBM2 file is a sequence of LZ4-compressed MessagePack chunks. Each chunk
# is preceded by an 11-byte header with the following layout (big-endian):
#
#   Bytes 0-5  : MessagePack Ext32 header
#                  byte 0    = 0xC9  (Ext32 format byte)
#                  bytes 1-4 = uint32 big-endian: LZ4 payload length + 5
#                  byte 5    = 0x63  (type code 99, MSDIAL custom)
#   Bytes 6-9  : int32 big-endian: uncompressed size of the MessagePack data
#   Byte  10   : first byte of the LZ4 payload (already read as part of the
#                5-byte "length prefix" block in the original C# code)
#
# Total: 6 (ext32 header) + 5 (int32 block: 1 format byte + 4 bytes) = 11 bytes
# before the actual LZ4 compressed data begins.

# Expected extension type code written by MSDIAL
.EXT_TYPE_CODE  <- 0x63L   # decimal 99
.EXT32_FORMAT   <- 0xC9L   # MessagePack Ext32 format byte
.INT32_FORMAT   <- 0xD2L   # MessagePack int32 format byte (used for length prefix)

#' Read one chunk header from a binary connection
#'
#' Reads exactly 11 bytes. Returns a named list with:
#' \describe{
#'   \item{lz4_length}{number of LZ4-compressed bytes that follow}
#'   \item{uncompressed_size}{expected decompressed size in bytes}
#' }
#' Returns \code{NULL} at end-of-file.
#' @keywords internal
.read_chunk_header <- function(con) {
  # --- byte 0: Ext32 format byte ---
  b <- readBin(con, what = "raw", n = 1L, endian = "big")
  if (length(b) == 0L) return(NULL)   # EOF

  if (as.integer(b) != .EXT32_FORMAT) {
    stop(sprintf(
      "LBM2 parse error: expected Ext32 format byte 0xC9, got 0x%02X at position %d",
      as.integer(b), seek(con) - 1L
    ))
  }

  # --- bytes 1-4: uint32 big-endian = LZ4 payload length + 5 ---
  len_bytes <- readBin(con, what = "raw", n = 4L, endian = "big")
  # Interpret as unsigned 32-bit big-endian integer
  payload_plus5 <- sum(as.integer(len_bytes) * c(16777216L, 65536L, 256L, 1L))

  # --- byte 5: type code (must be 99 / 0x63) ---
  type_byte <- readBin(con, what = "raw", n = 1L, endian = "big")
  if (as.integer(type_byte) != .EXT_TYPE_CODE) {
    stop(sprintf(
      "LBM2 parse error: expected type code 0x63 (99), got 0x%02X",
      as.integer(type_byte)
    ))
  }

  # --- bytes 6-9: int32 big-endian = uncompressed size ---
  # Written by C# as MessagePackBinary.WriteInt32ForceInt32Block which writes
  # format byte 0xD2 followed by 4 big-endian bytes.
  fmt_byte <- readBin(con, what = "raw", n = 1L, endian = "big")
  if (as.integer(fmt_byte) != .INT32_FORMAT) {
    stop(sprintf(
      "LBM2 parse error: expected Int32 format byte 0xD2, got 0x%02X",
      as.integer(fmt_byte)
    ))
  }
  unc_bytes <- readBin(con, what = "raw", n = 4L, endian = "big")
  uncompressed_size <- sum(as.integer(unc_bytes) * c(16777216L, 65536L, 256L, 1L))

  # The lz4 payload length = (payload_plus5 - 5)
  lz4_length <- payload_plus5 - 5L

  list(lz4_length = lz4_length, uncompressed_size = uncompressed_size)
}


#' Write one chunk header to a binary connection
#'
#' @param con    A writable binary connection.
#' @param lz4_length        Length of the LZ4 payload in bytes.
#' @param uncompressed_size Decompressed size in bytes.
#' @keywords internal
.write_chunk_header <- function(con, lz4_length, uncompressed_size) {
  # byte 0: Ext32 format byte
  writeBin(as.raw(.EXT32_FORMAT), con)

  # bytes 1-4: uint32 big-endian = lz4_length + 5
  val <- lz4_length + 5L
  writeBin(as.raw(c(
    bitwAnd(bitwShiftR(val, 24L), 0xFFL),
    bitwAnd(bitwShiftR(val, 16L), 0xFFL),
    bitwAnd(bitwShiftR(val,  8L), 0xFFL),
    bitwAnd(val,                  0xFFL)
  )), con)

  # byte 5: type code 99
  writeBin(as.raw(.EXT_TYPE_CODE), con)

  # byte 6: Int32 format byte
  writeBin(as.raw(.INT32_FORMAT), con)

  # bytes 7-10: int32 big-endian = uncompressed_size
  writeBin(as.raw(c(
    bitwAnd(bitwShiftR(uncompressed_size, 24L), 0xFFL),
    bitwAnd(bitwShiftR(uncompressed_size, 16L), 0xFFL),
    bitwAnd(bitwShiftR(uncompressed_size,  8L), 0xFFL),
    bitwAnd(uncompressed_size,             0xFFL)
  )), con)
}
