// lz4_wrapper.cpp
// Thin Rcpp wrappers around the system liblz4 library.
// Exposes lz4_compress() and lz4_decompress() to R.

// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <lz4.h>

using namespace Rcpp;

//' Compress a raw vector using LZ4
//'
//' @param input A \code{raw} vector to compress.
//' @return A \code{raw} vector of LZ4-compressed bytes.
//' @keywords internal
// [[Rcpp::export]]
RawVector lz4_compress(RawVector input) {
  int src_size = input.size();
  int max_dst_size = LZ4_compressBound(src_size);

  RawVector output(max_dst_size);

  int compressed_size = LZ4_compress_default(
    reinterpret_cast<const char*>(RAW(input)),
    reinterpret_cast<char*>(RAW(output)),
    src_size,
    max_dst_size
  );

  if (compressed_size <= 0) {
    stop("LZ4 compression failed");
  }

  output.erase(output.begin() + compressed_size, output.end());
  return output;
}

//' Decompress an LZ4-compressed raw vector
//'
//' @param input A \code{raw} vector of LZ4-compressed bytes.
//' @param uncompressed_size The expected size of the decompressed output in bytes.
//' @return A \code{raw} vector of decompressed bytes.
//' @keywords internal
// [[Rcpp::export]]
RawVector lz4_decompress(RawVector input, int uncompressed_size) {
  RawVector output(uncompressed_size);

  int result = LZ4_decompress_safe(
    reinterpret_cast<const char*>(RAW(input)),
    reinterpret_cast<char*>(RAW(output)),
    input.size(),
    uncompressed_size
  );

  if (result < 0) {
    stop("LZ4 decompression failed (error code %d)", result);
  }
  if (result != uncompressed_size) {
    stop("LZ4 decompression produced %d bytes, expected %d", result, uncompressed_size);
  }

  return output;
}
