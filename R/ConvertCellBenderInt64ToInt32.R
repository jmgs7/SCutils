#' @title ConvertCellBenderInt64ToInt32
#'
#' @description
#' Some CellBender HDF5 outputs store the count values in `/matrix/data` as
#' 64-bit integers instead of the 32-bit integers used by the standard 10x
#' Genomics HDF5 format. `BPCells::open_matrix_10x_hdf5()` cannot read such
#' files and fails with:
#' \preformatted{
#' get10xMatrixType: unrecognized type for group matrix
#' }
#' This has been observed for some CellBender outputs, while every
#' other `/matrix` dataset (`indices`, `indptr`, `shape`) already
#' uses standard, BPCells-compatible types.
#'
#' `ConvertCellbenderInt64ToInt32()` repairs only the affected dataset. It
#' never touches the input file; it writes a new file next to it (or
#' wherever `output.h5` points) and performs an exhaustive, whole-file
#' checksum comparison before returning, so the caller can trust that no
#' values were altered besides the intended dtype change.
#'
#' @details
#' The function proceeds in four stages, each explained in the inline
#' comments of the function body:
#' \enumerate{
#'   \item \strong{Cheap pre-check.} Read only the HDF5 *datatype* (not the
#'     values) of `/matrix/data`. If it is already 4 bytes (int32),
#'     there is nothing to fix: the function prints an informative
#'     message and returns the input file path as invisible output
#'     \emph{without raising an error}, so it is safe to call unconditionally
#'     inside a loop or a pipeline over many files without special-casing
#'     "already OK" files.
#'   \item \strong{Copy + targeted rewrite.} Copy the whole input file
#'     byte-for-byte (so every group/dataset/attribute other than
#'     `/matrix/data` is preserved exactly), then read just the
#'     `/matrix/data` values, validate that they are safe to store as
#'     32-bit integers, and rewrite that one dataset in place using the
#'     *same* chunk size, gzip level and shuffle filter as the original
#'     dataset (so compression behaviour is not silently changed).
#'   \item \strong{Whole-file verification.} Every dataset in the file,
#'     including ones BPCells never reads (`droplet_latents`,
#'     `global_latents`, `metadata`, ...), is read back from both the input
#'     and the output file and compared via MD5 digest of its decoded
#'     values. Any missing dataset or value mismatch stops execution with a
#'     detailed error naming the offending dataset(s).
#'   \item \strong{Return.} The path to the new file is returned invisibly.
#' }
#'
#' @section Why the shuffle filter does not threaten data integrity:
#' HDF5's byte-shuffle filter is a lossless, fully reversible
#' rearrangement of bytes applied \emph{before} compression. For an array
#' of same-sized numbers, it groups together the 1st byte of every
#' element, then the 2nd byte of every element, and so on; this tends to
#' make the byte stream more compressible by gzip. On read, HDF5
#' transparently un-shuffles the bytes before handing values back to the
#' caller, so the decoded values are always bit-for-bit identical to what
#' was written, whether or not shuffle was used. It is a storage-layer
#' optimisation, not a transformation of the data itself, and it cannot by
#' itself introduce precision loss or corruption. The whole-file MD5
#' verification performed by this function (stage 3 above) checks decoded
#' values, which empirically confirms this for every run.
#'
#' @param input.h5 Character scalar; path to the original CellBender HDF5
#'   file. This file is only ever opened for reading and is never
#'   modified.
#' @param output.h5 Character scalar; path where the repaired copy should
#'   be written. Must not already exist unless `overwrite = TRUE`.
#' @param overwrite Logical scalar; if `TRUE`, silently replace
#'   `output.h5` if it already exists. Default `FALSE`.
#' @param verify Logical scalar; if `TRUE` (the default), perform the
#'   whole-file, all-datasets MD5 verification (stage 3) after rewriting
#'   and stop with an error if anything does not match. Set to `FALSE`
#'   only for quick, uncontrolled experiments; it is strongly recommended
#'   to keep this `TRUE` for anything used in an analysis pipeline.
#'
#' @return
#' * If `/matrix/data` in `input.h5` is already stored as 32 bits: `invisible(input.h5)`.
#'   No file is written and no error is raised.
#' * Otherwise: `invisible(output.h5)`, the path to the newly written,
#'   repaired HDF5 file.
#'
#' @examples
#' \dontrun{
#' # Single file.
#' ConvertCellBenderInt64ToInt32(
#'   input.h5  = "path/to/input.h5",
#'   output.h5 = "~path/to/output.h5"
#' )
#'
#' # Looping over every file in a directory: files that are already int32
#' # are skipped automatically (no error), so no filtering is needed
#' # beforehand.
#' h5.files <- list.files(
#'   "datadir/", pattern = "\\.h5$", full.names = TRUE
#' )
#' for (f in h5.files) {
#'   ConvertCellBenderInt64ToInt32(
#'     input.h5  = f,
#'     output.h5 = sub("\\.h5$", "_int32.h5", f)
#'   )
#' }
#' }
#'
#' @seealso `.VerifyAllDatasetsIdentical()` for the standalone verification
#'   step; `.GetH5DatasetStorageSettings()` and `.H5DatasetIntSize()` for the
#'   low-level HDF5 introspection helpers used internally.
#'
#' @importFrom rhdf5 h5read h5delete h5createDataset h5write h5closeAll
#' @export
ConvertCellBenderInt64ToInt32 <- function(
  input.h5,
  output.h5,
  overwrite = FALSE,
  verify = TRUE
) {
  # --- Argument checks ----------------------------------------------------
  # Fail fast and clearly if the input does not exist, or if we would
  # silently clobber an existing output file.
  stopifnot(file.exists(input.h5))
  if (file.exists(output.h5) && !overwrite) {
    stop("`output.h5` already exists; set overwrite = TRUE to replace it.")
  }

  # --- Stage 1: cheap pre-check (datatype only, no values read) ----------
  # .H5DatasetIntSize() opens the dataset and asks HDF5 for the *storage
  # type's* byte width (e.g. 4 for int32, 8 for int64) without reading a
  # single value. This is effectively free even for huge matrices, which is
  # what lets us call this function unconditionally on every file in a
  # batch without a separate "is this file already OK?" check beforehand.
  current.size <- .H5DatasetIntSize(input.h5, "/matrix/data")
  if (current.size == 4L) {
    # Already int32 (4 bytes): nothing to do. We
    # deliberately use message(), not warning() or stop(), and we return
    # invisible(NULL) rather than raising a condition, so that calling this
    # function inside a lapply()/for loop over many files never interrupts
    # the loop for files that do not need fixing.
    message(
      "`/matrix/data` in ",
      input.h5,
      " is already 32-bit; skipping conversion."
    )
    return(invisible(input.h5))
  }

  # Ensure any HDF5 file handles opened below are released even if the
  # function exits early via stop().
  on.exit(rhdf5::h5closeAll(), add = TRUE)

  # --- Stage 2a: copy the whole file as-is --------------------------------
  # file.copy() duplicates every byte of the original file, which means
  # every group, dataset and attribute other than /matrix/data (which we
  # are about to rewrite) is preserved exactly, with no need to re-read or
  # re-write them individually.
  file.copy(input.h5, output.h5, overwrite = overwrite)

  # --- Stage 2b: learn how the ORIGINAL dataset was stored ----------------
  # We want the rewritten (int32) dataset to use the same chunk size, gzip
  # level and shuffle setting as the original (int64) dataset, rather than
  # imposing arbitrary new values. .GetH5DatasetStorageSettings() inspects
  # the dataset's HDF5 "creation property list" to recover these settings.
  original.settings <- .GetH5DatasetStorageSettings(input.h5, "/matrix/data")

  # --- Stage 2c: read the int64 values ------------------------------------
  # h5read() cannot represent a 64-bit integer as a native R integer (R
  # integers are 32-bit), so we ask it to hand the values back as doubles
  # instead. A double can exactly represent any integer up to 2^53, which
  # comfortably covers the uint32 range (up to 2^32 - 1) we are about to
  # cast into, so no precision is lost in this intermediate step.
  data.as.double <- rhdf5::h5read(
    output.h5,
    "/matrix/data",
    bit64conversion = "double"
  )
  n.values <- length(data.as.double)

  # --- Stage 2d: validate it is safe to narrow to (u)int32 ----------------
  # We only ever want to shrink the storage type, never to change the
  # values. These checks make sure that doing so cannot silently lose
  # information: no missing values, nothing negative (10X/BPCells counts
  # are non-negative), nothing above the uint32 maximum, and nothing with a
  # fractional part (which would indicate these are not raw counts).
  if (anyNA(data.as.double)) {
    stop("`/matrix/data` contains NA values; aborting conversion.")
  }
  if (min(data.as.double) < 0 || max(data.as.double) > (2**32 - 1)) {
    stop("`/matrix/data` values fall outside the uint32 range [0, 2^32 - 1].")
  }
  if (any(data.as.double != trunc(data.as.double))) {
    stop("`/matrix/data` contains non-integer values; aborting conversion.")
  }

  # Now that we know it is safe, cast to a native R integer (32-bit),
  # which is what will be written to HDF5 as a 32-bit integer type.
  data.as.integer <- as.integer(data.as.double)
  # Free the double-precision copy immediately: for large matrices this can
  # be a sizeable chunk of memory we no longer need.
  rm(data.as.double)

  # --- Stage 2e: replace the dataset on disk ------------------------------
  # h5delete() removes the existing (int64) "/matrix/data" dataset from the
  # copied file. h5createDataset() then creates a fresh dataset at the same
  # path, but as a 32-bit integer, reusing the chunk size, gzip level and
  # shuffle flag recovered from the original file in Stage 2b. h5write()
  # finally fills it with the validated, narrowed values.
  rhdf5::h5delete(output.h5, "/matrix/data")
  rhdf5::h5createDataset(
    output.h5,
    "/matrix/data",
    dims = n.values,
    storage.mode = "integer",
    chunk = min(original.settings$chunk, n.values),
    level = original.settings$gzip.level,
    shuffle = original.settings$shuffle
  )
  rhdf5::h5write(data.as.integer, output.h5, "/matrix/data")
  # Close all open HDF5 handles now so the verification step below (Stage
  # 3) is reading from a fully flushed, closed file rather than one still
  # held open by this process.
  rhdf5::h5closeAll()
  rm(data.as.integer)

  # --- Stage 3: whole-file verification -----------------------------------
  # Confirm that *every* dataset in the file - not just /matrix/data, and
  # not just the ones BPCells reads - still has exactly the same decoded
  # values as in the original file. This is the safety net that lets us
  # trust the narrowing we performed above did not alter anything.
  if (verify) {
    .VerifyAllDatasetsIdentical(input.h5, output.h5)
  }

  # --- Stage 4: return -----------------------------------------------------
  invisible(output.h5)
}
