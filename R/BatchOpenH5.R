#' Batch-open multiple .h5ad files and convert features to gene symbols
#'
#' For each file in `files.set` this function:
#' - constructs the full path using `files.dir`,
#' - if a directory with suffix `"_BP"` doesn't exist it reads the HDF5 matrix with
#'   `open_matrix_10x_hdf5()` and writes the matrix directory using `write_matrix_dir()`,
#' - loads the BP matrix directory with `open_matrix_dir()` and converts ENSEMBL IDs
#'   to gene symbols via `Azimuth:::ConvertEnsembleToSymbol()`.
#'
#' Processing is parallelized with `parallel::mclapply()` on Unix-like systems.
#' On Windows, it falls back to `lapply()`.
#'
#' @param files Character vector. Filepaths (including the .h5* extension) to process.
#' @param data.dir Character. Directory where the BP matrices will be saved
#'   (without the `"_BP"` suffix). If `NULL`, uses the first file's path.
#' @param mc.cores Integer. Number of cores for `mclapply()` (default `4`).
#'
#' @return A named list of loaded/converted matrices (one per file).
#'
#' @examples
#' \dontrun{
#' BatchOpenH5(
#'   files = c("/path/to/h5ads/sample1.h5ad", "/path/to/h5ads/sample2.h5ad"),
#'   data.dir = "/path/to/data",
#'   mc.cores = 4
#' )
#' }
#'
#' @importFrom parallel mclapply
#' @importFrom Azimuth ConvertEnsembleToSymbol
#' @import BPCells
#' @export

BatchOpenH5 <- function(files, BP.data.dir = NULL, mc.cores = length(files)) {
  if (is.null(BP.data.dir)) {
    BP.data.dir <- dirname(files[1])
  }

  process_one_file <- function(file) {
    open.path <- file
    save.path <- file.path(BP.data.dir, paste0(gsub(".h5*", "", file), "_BP"))

    if (!dir.exists(save.path)) {
      data <- open_matrix_10x_hdf5(open.path)
      write_matrix_dir(mat = data, dir = save.path)
    }

    mat <- open_matrix_dir(dir = save.path)
    Azimuth:::ConvertEnsembleToSymbol(mat = mat, species = "human")
  }

  if (.Platform$OS.type == "windows") {
    data.list <- lapply(files.set, process_one_file)
  } else {
    data.list <- parallel::mclapply(
      files,
      process_one_file,
      mc.cores = mc.cores
    )
  }

  names(data.list) <- gsub(".h5*", "", basename(files))
  return(data.list)
}
