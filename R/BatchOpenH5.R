#' Batch-open multiple .h5ad files, generates the BP matrices and simple metadata, and converts features to gene symbols
#'
#' For each file in `files.set` this function:
#' - constructs the full path using `files.dir`,
#' - if a directory with suffix `"_BP"` doesn't exist it reads the HDF5 matrix with
#'   `open_matrix_10x_hdf5()` and writes the matrix directory using `write_matrix_dir()`,
#' - loads the BP matrix directory with `open_matrix_dir()` and converts ENSEMBL IDs
#'   to gene symbols via `Azimuth:::ConvertEnsembleToSymbol()` if `ensembl.to.symbol` is `TRUE`, and
#' - generates a basic metadata table with the procedence of the sample if `generate.metadata` is `TRUE`.
#'
#' Processing is parallelized with `parallel::mclapply()` on Unix-like systems.
#' On Windows, it falls back to `lapply()`.
#'
#' @param files Character vector. Filepaths (including the .h5* extension) to process.
#' @param BP.data.dir Character. Directory where the BP matrices will be saved
#'   (without the `"_BP"` suffix). If `NULL`, uses the first file's path.
#' @param platform Character. Platform of the input files, either `"10X"` or `"anndata"`. Default is `"10X"`.
#' @param ensembl.to.symbol Logical. Whether to convert ENSEMBL IDs to gene symbols.
#' @param species Character. Species for gene symbol conversion (default `"human"`).
#' @param generate.metadata Logical. Whether to generate a basic metadata table with the procedence of the sample (default `FALSE`).
#' @param mc.cores Integer. Number of cores for `mclapply()` (default is the number of files).
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
#' @import BPCells
#' @export

BatchOpenH5 <- function(
  files,
  BP.data.dir = NULL,
  platform = "10X",
  ensembl.to.symbol = TRUE,
  species = "human",
  generate.metadata = FALSE,
  mc.cores = length(files)
) {
  # Use the directory of the first file if BP.data.dir is not provided
  if (is.null(BP.data.dir)) {
    BP.data.dir <- dirname(files[1])
  }

  # Windows does not allow parallel processing with mclapply, so we set mc.cores to 1
  if (.Platform$OS.type == "windows") {
    mc.cores <- 1
  }

  # Internal function to process a single file which will be parallelized with mclapply
  process_one_file <- function(file, platform) {
    open.path <- file
    save.path <- file.path(
      BP.data.dir,
      paste0(gsub(".h5*", "", basename(file)), "_BP")
    ) # Create save path with _BP suffix

    if (!dir.exists(save.path)) {
      if (platform == '10X') {
        data <- open_matrix_10x_hdf5(open.path)
        write_matrix_dir(mat = data, dir = save.path)
      } else if (platform == 'anndata') {
        data <- open_matrix_anndata_hdf5(open.path)
        write_matrix_dir(mat = data, dir = save.path)
      } else {
        stop("Unsupported platform. Please specify '10X' or 'anndata'.")
      }
    }

    mat <- open_matrix_dir(dir = save.path)
    if (ensembl.to.symbol) {
      mat <- ConvertEnsembleToSymbol2(
        # TODO: Integrate this functionality in this module for better performance, and remove redundant database downloading process.
        mat = mat,
        species = species,
        mirror = "useast"
      )
    }
    return(mat)
  }

  # Multicore processing of files with mclapply. On Windows, this will run sequentially due to mc.cores being set to 1.
  data.list <- parallel::mclapply(
    files,
    process_one_file,
    platform = platform,
    mc.cores = mc.cores
  )
  # Set the names of the list to the base names of the files without the .h5* extension
  names(data.list) <- gsub(".h5*", "", basename(files))

  # Generate simple metadata indicating the procedence of the sample if requested.
  if (generate.metadata) {
    meta.list <- parallel::mclapply(
      names(data.list),
      function(sample) {
        cells <- colnames(data.list[[sample]])
        data.frame(
          cell.tag = cells,
          sample.procedence = sample,
          row.names = cells,
          stringsAsFactors = FALSE
        )
      },
      mc.cores = mc.cores
    )

    metadata <- Reduce(rbind, meta.list)
    output <- list(data.list = data.list, metadata = metadata)
    return(output)
  }

  return(data.list)
}
