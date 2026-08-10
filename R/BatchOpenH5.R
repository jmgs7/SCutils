#' Batch-open HDF5 single-cell matrices into BPCells format
#'
#' For each file in `files`, this function:
#' - creates a BPCells matrix directory with `"_BP"` suffix in `BP.data.dir`
#'   (if it does not already exist),
#' - reads the matrix from either 10X HDF5 (`platform = "10X"`) or AnnData HDF5
#'   (`platform = "anndata"`),
#' - optionally replaces 10X feature IDs with `/matrix/features/name` when
#'   `use.names = TRUE`,
#' - reopens the saved BPCells matrix directory,
#' - optionally rewrites `mat@dir` to a portable relative path when
#'   `relative = TRUE`, and
#' - optionally converts ENSEMBL IDs to symbols with `ConvertEnsembleToSymbol2()`
#'   when `ensembl.to.symbol = TRUE` and `use.names = FALSE`.
#'
#' Processing uses is parallelize using the future() framework. Future options
#' must be set in the global environment before calling this function. For example, to use 4 cores:
#'
#' \dontrun{
#' library(future)
#' plan(multisession, workers = 4)
#' }
#'
#' @param files Character vector of input `.h5` file paths.
#' @param relative Logical; if `TRUE`, store matrix directory paths as
#'   `./<basename(BP.data.dir)>/<matrix_dir>` in each matrix object's `@dir` slot.
#'   This is useful for portability if the BP directory is available from the
#'   working directory. Default is `TRUE`.
#' @param BP.data.dir Character scalar directory where `*_BP` matrix folders are
#'   created. If `NULL`, uses `dirname(files[[1]])`.
#' @param platform Character scalar, either `"10X"` or `"anndata"`. Default is
#'   `"10X"`.
#' @param use.names Logical; only used for `platform = "10X"`. If `TRUE`, replaces
#'   row names with `/matrix/features/name` from the HDF5 file. Default is `TRUE`.
#' @param ensembl.to.symbol Logical; if `TRUE` and `use.names = FALSE`, converts
#'   ENSEMBL IDs to gene symbols with `ConvertEnsembleToSymbol2()`. Default is
#'   `FALSE`.
#' @param species Character scalar species passed to `ConvertEnsembleToSymbol2()`.
#'   Default is `"human"`.
#' @param generate.metadata Logical; if `TRUE`, also returns a per-cell metadata
#'   table with `cell.tag` and `sample.procedence`. Default is `FALSE`.
#'
#' @return If `generate.metadata = FALSE`, a named list of BPCells matrices. Names
#'   are derived from input basenames without `.h5*` suffix.
#'   If `generate.metadata = TRUE`, a list with:
#'   - `data.list`: the named matrix list
#'   - `metadata`: a `data.frame` with per-cell sample provenance.
#'
#' @examples
#' \dontrun{
#' mats <- BatchOpenH5(
#'   files = c("/path/to/sample1.h5", "/path/to/sample2.h5"),
#'   BP.data.dir = "/path/to/bp_matrices",
#'   platform = "10X",
#'   mc.cores = 2
#' )
#'
#' out <- BatchOpenH5(
#'   files = c("/path/to/sample1.h5"),
#'   generate.metadata = TRUE
#' )
#' }
#'
#' @importFrom rhdf5 h5read
#' @import BPCells
#' @import futurize
#' @export

BatchOpenH5 <- function(
  files,
  relative = TRUE,
  BP.data.dir = NULL,
  platform = "10X",
  use.names = TRUE,
  ensembl.to.symbol = FALSE,
  species = "human",
  generate.metadata = FALSE
) {
  # Use the directory of the first file if BP.data.dir is not provided
  if (is.null(BP.data.dir)) {
    BP.data.dir <- dirname(files[[1]])
  }

  # Internal function to process a single file which will be parallelized with mclapply
  process_one_file <- function(
    file,
    BP.data.dir,
    relative,
    platform,
    use.names,
    ensembl.to.symbol,
    species
  ) {
    open.path <- file
    save.path <- file.path(
      BP.data.dir,
      paste0(gsub(".h5*", "", basename(file)), "_BP")
    ) # Create save path with _BP suffix

    if (!dir.exists(save.path)) {
      if (platform == '10X') {
        message(paste("Processing 10X file:", open.path))
        data <- open_matrix_10x_hdf5(open.path)
        if (use.names) {
          feature.names <- rhdf5::h5read(open.path, "/matrix/features/name")
          message(paste("Changed gene names for 10X file:", open.path))
          rownames(data) <- feature.names
        }
      } else if (platform == 'anndata') {
        message(paste("Processing anndata file:", open.path))
        data <- open_matrix_anndata_hdf5(open.path)
      } else {
        stop("Unsupported platform. Please specify '10X' or 'anndata'.")
      }
      write_matrix_dir(mat = data, dir = save.path)
    }

    mat <- open_matrix_dir(dir = save.path)
    message(paste("Loaded matrix file:", open.path))

    if (relative) {
      relative.path <- file.path(
        ".",
        basename(BP.data.dir),
        basename(save.path)
      )
      mat@dir <- relative.path
      message(paste0(
        "Set relative matrix directory to: ",
        mat@dir,
        " for file: ",
        open.path
      ))
    }

    if (ensembl.to.symbol && !use.names) {
      message(paste(
        "Converting ENSEMBL IDs to gene symbols for file:",
        save.path
      ))
      mat <- ConvertEnsembleToSymbol2(
        mat = mat,
        species = species,
        mirror = "useast"
      )
      message(paste(
        "Converted ENSEMBL IDs to gene symbols for file:",
        save.path
      ))
    }
    return(mat)
  }

  # Multicore processing of files with mclapply. On Windows, this will run sequentially due to mc.cores being set to 1.
  data.list <- lapply(
    files,
    process_one_file,
    BP.data.dir = BP.data.dir,
    relative = relative,
    platform = platform,
    use.names = use.names,
    ensembl.to.symbol = ensembl.to.symbol,
    species = species
  ) |>
    futurize::futurize()
  # Set the names of the list to the base names of the files without the .h5* extension
  names(data.list) <- gsub(".h5*", "", basename(files))

  # Generate simple metadata indicating the procedence of the sample if requested.
  if (generate.metadata) {
    meta.list <- lapply(
      names(data.list),
      function(sample) {
        cells <- colnames(data.list[[sample]])
        data.frame(
          cell.tag = cells,
          sample.procedence = sample,
          row.names = cells,
          stringsAsFactors = FALSE
        )
      }
    ) |>
      futurize::futurize()

    metadata <- Reduce(rbind, meta.list)
    output <- list(data.list = data.list, metadata = metadata)
    return(output)
  }

  return(data.list)
}
