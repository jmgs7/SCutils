#' Batch-open multiple .h5ad files, generates the BP matrices and simple metadata, and converts features to gene symbols
#'
#' For each file in `files.set` this function:
#' - constructs the full path using `files.dir`,
#' - if a directory with suffix `"_BP"` doesn't exist it reads the HDF5 matrix with
#'   `open_matrix_10x/anndata_hdf5()` and writes the matrix directory using `write_matrix_dir()`,
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
#' @param relative Character. If `TRUE`, the paths in the matrix metadata will be relative to `BP.data.dir` (`./basename(BP.data.dir)/matrix_folder`, so the BP matrices folder is portable. Default is `TRUE`. NOTE: The resulting folder must be in the working directory for Seurat to find it.
#' @param platform Character. Platform of the input files, either `"10X"` or `"anndata"`. Default is `"10X"`.
#' @param use.names Logical. Whether to use the gene symbols rather the gene ids. WARNING: Only with 10X matrices. The gene names should be under the default /matrix/features/name directory of an 10X h5 file. Default is `TRUE`.
#' @param ensembl.to.symbol Logical. Whether to convert ENSEMBL IDs to gene symbols. Won't work if use.names is `TRUE`. Default is `FALSE`.
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
#' @importFrom rhdf5 h5read
#' @import BPCells
#' @export

BatchOpenH5 <- function(
  files,
  relative = TRUE,
  BP.data.dir = NULL,
  platform = "10X",
  use.names = TRUE,
  ensembl.to.symbol = FALSE,
  species = "human",
  generate.metadata = FALSE,
  mc.cores = length(files)
) {
  # Use the directory of the first file if BP.data.dir is not provided
  if (is.null(BP.data.dir)) {
    BP.data.dir <- dirname(files[[1]])
  }

  # Windows does not allow parallel processing with mclapply, so we set mc.cores to 1
  if (.Platform$OS.type == "windows") {
    mc.cores <- 1
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
  data.list <- parallel::mclapply(
    files,
    process_one_file,
    BP.data.dir = BP.data.dir,
    relative = relative,
    platform = platform,
    use.names = use.names,
    ensembl.to.symbol = ensembl.to.symbol,
    species = species,
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
