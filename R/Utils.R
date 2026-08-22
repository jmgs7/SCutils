#' @title ExtractFeatureTestResults
#'
#' @description
#' Retrieves and summarises per-cell feature-test results stored in
#' \code{meta.data} and collapses them into a batch-level \code{data.frame}
#' that is easy to inspect, export, or use for downstream QC decisions.
#'
#' The function expects at least the \code{<feature>.pass} column, and can
#' optionally extract \code{<feature>.threshold} when available.
#'
#' @param SeuratObject A \code{\link[SeuratObject]{Seurat}} object whose
#'   \code{meta.data} slot contains the columns produced by feature-testing
#'   workflows (for example \code{.CalculateFeatureThresholdSeurat()}).
#' @param batch.col A single character string specifying the name of the
#'   \code{meta.data} column that identifies batches (e.g. sample IDs,
#'   library IDs, or any grouping variable used in the upstream workflow).
#' @param feature A single character string giving the name of the feature
#'   (gene or metadata variable) whose test results should be extracted.
#'   Defaults to \code{"MALAT1"}.
#' @param extract.threshold Logical. If \code{TRUE}, include
#'   \code{<feature>.threshold} in the output if that column exists.
#'   If the threshold column is missing, the function warns and continues
#'   with pass/fail summary only.
#'
#' @return A \code{data.frame} with one row per unique batch level, preserving
#'   the order in which batches first appear in \code{meta.data}. Columns are:
#'   \describe{
#'     \item{\code{<batch.col>}}{Batch identifier (character or factor), named
#'       after the input \code{batch.col} argument.}
#'     \item{\code{nCells.total}}{Integer. Total number of cells belonging to
#'       the batch.}
#'     \item{\code{nCells.pass}}{Numeric. Number of cells that passed the
#'       threshold test (\code{feature.pass == TRUE}).}
#'     \item{\code{pass.rate}}{Numeric in \code{[0, 1]} when at least one
#'       non-missing pass value exists in the batch; \code{NA} otherwise.}
#'     \item{\code{<feature>.threshold}}{The threshold value for the batch,
#'       returned only when \code{extract.threshold = TRUE} and the threshold
#'       column exists in \code{meta.data}.}
#'   }
#'   Row names are set to the batch identifiers for convenient subsetting.
#'
#' @details
#' \code{NA} values in \code{<feature>.pass} are excluded from \code{nCells.pass}
#' and \code{pass.rate}. If all pass values are \code{NA} in a batch,
#' \code{pass.rate} is returned as \code{NA}.
#'
#' The row order of the returned \code{data.frame} matches the order in which
#' batches first appear in \code{meta.data}.
#'
#' @seealso
#' \code{.CalculateFeatureThresholdSeurat()} for an upstream function that can
#' populate \code{<feature>.threshold} and \code{<feature>.pass} in
#' \code{meta.data}.
#'
#' @importFrom dplyr group_by summarise first n
#'
#' @examples
#' \dontrun{
#' # Assuming `seu` already has MALAT1.threshold and MALAT1.pass in meta.data:
#' results <- ExtractFeatureTestResults(
#'   SeuratObject = seu,
#'   batch.col    = "orig.ident",
#'   feature      = "MALAT1"
#' )
#' print(results)
#' #   orig.ident  MALAT1.threshold nCells.total nCells.pass pass.rate
#' # s1         s1             0.02         3000        2850      0.95
#' # s2         s2             0.03         2500        2200      0.88
#' }
#'
#' @export
ExtractFeatureTestResults <- function(
  SeuratObject,
  batch.col,
  feature = "MALAT1",
  extract.threshold = FALSE
) {
  # ---------------------------------------------------------------------------
  # 1. Validate user input and build expected meta.data column names.
  # ---------------------------------------------------------------------------

  # Basic argument checks to fail early with informative messages.
  if (!is.character(batch.col) || length(batch.col) != 1) {
    stop("batch.col must be a single character string")
  }
  if (!is.character(feature) || length(feature) != 1) {
    stop("feature must be a single character string")
  }
  if (!is.logical(extract.threshold) || length(extract.threshold) != 1) {
    stop("extract.threshold must be a single logical value")
  }

  # Build the expected metadata column names from the feature name.
  threshold.col <- paste0(feature, ".threshold")
  pass.col <- paste0(feature, ".pass")

  # Cache the meta.data column names to avoid repeated slot access.
  meta.data.cols <- colnames(SeuratObject@meta.data)

  # ---------------------------------------------------------------------------
  # 2. Validate required metadata columns.
  # ---------------------------------------------------------------------------

  # The pass/fail result column must always exist.
  if (!(pass.col %in% meta.data.cols)) {
    stop(pass.col, " column not found in meta.data")
  }

  # The user-provided batch column must always exist.
  if (!(batch.col %in% meta.data.cols)) {
    stop(batch.col, " column not found in meta.data")
  }

  # Threshold extraction is optional; if missing, warn and continue.
  if (extract.threshold && !(threshold.col %in% meta.data.cols)) {
    warning(
      threshold.col,
      " column not found in meta.data, skipping extraction of threshold values"
    )
    extract.threshold <- FALSE
  }

  # ---------------------------------------------------------------------------
  # 3. Build a cell-level table and collapse to batch-level summary.
  # ---------------------------------------------------------------------------

  feature.test <- data.frame(
    # One row per cell: the batch identifier for that cell.
    batch = SeuratObject@meta.data[[batch.col]],
    # Whether the cell passed the threshold test.
    feature.pass = SeuratObject@meta.data[[pass.col]]
  ) |>

    # Group all rows that belong to the same batch together.
    dplyr::group_by(batch) |>

    dplyr::summarise(
      # Count total cells in the batch (including those with NA in pass col).
      nCells.total = dplyr::n(),

      # Count cells that explicitly passed (TRUE); NAs are excluded by na.rm.
      nCells.pass = sum(feature.pass, na.rm = TRUE),

      # Compute pass rate on non-missing pass values; if all are missing,
      # return NA instead of NaN for downstream safety.
      pass.rate = if (all(is.na(feature.pass))) {
        NA_real_
      } else {
        mean(feature.pass, na.rm = TRUE)
      },

      # Drop grouping structure and return a plain summary table.
      .groups = "drop"
    ) |>
    as.data.frame()

  # ---------------------------------------------------------------------------
  # 4. Optionally extract batch-level threshold values.
  # ---------------------------------------------------------------------------

  if (extract.threshold) {
    feature.threshold <- SeuratObject@meta.data |>
      dplyr::select(
        batch = dplyr::all_of(batch.col),
        threshold = dplyr::all_of(threshold.col)
      ) |>
      dplyr::group_by(batch) |>
      dplyr::summarise(
        n.threshold = dplyr::n_distinct(threshold, na.rm = TRUE),
        threshold = if (all(is.na(threshold))) {
          NA_real_
        } else {
          dplyr::first(threshold[!is.na(threshold)])
        },
        .groups = "drop"
      ) |>
      as.data.frame()

    # Warn if more than one distinct threshold is found within a batch.
    inconsistent.batch <- feature.threshold$batch[
      feature.threshold$n.threshold > 1
    ]
    if (length(inconsistent.batch) > 0) {
      warning(
        "Multiple non-missing threshold values detected within batch(es): ",
        paste(inconsistent.batch, collapse = ", "),
        ". Using first non-missing threshold per batch."
      )
    }

    feature.threshold$n.threshold <- NULL

    # Attach threshold column while preserving all summary rows.
    feature.test <- dplyr::left_join(
      x = feature.test,
      y = feature.threshold,
      by = "batch"
    )
  }

  # ---------------------------------------------------------------------------
  # 5. Restore original batch order and rename output columns.
  # ---------------------------------------------------------------------------

  # Reorder rows to match first appearance of batch values in meta.data.
  batch.order <- unique(SeuratObject@meta.data[[batch.col]])
  feature.test <- feature.test[
    match(batch.order, feature.test$batch),
    ,
    drop = FALSE
  ]

  # Replace generic interim column names with user-facing names.
  if (extract.threshold) {
    names(feature.test) <- c(
      batch.col,
      "nCells.total",
      "nCells.pass",
      "pass.rate",
      threshold.col
    )

    # Reorder columns to keep threshold next to batch identifier.
    feature.test <- feature.test[, c(
      batch.col,
      threshold.col,
      "nCells.total",
      "nCells.pass",
      "pass.rate"
    )]
  } else {
    names(feature.test) <- c(
      batch.col,
      "nCells.total",
      "nCells.pass",
      "pass.rate"
    )
  }

  # Set stable row names for convenient label-based subsetting.
  row.names(feature.test) <- make.unique(ifelse(
    is.na(feature.test[[batch.col]]),
    "NA",
    as.character(feature.test[[batch.col]])
  ))

  # ---------------------------------------------------------------------------
  # 6. Return tidy batch-level summary.
  # ---------------------------------------------------------------------------
  return(as.data.frame(feature.test))
}


#' @title .SetCommonScales
#'
#' @description
#' Compute common x and y axis limits across all subplots in a patchwork object
#' and reapply those limits so that every panel is displayed on the same scale.
#' This is useful when comparing multiple ggplot2 panels that were assembled
#' with patchwork and need a consistent visual frame after plotting.
#'
#' @details
#' The function inspects each subplot inside the patchwork object with
#' `ggplot2::ggplot_build()` and extracts the computed panel ranges from the
#' first panel. It then calculates global minimum and maximum values across all
#' subplots and applies them back to the assembled patchwork using
#' `ggplot2::scale_x_continuous()` and `ggplot2::scale_y_continuous()`.
#'
#' This approach is appropriate for standard single-panel ggplot objects.
#' If a subplot uses faceting, only the first panel range is used.
#'
#' @param plot A patchwork object containing two or more ggplot2 subplots.
#'
#' @return A patchwork object with common x and y axes limits applied across all
#' subplots.
#'
#' @section Notes:
#' \itemize{
#'   \item The function assumes continuous x and y axes.
#'   \item For faceted plots, only the first panel is used to extract limits.
#'   \item Using scale limits will constrain the displayed range; it may not be
#'   equivalent to a visual zoom.
#' }
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(patchwork)
#'
#' p1 <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
#' p2 <- ggplot(mtcars, aes(wt, qsec)) + geom_point()
#' p3 <- ggplot(mtcars, aes(wt, disp)) + geom_point()
#'
#' pw <- wrap_plots(p1, p2, p3)
#' .SetCommonScales(pw)
#' }
#'
#' @importFrom ggplot2 ggplot_build coord_cartesian
#' @import patchwork
#'
#' @noRd

.SetCommonScales <- function(plot, collect.axes = FALSE) {
  # 1. Map over all sub-plots inside the patchwork structure
  all.limits.list <- lapply(seq_along(plot), function(plot.index) {
    # Build the individual sub-plot
    built <- ggplot2::ggplot_build(plot[[plot.index]])

    # Grab calculated panel ranges
    # (Handles standard plots. If faceting is used, [[1]] gets the first panel)
    x.range <- built$layout$panel_params[[1]]$x.range
    y.range <- built$layout$panel_params[[1]]$y.range

    # Return structured as a single row data frame
    data.frame(
      plot_index = plot.index,
      x_min = x.range[1],
      x_max = x.range[2],
      y_min = y.range[1],
      y_max = y.range[2]
    )
  })

  # 2. Combine rows into one comprehensive table
  plot.limits <- do.call(rbind, all.limits.list)

  # 3. Find absolute X and Y boundaries across the entire patchwork grid
  global.x.min <- min(plot.limits$x_min)
  global.x.max <- max(plot.limits$x_max)

  global.y.min <- min(plot.limits$y_min)
  global.y.max <- max(plot.limits$y_max)

  # 4. Redraw the plot with the new limits.
  plot <- (plot) &
    ggplot2::coord_cartesian(
      xlim = c(global.x.min, global.x.max),
      ylim = c(global.y.min, global.y.max)
    )

  return(plot)
}

#' @title .collapse_boolean_columns
#'
#' @description
#' Combine repeated Boolean columns after multiple full_join() operations
#'
#' This function is designed to process the result of a pipeline in which
#' multiple data frames have been combined using full_join().
#'
#' The input data frames may contain different test results obtained from
#' applying several thresholds to the same sample.
#'
#' For example, after several full_join() operations, the same metric may
#' appear under the following names:
#'
#'   percent.mt.pass.x
#'   percent.mt.pass.y
#'   percent.mt.pass
#'
#' The function recognises that all these columns belong to the same family
#' and combines them using the logical AND operator:
#'
#'   percent.mt.x & percent.mt.y & percent.mt
#'
#' The result is stored in a single column called percent.mt.
#'
#' The function can process several independent column families at the same
#' time, for example:
#'
#'   percent.mt.pass.x       percent.mt.pass.y       percent.mt.pass
#'   percent.ribo.pass.x     percent.ribo.pass.y     percent.ribo.pass
#'   nFeature_RNA.pass.x     nFeature_RNA.pass.y
#'
#' Each family is combined independently.
#'
#' The function does not need to know:
#'
#'   * The names of the metrics.
#'   * The number of input data frames.
#'   * The number of thresholds that were applied.
#'   * The column that was used as the full_join() key.
#'
#' The function collapses only repeated families where all columns are logical
#' values, such as TRUE and FALSE. Repeated families that contain non-logical
#' columns are skipped.
#'
#' If at least one NA value is found in any column, the function emits a
#' warning. It does not automatically convert NA values to FALSE, because the
#' appropriate strategy for handling missing values should be decided explicitly
#' before combining the Boolean results.
#'
#' @param data
#' A data frame resulting from one or more full_join() operations.
#'
#' @return
#' The same type of object as data, with each repeated column family reduced
#' to a single column per metric.
#'
#' @details
#' dplyr::full_join() normally uses the suffixes ".x" and ".y" when it finds
#' columns with the same name that are not used as join keys.
#'
#' When several full_join() operations are performed successively, the same
#' metric may end up with names such as:
#'
#'   metric.x
#'   metric.y
#'   metric
#'
#' or even:
#'
#'   metric.x.x
#'   metric.x.y
#'   metric.y
#'
#' The function removes all final ".x" and ".y" suffix sequences in order to
#' obtain the base name of each column.
#'
#' For example:
#'
#'   metric.x       -> metric
#'   metric.y       -> metric
#'   metric.x.x     -> metric
#'   metric.x.y     -> metric
#'   metric         -> metric
#'
#' The columns are then grouped according to their base name. For each
#' repeated family where all columns are logical, values are combined using
#' \code{purrr::reduce()} and the logical AND operator. Families that are
#' not fully logical are left unchanged.
#'
#' @examples
#' # Create a table containing the results of a first threshold.
#' first_test <- tibble::tibble(
#'   sample = c("S1", "S2", "S3"),
#'   percent.mt = c(TRUE, TRUE, FALSE)
#' )
#'
#' # Create a table containing the results of a second threshold.
#' second_test <- tibble::tibble(
#'   sample = c("S1", "S2", "S3"),
#'   percent.mt = c(TRUE, FALSE, TRUE)
#' )
#'
#' # Create a table containing a different metric.
#' ribosomal_test <- tibble::tibble(
#'   sample = c("S1", "S2", "S3"),
#'   percent.ribo = c(TRUE, TRUE, FALSE)
#' )
#'
#' # Store the tables in a list.
#' data.frames <- list(
#'   first_test,
#'   second_test,
#'   ribosomal_test
#' )
#'
#' # Combine all tables using full_join().
#' joined_data <- data.frames |>
#'   purrr::reduce(
#'     .f = dplyr::full_join,
#'     by = "sample"
#'   )
#'
#' # Reduce repeated columns using logical AND.
#' result <- joined_data |>
#'   .collapse_boolean_columns()
#'
#' @importFrom dplyr select any_of
#' @importFrom purrr reduce
#'
#' @noRd

.collapse_boolean_columns <- function(data) {
  # Capture original column names once so grouping logic is deterministic.
  column.names <- names(data)

  # Identify columns containing at least one NA value.
  #
  # The function keeps existing behavior: it warns but does not stop.
  na.columns <- column.names[
    vapply(
      X = data,
      FUN = anyNA,
      FUN.VALUE = logical(1)
    )
  ]

  # Emit a single warning that lists all columns containing NA values.
  if (length(na.columns) > 0L) {
    warning(
      paste(
        "NA values were found in the following columns:",
        paste(na.columns, collapse = ", ")
      )
    )
  }

  # Compute each column's base name by removing one-or-more trailing .x/.y suffixes.
  #
  # Examples:
  #   metric.x     -> metric
  #   metric.y     -> metric
  #   metric.x.y   -> metric
  #   metric       -> metric
  base.names <- sub(
    pattern = "(\\.x|\\.y)+$",
    replacement = "",
    x = column.names
  )

  # Group original column names by their base name.
  #
  # This gives a named list where each element contains all versions of the same
  # logical metric (e.g., metric, metric.x, metric.y).
  grouped.names <- split(
    x = column.names,
    f = base.names
  )

  # Keep only groups that truly have repeated columns.
  #
  # Single-name groups do not need reduction and are left untouched.
  repeated.groups <- grouped.names[
    lengths(grouped.names) > 1L
  ]

  # Iterate over each repeated base-name family.
  for (base.name in names(repeated.groups)) {
    # Retrieve the concrete column names for this family.
    repeated.names <- repeated.groups[[base.name]]

    # Extract each column as a vector using [[ ]].
    #
    # This avoids `[.data.table` join semantics and is stable for both
    # data.frame and data.table inputs.
    repeated.columns <- lapply(
      X = repeated.names,
      FUN = function(column.name) data[[column.name]]
    )

    # Determine whether each repeated column is logical.
    #
    # This function must collapse only Boolean families.
    column.is.logical <- vapply(
      X = repeated.columns,
      FUN = is.logical,
      FUN.VALUE = logical(1)
    )

    # Skip this family if any repeated column is not logical.
    #
    # This guarantees we apply '&' only to explicit Boolean vectors.
    if (!all(column.is.logical)) {
      next
    }

    # Reduce the logical family to one vector using pairwise AND.
    #
    # This preserves previous semantics exactly:
    #   metric = metric.x & metric.y (& metric if present, etc.)
    data[[base.name]] <- purrr::reduce(
      .x = repeated.columns,
      .f = `&`
    )

    # Identify only suffixed versions that should be dropped after reduction.
    #
    # If an unsuffixed base column exists, it is retained (now overwritten with
    # the reduced values). If it does not exist, it has just been created above.
    suffixed.names <- repeated.names[
      grepl(
        pattern = "(\\.x|\\.y)+$",
        x = repeated.names
      )
    ]

    # Remove auxiliary suffixed columns from the result object.
    data <- dplyr::select(
      .data = data,
      -dplyr::any_of(suffixed.names)
    )
  }

  # Return the collapsed table with one boolean column per base metric.
  return(data)
}


#' @title FilterVariableFeatures
#'
#' @description
#' Filters the variable features of a Seurat object to remove TCR/BCR and IG genes,
#' mitochondrial and reibosomal genes, and non-coding and antisense genes.
#'
#' @param SeuratObject A Seurat object containing variable features.
#' @param pattern A character vector of regular expression patterns to filter variable features.
#'   If NULL, default patterns will be used to remove TCR/BCR and IG genes,
#'   HLA genes, mitochondrial and ribosomal genes, and non-coding and antisense genes.
#' @param assay The assay to filter variable features from. Defaults to "RNA".
#' @param verbose Logical. If TRUE, prints the number of features removed and remaining.
#' @return A Seurat object with filtered variable features.
#'
#' @import Seurat
#' @import SeuratObject
#'
#' @export
FilterVariableFeatures <- function(
  SeuratObject,
  pattern = NULL,
  assay = "RNA",
  verbose = TRUE
) {
  # VALIDATE INPUTS
  if (!inherits(SeuratObject, "Seurat")) {
    stop("SeuratObject must be a Seurat object.")
  }
  if (!is.character(pattern) && !is.null(pattern)) {
    stop(
      "Pattern must be a character vector or NULL."
    )
  }
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("Verbose must be a single logical value.")
  }
  if (!assay %in% names(SeuratObject@assays)) {
    stop("Assay '", assay, "' not found in Seurat object.")
  }

  # Get the variable features
  hvg <- Seurat::VariableFeatures(SeuratObject)

  if (is.null(pattern)) {
    # Otherwise, use the default patterns to remove
    # Define patterns for genes to remove
    pattern <- c(
      "^TR[ABDG][VDJ]", # TCR genes
      "^IG[HKL][VDJ]", # BCR/IG genes
      "^HLA-", # HLA genes (avoids multi-donor issues)
      "^MT-", # Mitochondrial genes
      "^RP[LS]", # Ribosomal protein genes
      "^MRP[LS]", # Mitochondrial ribosomal protein genes
      "^LINC", # Long intergenic non-coding RNAs
      "^LOC", # Eliminate LOC genes (uncharacterized)
      "\\.[0-9]+$", # Versioned genes
      "-DT$", # Divergent transcripts
      "-AS[0-9]*$" # Antisense transcripts
    )
  }

  # Filter the variable features based on the patterns
  for (regexp in pattern) {
    hvg <- grep(regexp, hvg, invert = TRUE, value = TRUE)
  }

  # Update the Seurat object with the filtered variable features
  Seurat::VariableFeatures(SeuratObject) <- hvg

  if (verbose) {
    message("Filtered variable features. Remaining features: ", length(hvg))
  }

  return(SeuratObject)
}

#' @title SelectPCs
#'
#' @description Select the number of principal components to use for downstream analysis.
#' The function uses the intrinsicDimension::maxLikGlobalDimEst to estimate the intrinsic
#' dimensionality of the PCA embeddings in the Seurat object. It handles large datasets by chunking
#' the data to avoid memory issues and overflow errors. If the intrinsic dimension cannot be
#' estimated, it falls back to a pointwise estimation method and takes the median of the valid
#' estimates as a robust global estimate.
#'
#' @param SeuratObject A Seurat object.
#'
#' @return Numerical A vector of selected principal components.
#'
#' @import Seurat
#' @import SeuratObject
#' @importFrom intrinsicDimension maxLikGlobalDimEst maxLikPointwiseDimEst
#' @export
#'
SelectPCs <- function(SeuratObject, seed = 42L) {
  # VALIDATE INPUTS
  if (!inherits(SeuratObject, "Seurat")) {
    stop("SeuratObject must be a Seurat object.")
  }
  if (!is.integer(seed) || seed < 0 || seed > .Machine$integer.max) {
    stop(
      "Seed must be a non-negative integer within the range of valid integers. Max value is ",
      .Machine$integer.max,
      "."
    )
  }

  # Set the random seed for reproducibility
  set.seed(seed)

  # Obtain the PCA embeddings from the Seurat object
  PCs <- SeuratObject@reductions$pca@cell.embeddings
  # Get the number of cells
  n.PCs <- nrow(PCs)

  # For large datasets, chunk the data to avoid memory issues and overflow errors.
  if (n.PCs >= 100000) {
    # Split the data into 10 chunks for processing
    chunks <- split(sample(1:n.PCs), rep(1:10, length.out = n.PCs))
    # Estimate the intrinsic dimensions for each chunk and store the results
    dim.estimates <- sapply(chunks, function(idx) {
      intrinsicDimension::maxLikGlobalDimEst(
        PCs[idx, ],
        k = 20,
        unbiased = TRUE,
        neighborhood.aggregation = "robust"
      )
    })
    # Print the average and maximum intrinsic dimensions across all chunks
    message(paste0("Average dimensions: ", mean(unlist(dim.estimates))))
    message(paste0("Max dimensions: ", max(unlist(dim.estimates))))
    # We consider the maximum intrinsic dimension across all chunks as the final estimate
    # We round up to the nearest integer.
    est.PC <- ceiling(max(unlist(dim.estimates)))

    # For smaller datasets, we can directly estimate the intrinsic dimensions without chunking.
  } else if (n.PCs < 100000) {
    int.dim <- intrinsicDimension::maxLikGlobalDimEst(
      SeuratObject@reductions$pca@cell.embeddings,
      k = 20,
      unbiased = TRUE,
      neighborhood.aggregation = 'robust'
    ) # this is teh best but crashes when run on the full dataset
    est.PC <- ceiling(int.dim[[1]])
  }

  # If it is still NA, fallback to the pointwise estimation method, and take the median of the valid estimates as a robust global estimate.
  if (exists("est.PC") == FALSE) {
    warning("Using alternative approach to estimate intrinsic dimensions")
    X <- SeuratObject@reductions$pca@cell.embeddings
    pt <- intrinsicDimension::maxLikPointwiseDimEst(X, k = 20, unbiased = TRUE)
    m <- pt$dim.est
    m2 <- m[is.finite(m) & m > 0] # keep only finite positive estimates
    int.dim <- median(m2) # robust global estimate = median of valid ones
    est.PC <- ceiling(int.dim[[1]])
  }

  message(paste0("Intrinsic dimensions: ", est.PC))
  return(1:est.PC)
}
