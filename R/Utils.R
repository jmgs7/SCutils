#' @title ExtractFeatureTestResults
#'
#' @description
#' Retrieves and summarises the per-cell results produced by
#' \code{.ComputeFeatureThresholdSeurat()} — which stores a computed threshold
#' value and a logical pass/fail flag in \code{meta.data} — and collapses them
#' into a tidy, batch-level \code{data.frame} that is easy to inspect, export,
#' or use for downstream QC decisions.
#'
#' The function looks for two columns in \code{meta.data} that follow the
#' naming convention established by \code{.ComputeFeatureThresholdSeurat}:
#' \itemize{
#'   \item \code{<feature>.threshold} — the numeric threshold assigned to each
#'     cell (constant within a batch).
#'   \item \code{<feature>.pass} — a logical vector indicating whether the cell
#'     passed (\code{TRUE}) or failed (\code{FALSE}) the threshold test.
#' }
#' It then groups cells by the user-supplied batch column, and for every batch
#' returns: the shared threshold value, the total number of cells, the number
#' of passing cells, and the pass rate (proportion of passing cells).
#'
#' @param SeuratObject A \code{\link[SeuratObject]{Seurat}} object whose
#'   \code{meta.data} slot contains the columns produced by
#'   \code{.ComputeFeatureThresholdSeurat()}.
#' @param batch.col A single character string specifying the name of the
#'   \code{meta.data} column that identifies the batches (e.g. sample IDs,
#'   library IDs, or any grouping variable used when calling
#'   \code{.ComputeFeatureThresholdSeurat()}).
#' @param feature A single character string giving the name of the feature
#'   (gene or metadata variable) whose test results should be extracted.
#'   Defaults to \code{"MALAT1"}, which is a common QC marker for ambient
#'   RNA contamination in single-cell experiments.
#'
#' @return A \code{data.frame} with one row per unique batch level, preserving
#'   the order in which batches first appear in \code{meta.data}. Columns are:
#'   \describe{
#'     \item{\code{<batch.col>}}{Batch identifier (character or factor), named
#'       after the input \code{batch.col} argument.}
#'     \item{\code{<feature>.threshold}}{The numeric threshold value computed
#'       for each batch by \code{.ComputeFeatureThresholdSeurat()}.}
#'     \item{\code{nCells.total}}{Integer. Total number of cells belonging to
#'       the batch.}
#'     \item{\code{nCells.pass}}{Numeric. Number of cells that passed the
#'       threshold test (\code{feature.pass == TRUE}).}
#'     \item{\code{pass.rate}}{Numeric in \code{[0, 1]}. Proportion of cells
#'       in the batch that passed the threshold test.}
#'   }
#'   Row names are set to the batch identifiers for convenient subsetting.
#'
#' @details
#' The function relies on \pkg{dplyr} for the grouping and summarisation step.
#' \code{NA} values in the pass column are silently excluded from the counts
#' and the pass-rate calculation (\code{na.rm = TRUE}); this mirrors the
#' behaviour of \code{.ComputeFeatureThresholdSeurat()} when a threshold cannot
#' be estimated for a cell.
#'
#' The row order of the returned \code{data.frame} matches the order in which
#' batches first appear in \code{meta.data}, so it is consistent with the
#' original sample ordering of the \code{SeuratObject}.
#'
#' @seealso
#' \code{.ComputeFeatureThresholdSeurat()} for the upstream function that
#' populates the \code{meta.data} columns consumed here.
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
  feature = "MALAT1"
) {
  # ---------------------------------------------------------------------------
  # 1. Build the expected meta.data column names from the feature name.
  #    .ComputeFeatureThresholdSeurat() always stores results under
  #    "<feature>.threshold" and "<feature>.pass".
  # ---------------------------------------------------------------------------

  # Column that holds the numeric threshold assigned per batch.
  threshold.col <- paste0(feature, ".threshold")

  # Column that holds the logical pass/fail flag per cell.
  pass.col <- paste0(feature, ".pass")

  # Cache the column names of meta.data to avoid repeated @-slot access.
  meta.data.cols <- colnames(SeuratObject@meta.data)

  # ---------------------------------------------------------------------------
  # 2. Input validation — fail early with an informative error if any of the
  #    three required columns is missing from meta.data.
  # ---------------------------------------------------------------------------

  # Check that the threshold column exists; stop if not found.
  if (!threshold.col %in% meta.data.cols) {
    stop(threshold.col, " column not found in meta.data")

    # Check that the pass/fail column exists; stop if not found.
  } else if (!pass.col %in% meta.data.cols) {
    stop(pass.col, " column not found in meta.data")

    # Check that the user-supplied batch column exists; stop if not found.
  } else if (!batch.col %in% meta.data.cols) {
    stop(batch.col, " column not found in meta.data")
  }

  # ---------------------------------------------------------------------------
  # 3. Build a long (cell-level) data.frame with only the three columns needed.
  #    Using [[...]] on meta.data preserves the vector class (numeric/logical)
  #    and avoids creating an unnecessary copy of the entire meta.data slot.
  # ---------------------------------------------------------------------------

  feature.test <- data.frame(
    # One row per cell: the batch identifier for that cell.
    batch = SeuratObject@meta.data[[batch.col]],
    # The threshold value assigned to the cell's batch.
    feature.threshold = SeuratObject@meta.data[[threshold.col]],
    # Whether the cell passed the threshold test.
    feature.pass = SeuratObject@meta.data[[pass.col]]
  ) |>

    # ---------------------------------------------------------------------------
    # 4. Collapse from cell level to batch level via dplyr.
    # ---------------------------------------------------------------------------

    # Group all rows that belong to the same batch together.
    dplyr::group_by(batch) |>

    dplyr::summarise(
      # The threshold is constant within a batch — take the first non-NA value.
      feature.threshold = dplyr::first(feature.threshold),

      # Count total cells in the batch (including those with NA in pass col).
      nCells.total = dplyr::n(),

      # Count cells that explicitly passed (TRUE); NAs are excluded by na.rm.
      nCells.pass = sum(feature.pass, na.rm = TRUE),

      # Compute the fraction of passing cells; NAs excluded from numerator
      # and denominator via na.rm = TRUE.
      pass.rate = mean(feature.pass, na.rm = TRUE),

      # Drop the grouping structure from the result to return a plain data.frame.
      .groups = "drop"
    ) |>

    # Convert the tibble returned by dplyr::summarise() to a base data.frame
    # for compatibility with downstream code that may not expect a tibble.
    as.data.frame()

  # ---------------------------------------------------------------------------
  # 5. Restore the original batch order.
  #    dplyr::summarise() sorts batches alphabetically; here we re-index the
  #    result using the order of first appearance in meta.data, which matches
  #    the sample ordering of the SeuratObject.
  # ---------------------------------------------------------------------------

  # Set batch names as row names to enable label-based subsetting in step below.
  row.names(feature.test) <- feature.test$batch

  # Reorder rows to match the first-appearance order of batches in meta.data.
  feature.test <- feature.test[unique(SeuratObject@meta.data[[batch.col]]), ]

  # ---------------------------------------------------------------------------
  # 6. Rename columns to match the original meta.data naming conventions so
  #    that the output is self-documenting and directly traceable back to the
  #    Seurat object.
  # ---------------------------------------------------------------------------

  # Replace generic interim column names with the original meta.data names.
  names(feature.test) <- c(
    batch.col, # e.g. "orig.ident"
    threshold.col, # e.g. "MALAT1.threshold"
    "nCells.total",
    "nCells.pass",
    "pass.rate"
  )

  # ---------------------------------------------------------------------------
  # 7. Return the tidy, batch-level summary data.frame.
  # ---------------------------------------------------------------------------
  return(feature.test)
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
#' The function assumes that the columns being combined contain logical
#' values, such as TRUE and FALSE.
#'
#' If at least one NA value is found in any column, the function stops
#' execution using stop(). It does not automatically convert NA values to
#' FALSE, because the appropriate strategy for handling missing values should
#' be decided explicitly before combining the Boolean results.
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
#' The columns are then grouped according to their base name, and the values
#' within each group are combined using purrr::reduce() and the logical AND
#' operator.
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

    # Reduce the family to one logical vector using pairwise AND.
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
  data
}

#' @title FeatureTest
#'
#' @description
#' A helper function to compute the threshold and pass/fail results for a given features.
#' FeatureTest work by fetching the feature data from the Seurat object,
#' computing the threshold based on the specified method (MAD,percentile, or user's threshold),
#' and determining whether each cell passes or fails the threshold test.
#' The results are then added to the Seurat object's metadata.
#' For each feature, the user must imput a threshold and a "operator". Operators determine the
#' comparison that it will be performed between the feature values and the threshold.
#' If multiple features are given, but only one threshold and operator are provided,
#' the same threshold and operator will be applied to all features.
#'
#' @param SeuratObject A Seurat object containing the data.
#' @param features The features (gene or metadata variable) to compute the threshold for.
#'   You can provide the same feature multiple times with different thresholds and operators
#'   to compute multiple tests for the same feature. FeatureTest will automatically handle
#'   the repeated features and compute the results for each test separately, and the final
#'   output will be TRUE if all tests for the same feature pass, and FALSE if any test fails.
#' @param assay The assay to use for fetching the feature data.
#' @param layers The layer type to use for fetching the feature data (counts, data...).
#'   If NULL, the default layer (meta.data) is used. (NAs, emtpy strings, and "null" are treated as NULL).
#'   In the case of a split SeuratObject, the function will automatically compute the threshold for each layer separately
#'   respecting the data structure.
#' @param thresholds A list of thresholds to use for computing the pass/fail results.
#'   To use Median Absolute Deviation (MAD) or percentile method, use mad or percentile respectively.
#' @param operators A list of operators to use for computing the pass/fail results.
#'   Valid operators are "upper", "lower", or "both" in the case of MAD or percentile method.
#'   For numeric thresholds, input the desired operators as strings. Valid operators: ">", "<", ">=", "<=", "==", "!=".
#'   Use params nmad and percentile to set the number of MADs or the quantile to consider as outliers.
#' @param nmad The number of MADs to use to compute the threshold.
#' @param percentile The quantile to consider as outliers. It must be an integer between 1 and 99.
#'   Default is 1 and computes the top and bottom 1% of the distribution.
#'
#' @return SeuratObject A Seurat object with the computed threshold and pass/fail results added to the meta.data slot.
#'
#' @examples
#'\dontrun{
#'
#' library(Seurat)
#' library(dplyr)
#' library(data.table)
#' library(purrr)
#'
#' seurat_object <- FeatureTest(
#'   SeuratObject = seurat_object,
#'   features = c("Gene1", "Gene2", "percent.mt"),
#'   assay = "RNA",
#'   layers = c("couts", "data", NA),
#'   thresholds = "mad",
#'   operators = "both",
#'   nmad = 3,
#'   percentile = 1
#' )
#' }
#'
#' @import Seurat SeuratObject
#' @importFrom stats median mad quantile
#' @importFrom dplyr full_join
#' @importFrom purrr reduce
#' @importFrom data.table rbindlist()
#' @export

FeatureTest <- function(
  SeuratObject,
  features,
  assay = "RNA",
  layers = NULL,
  thresholds = "mad",
  operators = "both",
  nmad = 3,
  percentile = 1
) {
  # ---------------------------------------------------------------------------
  # 1. INPUT VALIDATION AND NORMALIZATION.
  # ---------------------------------------------------------------------------

  # Validate Seurat object: required for FetchData semantics and structure.
  if (!inherits(SeuratObject, "Seurat")) {
    stop("'SeuratObject' must be a Seurat object.")
  }

  # Validate assay: must be a character string and exist in the Seurat object.
  if (
    !is.character(assay) ||
      length(assay) != 1 ||
      !(assay %in% names(SeuratObject@assays))
  ) {
    stop(
      "'assay' must be a single character string corresponding to an existing assay in the Seurat object."
    )
  }

  # ---------------------------------------------------------------------------
  # LAYER VALIDATION
  # ---------------------------------------------------------------------------
  # Validate layers: must be NULL or a character vector of valid layer names in the Seurat object.
  # Normalize layer input to one layer assignment per requested feature.
  # Metadata variables will ignore the layer later, but keeping the vector
  # aligned to features simplifies the rest of the code.
  if (is.null(layers)) {
    layer.per.feature <- rep(list(NULL), length(features))
  } else if (is.character(layers)) {
    if (length(layers) == 1L) {
      layer.per.feature <- rep(list(layers), length(features))
    } else if (length(layers) == length(features)) {
      layer.per.feature <- lapply(layers, function(x) {
        if (is.na(x) || tolower(x) == "null" || tolower(x) == "") {
          return(NULL)
        }
        return(x)
      })
    } else {
      stop(sprintf(
        "'layers' has length %d but 'features' has length %d; it must be length 1 or length(features).",
        length(layers),
        length(features)
      ))
    }
  } else {
    stop("'layer' must be NULL or a character vector.")
  }
  names(layer.per.feature) <- features
  # After normalization, layer.per.feature is a list of the same length as features, where each element is either NULL or a character string specifying the layer to use for that feature.
  # Now we check if the specified layers exist in the Seurat object, if not NULL.
  if (!is.null(layers)) {
    # Check that each specified layer exists in the Seurat object.
    object.layers <- SeuratObject::Layers(SeuratObject, assay = assay)
    lapply(layer.per.feature, function(layer) {
      if (!is.null(layer)) {
        if (!any(grepl(layer, object.layers))) {
          stop(
            "Specified layer '",
            layer,
            "' is not present in the Seurat object for assay '",
            assay,
            "'. Available layers: ",
            paste(object.layers, collapse = ", ")
          )
        }
      }
    })
  }

  # Validate features: must be a character vector or character string.
  if (!is.character(features) || length(features) == 0) {
    stop("'features' must be a non-empty character vector.")
  }

  # ---------------------------------------------------------------------------
  # THRESHOLD VALIDATION
  # ---------------------------------------------------------------------------
  # Validate thresholds.
  # We will recycle the normalizeVline framework from other functions of the package.

  # DEFINE SINGLE ENTRY NORMALIZATION FUNCTION
  normalizeThresholdEntry <- function(
    x,
    valid.keywords = c("mad", "percentile")
  ) {
    # Handle NULL, NA, or empty string inputs by outputting an error message.
    if (
      is.null(x) ||
        is.na(x) ||
        (is.character(x) && tolower(x) %in% c("null", "na", ""))
    ) {
      stop(
        "Invalid threshold value. Threshold values can be either a supported keyword ('mad', 'percentile') or a numeric value."
      )
    }
    # Character inputs can be either a supported keyword or a numeric value
    # encoded as text.
    if (is.character(x)) {
      x.lower <- tolower(x)

      if (x.lower %in% valid.keywords) {
        return(x.lower)
      }

      x.numeric <- suppressWarnings(as.numeric(x))
      if (!is.na(x.numeric)) {
        return(x.numeric)
      }

      stop(sprintf(
        "threshold entry '%s' is not a valid keyword or numeric value.",
        x
      ))
    }
    # Numeric inputs must be scalar and non-missing.
    if (is.numeric(x)) {
      if (length(x) == 1L && !is.na(x)) {
        return(x)
      }

      stop("numeric thresholds must be a single non-NA value.")
    }
    # Any other type is unsupported.
    stop(sprintf(
      "threshold entry has unexpected type: %s",
      paste(class(x), collapse = ", ")
    ))
  }

  # DEFINE FULL VECTOR NORMALIZATION FUNCTION
  # Normalize a full threshold input vector against an expected target length.
  # The target length represents the number of requested features.
  normalizeThresholdVector <- function(
    thresholds,
    target.length = length(features),
    valid.keywords = c("mad", "percentile"),
    context = "features"
  ) {
    # A NULL threshold is not valid, and will output an error.
    if (is.null(thresholds)) {
      stop(
        "Thresholds cannot be NULL. Please provide a valid threshold value or keyword ('mad', 'percentile')."
      )
    }
    # Numeric inputs are accepted as either one value recycled to all targets
    # or one value per target.
    if (is.numeric(thresholds)) {
      if (length(thresholds) == 1L) {
        entry <- normalizeThresholdEntry(thresholds, valid.keywords)
        return(rep(list(entry), target.length))
      }
      if (length(thresholds) == target.length) {
        return(lapply(
          as.list(thresholds),
          normalizeThresholdEntry,
          valid.keywords = valid.keywords
        ))
      }
      stop(sprintf(
        "'thresholds' numeric input has length %d but %s has length %d; it must be length 1 or %d.",
        length(thresholds),
        context,
        target.length,
        target.length
      ))
    }

    # Character inputs follow the same contract as numeric inputs, but each
    # element may be either a keyword or a numeric value encoded as a string.
    if (is.character(thresholds)) {
      if (length(thresholds) == 1L) {
        entry <- normalizeThresholdEntry(thresholds, valid.keywords)
        return(rep(list(entry), target.length))
      }

      if (length(thresholds) == target.length) {
        return(lapply(
          thresholds,
          normalizeThresholdEntry,
          valid.keywords = valid.keywords
        ))
      }

      stop(sprintf(
        "'thresholds' character input has length %d but %s has length %d; it must be length 1 or %d.",
        length(thresholds),
        context,
        target.length,
        target.length
      ))
    }

    stop("'thresholds' must be numeric, or character.")
  }

  # Apply the normalization function to the input thresholds vector, ensuring it matches the expected length of features and contains valid keywords or numeric values.
  thresholds.per.feature <- normalizeThresholdVector(
    thresholds,
    target.length = length(features),
    valid.keywords = c("mad", "percentile"),
    context = "features"
  )
  names(thresholds.per.feature) <- features

  # ---------------------------------------------------------------------------
  # OPERATOR VALIDATION
  # ---------------------------------------------------------------------------
  # Validate operators parameter: must be one of the valid operators.
  if (is.character(operators) && length(operators) > 0) {
    if (
      !all(
        tolower(operators) %in%
          c("upper", "lower", "both", ">", "<", ">=", "<=", "==", "!=")
      )
    ) {
      stop(
        "'operators' must be one of the valid operators: upper, lower, both, >, <, >=, <=, ==, !=."
      )
    } else if (length(operators) == 1L) {
      operators <- rep(operators, length(features))
      operators <- tolower(operators)
    } else if (length(operators) == length(features)) {
      operators <- tolower(operators)
    } else {
      stop(sprintf(
        "'operators' has length %d but 'features' has length %d; it must be length 1 or length(features).",
        length(operators),
        length(features)
      ))
    }
  } else {
    stop("'operators' must be a character vector.")
  }
  operators.per.feature <- as.list(operators)
  names(operators.per.feature) <- features

  # Validate nmad parameter: must be a positive numeric value.
  if (!is.numeric(nmad) || length(nmad) != 1 || nmad <= 0) {
    stop("'nmad' must be a single positive numeric value.")
  }

  # Validate percentile parameter: must be a numeric value between 1 and 99.
  if (
    !is.numeric(percentile) ||
      length(percentile) != 1 ||
      percentile < 1 ||
      percentile > 99
  ) {
    stop("'percentile' must be a single numeric value between 1 and 99.")
  }

  # ---------------------------------------------------------------------------
  # 2. DEFINE FUNCTION TO CALCULATE THRESHOLDS BASED ON MAD OR PERCENTILE
  # ---------------------------------------------------------------------------
  calculateStats <- function(data, stat, nmad = 3, percentile = 1) {
    if (stat == "mad") {
      median <- stats::median(data, na.rm = TRUE)
      mad <- stats::mad(data, na.rm = TRUE)
      lower.threshold <- median - nmad * mad
      upper.threshold <- median + nmad * mad
      return(list(lower = lower.threshold, upper = upper.threshold))
    } else if (stat == "percentile") {
      percentile <- percentile / 100
      lower.threshold <- quantile(data, probs = percentile, na.rm = TRUE)
      upper.threshold <- quantile(
        data,
        probs = 1 - percentile,
        na.rm = TRUE
      )
      return(list(lower = lower.threshold, upper = upper.threshold))
    }
  }

  # ---------------------------------------------------------------------------
  # 3. COMPUTE THRESHOLDS AND PASS/FAIL RESULTS FOR EACH FEATURE
  # ---------------------------------------------------------------------------
  # Iterate by integer position rather than by feature name.
  #
  # Using seq_along(features) gives each entry a unique index `i` even when
  # the same feature name appears more than once in `features`. This makes
  # positional lookup in the per-feature tables (layer, threshold, operator)
  # unambiguous: `list[[i]]` always retrieves the i-th element, whereas
  # `list[[name]]` returns only the first element matching that name.
  new.metadata <- lapply(seq_along(features), function(i) {
    # Retrieve the feature name for this position.
    feature <- features[[i]]

    # ---------------------------------------------------------------------------
    # 3.1. Determine the layer to use for the current feature and fetch the
    #      corresponding data from the Seurat object.
    # ---------------------------------------------------------------------------

    # Look up the layer assigned to this specific entry by position, not by name.
    # This ensures that when the same feature appears twice with different layers,
    # each occurrence retrieves its own intended layer.
    layer <- layer.per.feature[[i]]

    # To respect the split layer structure, we will fetch the cells id of each layer and use them to subset the Seurat object for the feature threshold computation.
    if (!is.null(layer)) {
      # If a specific layer is provided, use it to obtain the split layers.
      split.layers <- SeuratObject::Layers(
        SeuratObject,
        assay = assay,
        search = layer
      )
    } else {
      # If no specific layer is provided, it will default to the meta.data (layer = NULL). To respect the split layer structure, we will use the counts layers to obtain the split layers.
      # If the object is not split, it will return a single layer containing all the cells of the SeuratObject anyway.
      split.layers <- SeuratObject::Layers(
        SeuratObject,
        assay = assay,
        search = "counts"
      )
    }

    # Generate a list with the cells id of each layer to subset the Seurat object for the feature threshold computation.
    # In the case the layer is NULL (use meta.data), we ensure to obtain the layer structure by using the counts layers.
    cells.list <- lapply(split.layers, function(split.layer) {
      cells <- SeuratObject::Cells(SeuratObject, layer = split.layer)
      return(cells)
    })
    names(cells.list) <- split.layers

    # Look up threshold and operator for this entry by position.
    #
    # Positional lookup mirrors the layer lookup above: each duplicate occurrence
    # of the same feature name retrieves its own independently validated setting.
    threshold <- thresholds.per.feature[[i]]

    # Get the operator for the current feature from the normalized operators vector.
    operator <- operators.per.feature[[i]]

    # ---------------------------------------------------------------------------
    # 3.2. Compute the threshold and pass/fail results for each split layer
    #      and combine them into a single data frame.
    # ---------------------------------------------------------------------------
    feature.metadata <- lapply(cells.list, function(split.cells) {
      # Fetch the feature data for the current split layer.
      feature.data <- SeuratObject::FetchData(
        SeuratObject,
        assay = assay,
        vars = feature,
        layer = layer,
        cells = split.cells
      )

      # ---------------------------------------------------------------------------
      # A. If the threshold is a character keyword (mad or percentile),
      #    compute the threshold based on the specified method and the provided
      #    nmad or percentile values.
      # ---------------------------------------------------------------------------
      if (is.character(threshold) && threshold %in% c("mad", "percentile")) {
        # If the threshold is a character keyword (mad or percentile), validate the operator and output an error message if invalid.
        if (!tolower(operator) %in% c("upper", "lower", "both")) {
          stop(paste(
            "Invalid operator",
            operator,
            "for feature",
            feature,
            ". Valid operators for MAD or percentile thresholds are 'upper', 'lower', or 'both'."
          ))
        }

        # Calculate the lower and upper thresholds based on the specified method (MAD or percentile) and the provided nmad or percentile values.
        stats <- calculateStats(
          feature.data[[feature]],
          threshold,
          nmad,
          percentile
        )

        # Determine whether each cell passes or fails the threshold test based on the specified operator and the computed thresholds.
        #
        # A default branch is included to stop explicitly if an unrecognised operator
        # somehow reaches this point. Without a default, switch() returns NULL silently,
        # which would propagate to results.df and produce an all-NA pass column with
        # no diagnostic message.
        feature.pass <- switch(
          operator,
          "upper" = feature.data[[feature]] < stats$upper,
          "lower" = feature.data[[feature]] > stats$lower,
          "both" = feature.data[[feature]] > stats$lower &
            feature.data[[feature]] < stats$upper,
          # Default: should never be reached because operator is validated above.
          # Included defensively so any bypass produces an informative error
          # rather than a silent NULL.
          stop(sprintf(
            "Unrecognised operator '%s' for feature '%s'. Valid operators for MAD/percentile thresholds are 'upper', 'lower', 'both'.",
            operator,
            feature
          ))
        )

        # ---------------------------------------------------------------------------
        # B. If the threshold is a character keyword (mad or percentile),
        #    compute the threshold based on the specified method and the provided
        #    nmad or percentile values.
        # ---------------------------------------------------------------------------
      } else {
        feature.pass <- switch(
          operator,
          ">" = feature.data[[feature]] > threshold,
          "<" = feature.data[[feature]] < threshold,
          ">=" = feature.data[[feature]] >= threshold,
          "<=" = feature.data[[feature]] <= threshold,
          "==" = feature.data[[feature]] == threshold,
          "!=" = feature.data[[feature]] != threshold,
          stop(paste(
            "Invalid operator",
            operator,
            "for feature",
            feature,
            ". Valid operators for numeric thresholds are '>', '<', '>=', '<=', '==', '!='."
          ))
        )
      }

      # ---------------------------------------------------------------------------
      # 3.3. Create a data frame with the cell IDs and the pass/fail results for the current
      # feature and split layer. The column name for the pass/fail results is renamed
      # to include the feature name for clarity.
      # ---------------------------------------------------------------------------
      results.df <- data.frame(
        cell.id = row.names(feature.data),
        feature.pass = feature.pass
      )
      # Rename the column of the results data frame to include the feature name for clarity.
      names(results.df)[2] <- paste0(feature, ".pass")
      return(results.df)
    }) |> # (END OF LAYER LAPPLY)
      # Combine the results from all layers into a single data frame.
      data.table::rbindlist()

    # ---------------------------------------------------------------------------
    # 3.4. Join all the results for each feature into a single data frame, ensuring
    # that the cell IDs are preserved and aligned with the Seurat object's metadata.
    # ---------------------------------------------------------------------------
    return(feature.metadata)
  }) |> # (END OF FEATURE LAPPLY)
    # Join the results of all features into a single data frame.
    purrr::reduce(dplyr::full_join, by = "cell.id") |>
    # Collapse multiple tests for the same feature into a single column using logical AND.
    .collapse_boolean_columns() |>
    # Convert the combined results to a data frame.
    as.data.frame()
  # Set the row names of the feature.metadata data frame to the cell IDs for proper alignment with the Seurat object's metadata.
  row.names(new.metadata) <- new.metadata$cell.id
  # Delete the cell.id column from feature.metadata as it is now redundant with the row names
  new.metadata$cell.id <- NULL

  # ---------------------------------------------------------------------------
  # 4. FINAL OUTPUT GENERATION.
  # ---------------------------------------------------------------------------
  SeuratObject <- Seurat::AddMetaData(
    object = SeuratObject,
    metadata = new.metadata
  )
  return(SeuratObject)
} # (END OF FUNCTION).
