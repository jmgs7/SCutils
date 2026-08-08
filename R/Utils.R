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

#' @title .collapse_value_columns
#'
#' @description
#' Collapse repeated non-Boolean columns after multiple full_join() operations
#'
#' This helper is the value-oriented counterpart of
#' \code{.collapse_boolean_columns()}.
#'
#' It identifies repeated column families generated by successive
#' \code{dplyr::full_join()} calls (for example, \code{metric},
#' \code{metric.x}, \code{metric.y}, \code{metric.x.y}), and collapses each
#' repeated family into a single list-column.
#'
#' In contrast to \code{.collapse_boolean_columns()}, this function does not
#' combine columns with the logical \code{&} operator. Instead, it creates a
#' vector of values for each row using \code{c()} and stores that vector in a
#' list-column.
#'
#' The function explicitly collapses only repeated families where all columns
#' are non-logical. Families that are fully logical are skipped so that Boolean
#' columns can be handled separately by \code{.collapse_boolean_columns()}.
#'
#' @param data
#' A data frame (or compatible tabular object) resulting from one or more
#' \code{dplyr::full_join()} operations.
#'
#' @return
#' The same object as \code{data}, where each repeated non-Boolean family is
#' reduced to one list-column named after the base column name. Each row of that
#' list-column contains a vector with the values coming from the repeated
#' columns in that row.
#'
#' @details
#' Base names are computed by removing one-or-more trailing \code{.x}/\code{.y}
#' suffix sequences.
#'
#' For example:
#'
#' \itemize{
#'   \item \code{metric.x}   \eqn{\rightarrow} \code{metric}
#'   \item \code{metric.y}   \eqn{\rightarrow} \code{metric}
#'   \item \code{metric.x.y} \eqn{\rightarrow} \code{metric}
#'   \item \code{metric}     \eqn{\rightarrow} \code{metric}
#' }
#'
#' Only groups with more than one column are considered repeated families.
#' Among those repeated families:
#' \itemize{
#'   \item If all columns are logical, the group is skipped.
#'   \item If columns mix logical and non-logical types, the group is skipped
#'   and a warning is emitted to avoid ambiguous coercion.
#'   \item If all columns are non-logical, each row is collapsed with
#'   \code{c()} into one vector.
#' }
#'
#' @examples
#' # Simulate two full_join()-derived value columns for the same metric.
#' value_data <- tibble::tibble(
#'   sample = c("S1", "S2"),
#'   metric.x = c(10, 20),
#'   metric.y = c(11, 21),
#'   pass.x = c(TRUE, FALSE),
#'   pass.y = c(TRUE, TRUE)
#' )
#'
#' # Collapse only non-Boolean repeated families.
#' collapsed <- value_data |>
#'   .collapse_value_columns()
#'
#' # `metric` becomes a list-column:
#' # row 1 -> c(10, 11)
#' # row 2 -> c(20, 21)
#' #
#' # Boolean family (`pass.x`, `pass.y`) is left unchanged.
#'
#' @importFrom dplyr select any_of
.collapse_value_columns <- function(data) {
  # Capture original column names once to keep grouping deterministic.
  column.names <- names(data)

  # Compute the base name of each column by removing all trailing .x/.y suffixes.
  #
  # Examples handled by this pattern:
  #   metric.x   -> metric
  #   metric.y   -> metric
  #   metric.x.y -> metric
  #   metric     -> metric
  base.names <- sub(
    pattern = "(\\.x|\\.y)+$",
    replacement = "",
    x = column.names
  )

  # Group original column names by their base name.
  #
  # Each list element contains one full family of related columns generated by
  # full_join() suffixing rules.
  grouped.names <- split(
    x = column.names,
    f = base.names
  )

  # Keep only families that truly contain repeated columns.
  #
  # Single-column families do not require any collapsing.
  repeated.groups <- grouped.names[
    lengths(grouped.names) > 1L
  ]

  # Iterate over each repeated family independently.
  for (base.name in names(repeated.groups)) {
    # Retrieve the concrete column names belonging to this family.
    repeated.names <- repeated.groups[[base.name]]

    # Extract each repeated column as a vector.
    repeated.columns <- lapply(
      X = repeated.names,
      FUN = function(column.name) data[[column.name]]
    )

    # Check whether each repeated column is logical (Boolean).
    column.is.logical <- vapply(
      X = repeated.columns,
      FUN = is.logical,
      FUN.VALUE = logical(length = 1L)
    )

    # Skip pure-Boolean families explicitly.
    #
    # This function is intentionally restricted to value columns.
    if (all(column.is.logical)) {
      next
    }

    # Skip mixed-type families (logical + non-logical) to avoid implicit
    # coercion and ambiguous semantics.
    if (any(column.is.logical) && !all(column.is.logical)) {
      warning(
        paste(
          "Skipping mixed logical/non-logical family:",
          base.name,
          "(columns:",
          paste(repeated.names, collapse = ", "),
          ")"
        )
      )
      next
    }

    # Collapse the repeated non-Boolean columns row-by-row.
    #
    # Map(c, col1, col2, ...) returns one vector per row:
    #   row_i -> c(col1[i], col2[i], ...)
    collapsed.values <- do.call(
      what = Map,
      args = c(
        f = c,
        repeated.columns
      )
    )

    # Store collapsed vectors in a single list-column with the base name.
    data[[base.name]] <- collapsed.values

    # Identify suffixed versions to remove after successful collapse.
    #
    # The unsuffixed base column (if present) is retained and overwritten above.
    suffixed.names <- repeated.names[
      grepl(
        pattern = "(\\.x|\\.y)+$",
        x = repeated.names
      )
    ]

    # Drop only auxiliary suffixed columns from the output object.
    data <- dplyr::select(
      .data = data,
      -dplyr::any_of(suffixed.names)
    )
  }

  # Return the data with non-Boolean repeated families collapsed.
  return(data)
}
