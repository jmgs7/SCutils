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
#' @importFrom ggplot2 ggplot_build scale_x_continuous scale_y_continuous
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
  # Supress the "already present scale" warning.
  plot <- (plot) &
    suppressWarnings(ggplot2::scale_x_continuous(
      limits = c(global.x.min, global.x.max)
    )) &
    suppressWarnings(ggplot2::scale_y_continuous(
      limits = c(global.y.min, global.y.max)
    ))

  return(plot)
}
