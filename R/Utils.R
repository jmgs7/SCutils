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

#' @title .FeatureTest
#'
#'  @description
#' A helper function to compute the threshold and pass/fail results for a given features.
#'
#' @param SeuratObject A Seurat object containing the data.
#' @param assay The assay to use for fetching the feature data.
#' @param layers The layer type to use for fetching the feature data (meta.data, counts, data...).
#'   If NULL, the default layer (meta.data) is used.
#'   In the case of a split SeuratObject, the function will automatically compute the threshold for each layer separately
#'   respecting the data structure.
#' @param features The features (gene or metadata variable) to compute the threshold for.
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
#'
#'

.FeatureTest <- function(
  SeuratObject,
  assay = "RNA",
  layers = NULL,
  split.metadata = NULL,
  features,
  thresholds = "mad",
  operators = "both",
  stat = "mad",
  nmad = 3,
  percentile = 1
) {
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
      layer.per.feature <- lapply(layer, function(x) {
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

  # Validate split.metadata: must be NULL or a character string corresponding to a column in meta.data.
  if (!is.null(split.metadata)) {
    if (!is.character(split.metadata) || length(split.metadata) != 1) {
      stop(
        "'split.metadata' must be a single character string corresponding to a column in meta.data."
      )
    }
    if (!split.metadata %in% colnames(SeuratObject@meta.data)) {
      stop(
        "'split.metadata' column '",
        split.metadata,
        "' not found in meta.data."
      )
    }
  }

  # Validate features: must be a character vector or character string.
  if (!is.character(features) || length(features) == 0) {
    stop("'features' must be a non-empty character vector.")
  }

  # Validate thresholds.
  # We will recycle the normalizeVline framework from other functions of the package.
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
        "Invalid threshold value. Threshold values can be either a supported keyword ('upper', 'lower', 'both') or a numeric value."
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

  # Normalize a full threshold input vector against an expected target length.
  # The target length represents either the number of requested features or the
  # number of group levels, depending on the plotting mode.
  normalizeThresholdVector <- function(
    thresholds,
    target.length = length(features),
    valid.keywords = c("mad", "percentile"),
    context = "features"
  ) {
    # A NULL threshold is not valid, and will output an error.
    if (is.null(thresholds)) {
      stop(
        "Thresholds cannot be NULL. Please provide a valid threshold value or keyword ('upper', 'lower', 'both')."
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

  # Validate operators parameter: must be one of the valid operators.
  if (is.character(operators) && length(operators) > 0) {
    if (
      !tolower(operators) %in%
        c("upper", "lower", "both", ">", "<", ">=", "<=", "==", "!=")
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
  names(operators) <- features

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

  calculateStats <- function(data, stat, nmad = 3, percentile = 1) {
    if (stat == "mad") {
      median <- stats::median(data, na.rm = TRUE)
      mad <- stats::mad(data, na.rm = TRUE)
      lower.threshold <- median - nmad * mad
      upper.threshold <- median + nmad * mad
      return(list(lower = lower.threshold, upper = upper.threshold))
    } else if (stat == "percentile") {
      percentile <- percentile / 100
      lower.threshold <- quantile(data, probs = percentile / 100, na.rm = TRUE)
      upper.threshold <- quantile(
        data,
        probs = 1 - percentile / 100,
        na.rm = TRUE
      )
      return(list(lower = lower.threshold, upper = upper.threshold))
    }
  }

  lapply(features, function(feature) {
    # For each feature, determine the layer to use (if any) and fetch the corresponding data from the Seurat object.
    layer <- layer.per.feature[[feature]]

    # To respect the split layer structure, we will fetch the data for each layer separately and compute the threshold for each layer independently.
    if (!is.null(layer)) {
      # If a specific layer is provided, use it to fetch the data.
      split.layers <- SeuratObject::Layers(
        SeuratObject,
        assay = assay,
        search = layer
      )
    } else {
      # If no specific layer is provided, it will default to the meta.data. To respect the split layer structure, we will fetch the data for each layer separately and compute the threshold for each layer independently using the split.metatadata column to identify the layers.
      split.layers <- layer
      if (!is.null(split.metadata)) {
        metadata.layers <- unique(SeuratObject@meta.data[[split.metadata]])
      } else {
        metadata.layers <- NULL
      }
    }

    # Get the threshold value for the current feature from the normalized thresholds list.
    threshold <- thresholds.per.feature[[feature]]

    # Get the operator for the current feature from the normalized operators vector.
    operator <- operators[[feature]]

    if (!is.null(layer)) {
      feature.metadata <- lapply(split.layers, function(split.layer) {
        # Fetch the feature data for the current split layer.
        feature.data <- SeuratObject::FetchData(
          SeuratObject,
          assay = assay,
          vars = feature,
          layer = split.layer
        )

        # Compute the nmad or percentile threshold if required.
        if (is.character(threshold) && threshold %in% c("mad", "percentile")) {
          if (!tolower(operator) %in% c("upper", "lower", "both")) {
            stop(paste(
              "Invalid operator",
              operator,
              "for feature",
              feature,
              ". Valid operators for MAD or percentile thresholds are 'upper', 'lower', or 'both'."
            ))
          }

          stats <- calculateStats(
            feature.data[[feature]],
            threshold,
            nmad,
            percentile
          )

          feature.pass <- switch(
            operator,
            "upper" = feature.data[[feature]] < stats$upper,
            "lower" = feature.data[[feature]] > stats$lower,
            "both" = feature.data[[feature]] > stats$lower &
              feature.data[[feature]] < stats$upper
          )
        } else {
          # If a numeric threshold is provided, use it directly.
          threshold.value <- threshold

          feature.pass <- switch(
            operator,
            ">" = feature.data[[feature]] > threshold.value,
            "<" = feature.data[[feature]] < threshold.value,
            ">=" = feature.data[[feature]] >= threshold.value,
            "<=" = feature.data[[feature]] <= threshold.value,
            "==" = feature.data[[feature]] == threshold.value,
            "!=" = feature.data[[feature]] != threshold.value,
            stop(paste(
              "Invalid operator",
              operator,
              "for feature",
              feature,
              ". Valid operators for numeric thresholds are '>', '<', '>=', '<=', '==', '!='."
            ))
          )
        }

        # Create a data.frame to store the threshold value and pass/fail results for each cell in the Seurat object.
        results.df <- data.frame(
          cell.id = row.names(feature.data),
          feature.pass = feature.pass
        )
        # Rename the column of the results data frame to include the feature name for clarity.
        names(results.df)[2] <- paste0(feature, ".pass")

        return(results.df)
      }) |> # Combine the results from all layers into a single data frame.
        data.table::rbindlist() |>
        as.data.frame()

      # Set the row names of the feature.metadata data frame to the cell IDs
      # for proper alignment with the Seurat object's metadata.
      row.names(feature.metadata) <- feature.metadata$cell.id
      # Delete the cell.id column from feature.metadata as it is now redundant with the row names.
      feature.metadata$cell.id <- NULL
    } else {
      feature.data <- SeuratObject::FetchData(
        SeuratObject,
        assay = assay,
        layer = layer,
        vars = c(feature, split.metadata)
      )

      if (!is.null(metadata.layers)) {}
    }

    # END OF FEATRUE LAPPLY
  })

  # END OF FUNCTION.
}

## TODO: Adapt the function to the new philosophy, it only matters if the layer is null or not. split.layer will either take as values the list of layers of the splited object, or a single layer if it is not splitted. If we are using the metadata layer, split.layer will be NULL, and automatically the function will use the metadata layer. In that case, we split the metadata with the split.metadata column. It that is aslo null, it will compute the threshold (if necesary) for the entire SeuratObject. If the layer is not null, it will compute the threshold for each layer separately, respecting the split structure of the SeuratObject.
