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
      # Rename result columns to include the feature name for clarity.
      names(results.df)[2] <- paste0(feature, ".pass")

      return(results.df)
    }) |> # (END OF LAYER LAPPLY)
      # Combine the results from all layers into a single data frame.
      #
      # `fill = TRUE` is a defensive guard in case an upstream branch ever
      # introduces optional columns. With the stable result schema above, this
      # should not change normal outputs.
      data.table::rbindlist(fill = TRUE)

    # ---------------------------------------------------------------------------
    # 3.4. Join all the results for each feature into a single data frame, ensuring
    # that the cell IDs are preserved and aligned with the Seurat object's metadata.
    # ---------------------------------------------------------------------------
    return(feature.metadata)
  }) |> # (END OF FEATURE LAPPLY)
    # Join the results of all features into a single data frame.
    purrr::reduce(dplyr::full_join, by = "cell.id")

  new.metadata <- new.metadata |>
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
