#' @title FeatureDensityPlot
#' @description
#' `FeatureDensityPlot()` draws density plots for features from a Seurat object
#' without splitting the object. It supports metadata columns, assay features
#' (from any layer), and dimensional reduction variables via `Seurat::FetchData()`.
#'
#' Grouping is resolved from the specified `group.by` variable, avoiding expensive
#' object partitioning for large datasets.
#'
#' The function supports:
#' - overlayed grouped densities,
#' - split (faceted-by-group) density panels,
#' - optional vertical reference lines (`vline`) in dashed red,
#' - optional independent median overlays (`plot.median`) in dashed black,
#' - custom plot titles via `plot.title`,
#' - multi-feature output as a named list (one plot per feature),
#' - per-feature layer selection via `layer`.
#'
#' The implementation uses sequential base `lapply()` for efficiency.
#'
#' @param SeuratObject A Seurat object.
#' @param features Character vector of features to plot (metadata columns, assay
#'   feature names, or reduction variable names like `"PC_1"`).
#' @param group.by Character scalar. Grouping variable: a metadata column name,
#'   `"ident"` (default, SeuratObject$active.ident), or `NULL` for no grouping.
#' @param split.plot Logical. If `TRUE` (default), creates one panel per group level
#'   for each feature. If `FALSE`, groups are overlaid in one panel per feature.
#' @param scale.colors Character scalar. Viridis palette option used for grouped
#'   density colors (`"viridis"`, `"magma"`, `"plasma"`, `"inferno"`,
#'   `"cividis"`, `"rocket"`, `"mako"`, `"turbo"`). Default is `"viridis"`.
#' @param ncol Integer or `NULL`. Number of columns for split panels within each
#'   feature plot when `split.plot = TRUE`. If `NULL`, inferred automatically.
#' @param vline Optional reference-line specification (drawn in dashed red).
#'   Accepted values: `NULL`, `"mean"`, `"median"`, `"upper"`, `"lower"`,
#'   `"both"`, a numeric value, or a numeric/character vector.
#'   For feature-level vectors, length must be 1 or `length(features)`, and
#'   each element may be a keyword or a numeric value (including numeric values
#'   encoded as strings).
#'   When `length(features) == 1`, `group.by` is not `NULL`, and
#'   `split.plot = TRUE`, `vline` may alternatively have length equal to the
#'   number of group levels, in which case each split panel receives its own
#'   specification.
#'   `"upper"`, `"lower"`, and `"both"` use median +/- `nmad` * MAD.
#'   Length-1 vectors are recycled to all features or groups.
#' @param layer Character, `NULL` or `NA`. Specifies which assay layer to extract feature
#'   values from (e.g., `"data"`, `"counts"`, `"scale.data"`). `NULL` uses the
#'   default layer. May be a single string (applied to all features) or a character
#'   vector of length equal to `features` (per-feature layer specification).
#' @param plot.median Logical. If `TRUE` (default), draws median line(s) in dashed
#'   black, independently of `vline`.
#' @param plot.title Optional custom title(s). `NULL` uses feature names as titles.
#'   A length-1 string applies to all features. A character vector with length equal
#'   to `length(features)` applies one title per feature.
#' @param nmad Numeric. Number of MADs used when `vline` is `"upper"`, `"lower"`,
#'   or `"both"`. Default is `2`.
#' @param alpha Numeric in `[0, 1]`. Fill alpha for density geometries.
#' @param pt.size Numeric. If `0` (default), no rug is drawn. If `> 0`, adds a rug
#'   (`geom_rug()`) with this line width.
#' @param common.scales Logical. If `TRUE`, all sub-plots share the same x and y axis limits. If `FALSE`,
#'   each sub-plot has its own axis limits. Default is `TRUE`.
#' @param collect.axes Logical. If `TRUE`, collects the axes across all sub-plots when
#'   `common.scales = TRUE`. Default is `FALSE`.
#'
#' @return If `length(features) == 1`, returns a `ggplot2`/`patchwork` plot object.
#'   If `length(features) > 1`, returns a named list of plot objects (one per feature;
#'   list names equal `features`).
#'
#' @details
#' **Data access**: Features are resolved via `Seurat::FetchData()`, which tries
#' metadata columns first, then keyed variables, then assay features, then the
#' special `ident` keyword.
#'
#' **Layer handling**: When `layer` is a vector of different values, each feature
#' is fetched with its assigned layer. Metadata columns ignore the layer parameter;
#' assay features use it to select the desired matrix.
#'
#' **Efficiency**: A top-level `lapply()` iterates over features. A second nested
#' `lapply()` is used only when `split.plot = TRUE` to iterate over group levels.
#' This structure keeps per-feature data local and reduces indexing overhead.
#'
#' **Vline semantics**: Keyword inputs ("mean", "median", "upper", "lower", "both")
#' compute reference positions from the empirical distribution per feature or per
#' group. Numeric inputs are used as explicit x-intercepts. When a per-group vector
#' is supplied for a single feature, each group panel uses its own element.
#'
#' @examples
#' \dontrun{
#' # Metadata-only (backward compatible)
#' plot.list <- FeatureDensityPlot(
#'   SeuratObject,
#'   features = c("percent.mt", "nCount_RNA"),
#'   group.by = "batch",
#'   split.plot = TRUE,
#'   vline = "upper",
#'   plot.median = TRUE
#' )
#'
#' # Assay features from log-normalized data
#' FeatureDensityPlot(
#'   SeuratObject,
#'   features = "CD3D",
#'   layer = "data"
#' )
#'
#' # Mixed per-feature layers
#' FeatureDensityPlot(
#'   SeuratObject,
#'   features = c("percent.mt", "CD3D", "nCount_RNA"),
#'   layer = c(NA, "counts", NA)  # CD3D from raw counts, others from metadata
#' )
#'
#' # Single feature with per-group vlines (one entry per batch level)
#' FeatureDensityPlot(
#'   SeuratObject,
#'   features = "percent.mt",
#'   group.by = "batch",
#'   split.plot = TRUE,
#'   vline = c("upper", "lower")
#' )
#'
#' # Single feature with numeric vlines per group
#' FeatureDensityPlot(
#'   SeuratObject,
#'   features = "CD3D",
#'   group.by = "ident",
#'   split.plot = TRUE,
#'   vline = c(0.5, 1.0, 1.5),
#'   collect.axes = TRUE
#' )
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @import patchwork
#' @import Seurat
#' @export
FeatureDensityPlot <- function(
  SeuratObject,
  features,
  group.by = "ident",
  split.plot = TRUE,
  scale.colors = "viridis",
  ncol = NULL,
  vline = NULL,
  layer = NULL,
  plot.median = TRUE,
  plot.title = NULL,
  nmad = 2,
  alpha = 0.3,
  pt.size = 0,
  common.scales = TRUE,
  collect.axes = FALSE
) {
  # ─────────────────────────────────────────────────────────────────────────────
  # 1) Input validation
  #
  # Validate all public inputs up front so that the plotting code below can
  # assume consistent, already-normalized objects and stay easier to read.
  # ─────────────────────────────────────────────────────────────────────────────
  if (!inherits(SeuratObject, "Seurat")) {
    stop("'SeuratObject' must be a Seurat object.")
  }

  if (!is.character(features) || length(features) == 0L) {
    stop("'features' must be a non-empty character vector.")
  }

  if (
    !is.logical(split.plot) || length(split.plot) != 1L || is.na(split.plot)
  ) {
    stop("'split.plot' must be a single logical value.")
  }

  if (
    !is.logical(plot.median) || length(plot.median) != 1L || is.na(plot.median)
  ) {
    stop("'plot.median' must be a single logical value.")
  }

  if (
    !is.numeric(alpha) ||
      length(alpha) != 1L ||
      is.na(alpha) ||
      alpha < 0 ||
      alpha > 1
  ) {
    stop("'alpha' must be a single numeric value between 0 and 1.")
  }

  if (
    !is.numeric(pt.size) ||
      length(pt.size) != 1L ||
      is.na(pt.size) ||
      pt.size < 0
  ) {
    stop("'pt.size' must be a single non-negative numeric value.")
  }

  if (!is.numeric(nmad) || length(nmad) != 1L || is.na(nmad) || nmad < 0) {
    stop("'nmad' must be a single non-negative numeric value.")
  }

  if (!is.null(plot.title)) {
    if (!is.character(plot.title)) {
      stop("'plot.title' must be NULL or a character vector.")
    }
    if (!(length(plot.title) == 1L || length(plot.title) == length(features))) {
      stop("'plot.title' must have length 1 or length(features).")
    }
  }

  if (!is.character(scale.colors) || length(scale.colors) != 1L) {
    stop("'scale.colors' must be a single character value.")
  }

  if (
    !is.logical(common.scales) ||
      length(common.scales) != 1L ||
      is.na(common.scales)
  ) {
    stop("'common.scales' must be a single logical value.")
  }

  if (
    !is.logical(collect.axes) ||
      length(collect.axes) != 1L ||
      is.na(collect.axes)
  ) {
    stop("'collect.axes' must be a single logical value.")
  }

  # Normalize a single vline entry to one of the internal representations used
  # downstream: NULL, a keyword string, or a numeric scalar.
  normalizeVlineEntry <- function(
    x,
    valid.keywords = c("mean", "median", "upper", "lower", "both")
  ) {
    # Treat explicit NULL-like character values as no reference line.
    if (
      is.null(x) || (is.character(x) && tolower(x) %in% c("null", "na", ""))
    ) {
      return(NULL)
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
        "vline entry '%s' is not a valid keyword or numeric value.",
        x
      ))
    }

    # Numeric inputs must be scalar and non-missing.
    if (is.numeric(x)) {
      if (length(x) == 1L && !is.na(x)) {
        return(x)
      }

      stop("numeric vline must be a single non-NA value.")
    }

    # Any other type is unsupported.
    stop(sprintf(
      "vline entry has unexpected type: %s",
      paste(class(x), collapse = ", ")
    ))
  }

  # Normalize a full vline input vector against an expected target length.
  # The target length represents either the number of requested features or the
  # number of group levels, depending on the plotting mode.
  normalizeVlineVector <- function(
    vline,
    target.length,
    valid.keywords = c("mean", "median", "upper", "lower", "both"),
    context = "features"
  ) {
    # A NULL vline means that no reference lines should be drawn for any target.
    if (is.null(vline)) {
      return(rep(list(NULL), target.length))
    }

    # Numeric inputs are accepted as either one value recycled to all targets
    # or one value per target.
    if (is.numeric(vline)) {
      if (length(vline) == 1L) {
        entry <- normalizeVlineEntry(vline, valid.keywords)
        return(rep(list(entry), target.length))
      }

      if (length(vline) == target.length) {
        return(lapply(
          as.list(vline),
          normalizeVlineEntry,
          valid.keywords = valid.keywords
        ))
      }

      stop(sprintf(
        "'vline' numeric input has length %d but %s has length %d; it must be length 1 or %d.",
        length(vline),
        context,
        target.length,
        target.length
      ))
    }

    # Character inputs follow the same contract as numeric inputs, but each
    # element may be either a keyword or a numeric value encoded as a string.
    if (is.character(vline)) {
      if (length(vline) == 1L) {
        entry <- normalizeVlineEntry(vline, valid.keywords)
        return(rep(list(entry), target.length))
      }

      if (length(vline) == target.length) {
        return(lapply(
          vline,
          normalizeVlineEntry,
          valid.keywords = valid.keywords
        ))
      }

      stop(sprintf(
        "'vline' character input has length %d but %s has length %d; it must be length 1 or %d.",
        length(vline),
        context,
        target.length,
        target.length
      ))
    }

    stop("'vline' must be NULL, numeric, or character.")
  }

  # Store allowed keyword vline specifications once so they can be reused both
  # for feature-level and group-level normalization.
  valid.vline.values <- c("mean", "median", "upper", "lower", "both")

  # Defer vline normalization until grouping has been resolved. This is required
  # because a single-feature split plot may accept one vline entry per group
  # level rather than one entry per feature.
  vline.per.feature <- NULL

  # Normalize layer input to one layer assignment per requested feature.
  # Metadata variables will ignore the layer later, but keeping the vector
  # aligned to features simplifies the rest of the code.
  if (is.null(layer)) {
    layer.per.feature <- rep(list(NULL), length(features))
  } else if (is.character(layer)) {
    if (length(layer) == 1L) {
      layer.per.feature <- rep(list(layer), length(features))
    } else if (length(layer) == length(features)) {
      layer.per.feature <- lapply(layer, function(x) {
        if (is.na(x) || tolower(x) == "null" || tolower(x) == "") {
          return(NULL)
        }
        x
      })
    } else {
      stop(sprintf(
        "'layer' has length %d but 'features' has length %d; it must be length 1 or length(features).",
        length(layer),
        length(features)
      ))
    }
  } else {
    stop("'layer' must be NULL or a character vector.")
  }

  # ─────────────────────────────────────────────────────────────────────────────
  # 2) Data extraction via Seurat::FetchData
  #
  # We construct per-feature internal keys so duplicated feature names remain
  # distinct internally. This is essential for cases like features = c("MALAT1",
  # "MALAT1") with different layer entries; name-based indexing would collapse
  # those entries and silently reuse one layer.
  # ─────────────────────────────────────────────────────────────────────────────
  cells.use <- colnames(SeuratObject)
  has.grouping <- !is.null(group.by)

  # Build one unique internal key per requested feature occurrence.
  feature.keys <- paste0(".feature_", seq_along(features))

  # Resolve grouping values once so every feature is fetched against the same
  # cell order and all downstream joins remain positionally aligned.
  if (has.grouping) {
    if (!is.character(group.by) || length(group.by) != 1L) {
      stop("'group.by' must be NULL or a single character value.")
    }

    group.fetch.var <- group.by
    group.data <- Seurat::FetchData(
      object = SeuratObject,
      vars = group.fetch.var,
      cells = cells.use,
      clean = FALSE
    )

    group.values <- as.character(group.data[[group.fetch.var]])
    group.label <- group.by
  } else {
    group.fetch.var <- NULL
    group.values <- rep("all", length(cells.use))
    group.label <- "all"
  }

  # Fetch each feature independently, keeping feature position as the primary
  # contract that links features, layer assignments, titles, and vline specs.
  feature.data <- lapply(seq_along(features), function(feature.id) {
    feature.name <- features[[feature.id]]
    feature.layer <- layer.per.feature[[feature.id]]

    fetched <- Seurat::FetchData(
      object = SeuratObject,
      vars = feature.name,
      cells = cells.use,
      layer = feature.layer,
      clean = FALSE
    )

    if (!feature.name %in% colnames(fetched)) {
      stop(sprintf(
        "The feature '%s' could not be fetched from the object. Check the feature name and layer.",
        feature.name
      ))
    }

    # Reindex explicitly to the input cell order so downstream columns remain
    # aligned even if FetchData changes row ordering internally in the future.
    fetched[cells.use, feature.name, drop = TRUE]
  })
  names(feature.data) <- feature.keys

  # Assemble the plotting data frame using internal keys so duplicate feature
  # names can coexist without overwriting one another.
  plot.data <- data.frame(row.names = cells.use)
  for (feature.id in seq_along(features)) {
    feature.key <- feature.keys[[feature.id]]
    plot.data[[feature.key]] <- feature.data[[feature.key]]
  }

  # Attach grouping information and drop cells without a valid grouping value.
  plot.data$.group <- group.values
  plot.data <- plot.data[!is.na(plot.data$.group), , drop = FALSE]

  # Density plots require numeric x values. Coercion is explicit, and any
  # non-convertible values become NA and are removed later on a per-feature basis.
  for (feature.id in seq_along(features)) {
    feature.key <- feature.keys[[feature.id]]
    if (!is.numeric(plot.data[[feature.key]])) {
      suppressWarnings(
        plot.data[[feature.key]] <- as.numeric(plot.data[[feature.key]])
      )
    }
  }

  # Build user-facing plot names. This preserves layer-specific distinctions for
  # repeated assay features while keeping metadata-derived names clean.
  metadata.colnames <- colnames(SeuratObject@meta.data)

  layer.labels <- vapply(
    layer.per.feature,
    function(x) {
      if (is.null(x)) {
        return("NULL")
      }
      as.character(x)
    },
    FUN.VALUE = character(1)
  )

  plot.names <- paste0(features, "_", layer.labels)

  # Metadata variables should not carry an assay-layer suffix in the returned
  # plot names because that suffix is not meaningful to the user in that case.
  is.metadata.feature <- features %in% metadata.colnames
  plot.names[is.metadata.feature] <- sub(
    "_[^_]*$",
    "",
    plot.names[is.metadata.feature]
  )

  # Retain legacy cleanup for edge cases involving empty, NA, or NULL-like
  # layer labels.
  plot.names <- gsub("_+$", "", plot.names)
  plot.names <- gsub("_NULL$", "", plot.names)
  plot.names <- gsub("_NA$", "", plot.names)

  # Resolve final plot titles in the same order as the requested features.
  feature.titles <- if (is.null(plot.title)) {
    plot.names
  } else if (length(plot.title) == 1L) {
    rep(plot.title, length(features))
  } else {
    plot.title
  }

  # Cache the encountered group levels once for split plotting.
  group.levels <- unique(plot.data$.group)
  if (has.grouping && length(group.levels) == 0L) {
    stop("No groups available to plot after filtering missing grouping values.")
  }

  # Normalize vline only after group levels are known. This allows a
  # single-feature split plot to accept one vline entry per group, while all
  # other modes keep the standard one-entry-per-feature contract.
  vline.per.group <- NULL
  if (has.grouping && split.plot && length(features) == 1L) {
    vline.per.group <- normalizeVlineVector(
      vline = vline,
      target.length = length(group.levels),
      valid.keywords = valid.vline.values,
      context = "groups"
    )

    # Keep a feature-level placeholder so indexing by feature position still
    # works in the main plotting loop.
    vline.per.feature <- list(NULL)
  } else {
    vline.per.feature <- normalizeVlineVector(
      vline = vline,
      target.length = length(features),
      valid.keywords = valid.vline.values,
      context = "features"
    )
  }

  # ─────────────────────────────────────────────────────────────────────────────
  # 3) Helper functions
  # ─────────────────────────────────────────────────────────────────────────────

  # Convert a normalized vline specification into one or more explicit
  # x-intercepts for the current numeric vector.
  computeVlinePositions <- function(values, vline.spec, nmad.value) {
    values <- values[!is.na(values)]

    if (length(values) == 0L || is.null(vline.spec)) {
      return(data.frame(xintercept = numeric(0), stringsAsFactors = FALSE))
    }

    # Numeric specifications are already concrete intercepts.
    if (is.numeric(vline.spec)) {
      return(data.frame(xintercept = vline.spec, stringsAsFactors = FALSE))
    }

    median.value <- stats::median(values)

    if (vline.spec == "mean") {
      return(data.frame(xintercept = mean(values), stringsAsFactors = FALSE))
    }

    if (vline.spec == "median") {
      return(data.frame(xintercept = median.value, stringsAsFactors = FALSE))
    }

    # MAD-based thresholds require a valid MAD estimate. If MAD is unavailable,
    # no keyword-based threshold line is returned.
    mad.value <- stats::mad(values, na.rm = TRUE)
    if (is.na(mad.value)) {
      return(data.frame(xintercept = numeric(0), stringsAsFactors = FALSE))
    }

    if (vline.spec == "upper") {
      return(data.frame(
        xintercept = median.value + nmad.value * mad.value,
        stringsAsFactors = FALSE
      ))
    }

    if (vline.spec == "lower") {
      return(data.frame(
        xintercept = median.value - nmad.value * mad.value,
        stringsAsFactors = FALSE
      ))
    }

    if (vline.spec == "both") {
      return(data.frame(
        xintercept = c(
          median.value - nmad.value * mad.value,
          median.value + nmad.value * mad.value
        ),
        stringsAsFactors = FALSE
      ))
    }

    return(data.frame(xintercept = numeric(0), stringsAsFactors = FALSE))
  }

  # Compute the single median x-intercept used for the optional black reference
  # line that is independent from the red vline specification.
  computeMedianPositions <- function(values) {
    values <- values[!is.na(values)]

    if (length(values) == 0L) {
      return(data.frame(xintercept = numeric(0), stringsAsFactors = FALSE))
    }

    data.frame(
      xintercept = stats::median(values),
      stringsAsFactors = FALSE
    )
  }

  # Build one density panel. This helper is reused for ungrouped plots and for
  # each per-group panel in split mode so the drawing logic lives in one place.
  buildSingleDensityPanel <- function(
    feature.df,
    x.label,
    panel.title,
    vline.spec
  ) {
    current.plot <- ggplot2::ggplot(feature.df, ggplot2::aes(x = value)) +
      ggplot2::geom_density(
        fill = "lightblue",
        alpha = alpha,
        color = "black"
      )

    # Add an optional rug to show the observed point distribution on the x-axis.
    if (pt.size > 0) {
      current.plot <- current.plot +
        ggplot2::geom_rug(
          sides = "b",
          linewidth = pt.size,
          alpha = 0.35
        )
    }

    # Add red reference line(s) derived from the normalized vline specification.
    vline.df <- computeVlinePositions(feature.df$value, vline.spec, nmad)
    if (nrow(vline.df) > 0L) {
      current.plot <- current.plot +
        ggplot2::geom_vline(
          data = vline.df,
          ggplot2::aes(xintercept = xintercept),
          color = "red",
          linetype = "dashed",
          linewidth = 0.6,
          show.legend = FALSE
        )
    }

    # Add black median line(s) independently when requested.
    if (plot.median) {
      median.df <- computeMedianPositions(feature.df$value)
      if (nrow(median.df) > 0L) {
        current.plot <- current.plot +
          ggplot2::geom_vline(
            data = median.df,
            ggplot2::aes(xintercept = xintercept),
            color = "black",
            linetype = "dashed",
            linewidth = 0.7,
            show.legend = FALSE
          )
      }
    }

    current.plot +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = panel.title,
        x = x.label,
        y = "Density"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5)
      )
  }

  # ─────────────────────────────────────────────────────────────────────────────
  # 4) Feature-level plotting
  # ─────────────────────────────────────────────────────────────────────────────
  feature.plots <- lapply(seq_along(features), function(feature.id) {
    # Resolve all feature-specific inputs by position so duplicate feature names
    # remain distinguishable internally.
    feature.name <- features[[feature.id]]
    feature.key <- feature.keys[[feature.id]]
    feature.title <- feature.titles[[feature.id]]
    feature.vline <- vline.per.feature[[feature.id]]

    # Extract the current feature column plus grouping information.
    feature.df <- plot.data[, c(feature.key, ".group"), drop = FALSE]
    names(feature.df)[1] <- "value"

    # Remove non-numeric or missing values before plotting.
    feature.df <- feature.df[!is.na(feature.df$value), , drop = FALSE]

    if (nrow(feature.df) == 0L) {
      stop(sprintf(
        "Feature '%s' has no non-missing numeric values to plot.",
        feature.name
      ))
    }

    # Ungrouped mode always returns a single panel.
    if (!has.grouping) {
      return(buildSingleDensityPanel(
        feature.df = feature.df,
        x.label = feature.name,
        panel.title = feature.title,
        vline.spec = feature.vline
      ))
    }

    # Overlay mode keeps all groups in one panel and computes red/black
    # reference lines independently within each group.
    if (!split.plot) {
      current.plot <- ggplot2::ggplot(
        feature.df,
        ggplot2::aes(x = value, color = .group, fill = .group)
      ) +
        ggplot2::geom_density(alpha = alpha)

      if (pt.size > 0) {
        current.plot <- current.plot +
          ggplot2::geom_rug(
            ggplot2::aes(color = .group),
            sides = "b",
            linewidth = pt.size,
            alpha = 0.35,
            show.legend = FALSE
          )
      }

      # In overlay mode, vline is interpreted at feature level and then evaluated
      # separately within each group.
      if (!is.null(feature.vline)) {
        vline.df <- feature.df |>
          dplyr::group_by(.group) |>
          dplyr::group_modify(
            ~ computeVlinePositions(.x$value, feature.vline, nmad)
          ) |>
          dplyr::ungroup()

        if (nrow(vline.df) > 0L) {
          current.plot <- current.plot +
            ggplot2::geom_vline(
              data = vline.df,
              ggplot2::aes(xintercept = xintercept, group = .group),
              color = "red",
              linetype = "dashed",
              linewidth = 0.55,
              show.legend = FALSE
            )
        }
      }

      if (plot.median) {
        median.df <- feature.df |>
          dplyr::group_by(.group) |>
          dplyr::summarise(
            xintercept = stats::median(value, na.rm = TRUE),
            .groups = "drop"
          )

        if (nrow(median.df) > 0L) {
          current.plot <- current.plot +
            ggplot2::geom_vline(
              data = median.df,
              ggplot2::aes(xintercept = xintercept, group = .group),
              color = "black",
              linetype = "dashed",
              linewidth = 0.65,
              show.legend = FALSE
            )
        }
      }

      return(
        current.plot +
          ggplot2::scale_color_viridis_d(
            option = scale.colors,
            name = group.label
          ) +
          ggplot2::scale_fill_viridis_d(
            option = scale.colors,
            name = group.label
          ) +
          ggplot2::theme_bw() +
          ggplot2::labs(
            title = feature.title,
            x = feature.name,
            y = "Density"
          ) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5)
          )
      )
    }

    # Split mode builds one panel per group. When exactly one feature is plotted,
    # each panel may receive its own vline specification via vline.per.group.
    group.plots <- lapply(seq_along(group.levels), function(group.index) {
      group.name <- group.levels[[group.index]]
      group.df <- feature.df[feature.df$.group == group.name, , drop = FALSE]

      panel.vline <- if (!is.null(vline.per.group) && length(features) == 1L) {
        vline.per.group[[group.index]]
      } else {
        feature.vline
      }

      buildSingleDensityPanel(
        feature.df = group.df,
        x.label = feature.name,
        panel.title = group.name,
        vline.spec = panel.vline
      )
    })

    # Use the user-provided ncol when available; otherwise infer a roughly square
    # layout for the split panels.
    ncol.groups <- if (!is.null(ncol)) {
      ncol
    } else {
      ceiling(sqrt(length(group.plots)))
    }

    # Combine the split panels into one feature-level plot with a shared title.
    # The title is centered and bolded for better visibility.
    combined <- patchwork::wrap_plots(group.plots, ncol = ncol.groups) +
      patchwork::plot_annotation(
        title = feature.title,
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
        )
      )

    # Apply a shared x and y axis scale if requested in common.scales option.
    if (common.scales) {
      combined <- .SetCommonScales(combined)
    }

    # Apply axes collection if requested in collect.axes option.
    if (collect.axes) {
      combined <- combined + patchwork::plot_layout(axes = "collect")
    }

    return(combined)
  })

  # ─────────────────────────────────────────────────────────────────────────────
  # 5) Return-shape contract
  # ─────────────────────────────────────────────────────────────────────────────

  # Preserve user-facing names in the returned list.
  names(feature.plots) <- plot.names

  # Return a single plot object for one feature to keep the API ergonomic and
  # backward compatible.
  if (length(feature.plots) == 1L) {
    return(feature.plots[[1L]])
  }

  # Return a named list for multi-feature calls.
  return(feature.plots)
}
