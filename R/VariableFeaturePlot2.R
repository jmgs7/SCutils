#' @title VariableFeaturePlot2
#'
#' @description
#' This function is a modified version of Seurat's VariableFeaturePlot.
#' It allows users to specify a custom list of highly variable features (hvf) for plotting,
#' rather than relying solely on the variable features stored in the Seurat object. It also forces
#' the use of the output of SeuratObject::VariableFeatures() to highligh the high variable features
#' when the hvf parameter is not provided. This can be useful for visualizing specific sets of
#' features, for example, those which have been filtered based on certain criteria.
#'
#' @inheritParams Seurat::FeatureScatter
#' @inheritParams SeuratObject::HVFInfo
#' @param cols Colors to specify non-variable/variable status
#' @param assay Assay to pull variable features from
#' @param log Plot the x-axis in log scale
#' @param raster Convert points to raster format, default is \code{NULL}
#' which will automatically use raster if the number of points plotted is greater than
#' 100,000
#' @param hvf The list of highly variable features to use for plotting. If \code{NULL}, the function
#'   will use the variable features stored in the Seurat object.
#' @param label Logical, whether to label the top variable features on the plot. Default is \code
#'   {FALSE}.
#' @param n.top.hvf Number of top variable features to label if \code{label} is \code{TRUE}.
#'   Default is 10.
#' @param custom.label A character vector of specific features to label on the plot. If provided,
#'   these features will be labeled instead of the top variable features. Default is \code{NULL}.
#' @param repel Logical, whether to use ggrepel for labeling points. Default is \code{TRUE}.
#'
#' @details
#' This function is designed to provide flexibility in visualizing variable features in a Seurat
#' object. By allowing users to specify a custom list of highly variable features, it enables the
#' visualization of specific features of interest, rather than being limited to the variable
#' features identified by Seurat's default methods.
#'
#' @return A ggplot object
#'
#' @import Seurat
#' @import SeuratObject
#' @importFrom ggplot2 labs scale_color_manual scale_x_log10
#' @importFrom dplyr filter arrange slice_head
#' @export
#'
#' @examples
#' \dontrun{
#'   data("pbmc_small")
#'   pbmc_small <- data("pbmc_small") |>
#'     Seurat::NormalizeData() |>
#'     Seurat::FindVariableFeatures() |>
#'     SCutils::FilterVariableFeatures() |>
#'     Seurat::ScaleData()
#'   hvf <- SeuratObject::VariableFeatures(pbmc_small)
#'   # Plot a unique variable feature plot for all layers but highligh
#'   # the consensus variable features across all layers instead of the r
#'   # recalculated variable features for each layer.
#'   pbmc_small |>
#'     JoinLayers() |>
#'     Seurat::FindVariableFeatures() |>
#'     SCutils::VariableFeaturePlot2(hvf = hvf)
#' }
#'
VariableFeaturePlot2 <- function(
  object,
  cols = c('black', 'red'),
  pt.size = 1,
  log = NULL,
  selection.method = NULL,
  assay = NULL,
  raster = NULL,
  raster.dpi = c(512, 512),
  hvf = NULL,
  label = FALSE,
  n.top.hvf = 10,
  custom.label = NULL,
  repel = TRUE
) {
  if (length(x = cols) != 2) {
    stop("'cols' must be of length 2")
  }
  hvf.df <- SeuratObject::HVFInfo(
    object = object,
    assay = assay,
    method = selection.method,
    status = TRUE
  )

  status.col <- colnames(hvf.df)[grepl("variable", colnames(hvf.df))][[1]]

  if (is.null(x = hvf)) {
    hvf <- SeuratObject::VariableFeatures(
      object = object,
      assay = assay,
      method = selection.method
    )
    # Ensure that the filtered variables have their status set to 'no'.
    hvf.df[!rownames(x = hvf.df) %in% hvf, status.col] <- FALSE
  } else {
    if (!all(hvf %in% rownames(x = hvf.df))) {
      stop("Some features in 'hvf' are not present in the hvf info dataframe.")
    }
    # Reset the status of all features to FALSE.
    hvf.df[rownames(x = hvf.df), status.col] <- FALSE
    # Set the status of the specified hvf features to TRUE.
    hvf.df[rownames(x = hvf.df) %in% hvf, status.col] <- TRUE
  }

  var.status <- c('no', 'yes')[unlist(hvf.df[[status.col]]) + 1]

  if (colnames(x = hvf.df)[3] == 'dispersion.scaled') {
    hvf.info <- hvf.df[, c(1, 2)]
  } else if (colnames(x = hvf.df)[3] == 'variance.expected') {
    hvf.info <- hvf.df[, c(1, 4)]
  } else {
    hvf.info <- hvf.df[, c(1, 3)]
  }
  axis.labels <- switch(
    EXPR = colnames(x = hvf.info)[2],
    'variance.standardized' = c('Average Expression', 'Standardized Variance'),
    'dispersion' = c('Average Expression', 'Dispersion'),
    'residual_variance' = c('Geometric Mean of Expression', 'Residual Variance')
  )
  log <- log %||%
    (any(
      c('variance.standardized', 'residual_variance') %in%
        colnames(x = hvf.info)
    ))

  plot <- Seurat::SingleCorPlot(
    data = hvf.info,
    col.by = var.status,
    pt.size = pt.size,
    raster = raster,
    raster.dpi = raster.dpi
  )
  if (length(x = unique(x = var.status)) == 1) {
    switch(
      EXPR = var.status[1],
      'yes' = {
        cols <- cols[2]
        labels.legend <- 'Variable'
      },
      'no' = {
        cols <- cols[1]
        labels.legend <- 'Non-variable'
      }
    )
  } else {
    labels.legend <- c('Non-variable', 'Variable')
  }
  plot <- plot +
    labs(title = NULL, x = axis.labels[1], y = axis.labels[2]) +
    scale_color_manual(
      labels = paste(labels.legend, 'count:', table(var.status)),
      values = cols
    )
  if (log) {
    plot <- plot + scale_x_log10()
  }

  if (label) {
    if (is.null(x = custom.label)) {
      top.hvf <- hvf.df |>
        dplyr::filter(.data[[status.col]] == TRUE) |>
        dplyr::arrange(.data[["rank"]]) |>
        dplyr::slice_head(n = n.top.hvf) |>
        rownames()
      plot <- LabelPoints(plot = plot, points = top.hvf, repel = repel)
    } else {
      # check all custom features are present
      present.features <- custom.label[custom.label %in% rownames(x = hvf.df)]
      if (length(x = present.features) != length(x = custom.label)) {
        warning(
          "Some features in 'custom.label' were not found in the hvf info dataframe and will not be labeled."
        )
      }
      plot <- LabelPoints(
        plot = plot,
        points = present.features,
        repel = repel
      )
    }
  }

  return(plot)
}
