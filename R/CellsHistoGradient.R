#' Create a Histogram with Gradient Colors for Cell Counts.
#' 
#' @name CellsHistoGradient
#' @description Generates a bar plot showing the number of cells per group in a Seurat object, with bars colored by a gradient based on cell count.
#'
#' @param SeuratObject A Seurat object containing single-cell data
#' @param group.by Character string specifying the metadata column to group cells by.
#'   If NULL (default), uses the active identity of the Seurat object.
#' @param scale.colors Color pallette to use in the gradient. Corresponds to the color palletes codes included in scale_gradient_viridis: "magma" (or "A"), "inferno" (or "B"), "plasma" (or "C"), "viridis" (or "D"), "cividis" (or "E"), "rocket" (or "F"), "mako" (or "G"), turbo" (or "H")
#' @param breaks Function or vector specifying the breaks for the y-axis.
#'   Default is scales::extended_breaks().
#'
#' @return A ggplot object displaying a bar plot of cell counts per group
#'
#' @examples
#' \dontrun{
#' # Use active identity
#' CellsHistoGradient(seurat_obj)
#'
#' # Group by specific metadata column
#' CellsHistoGradient(seurat_obj, group.by = "cell_type")
#'
#' # Customize color scale
#' CellsHistoGradient(seurat_obj, group.by = "cluster", scale.colors = "plasma")
#'
#' # Customize breaks
#' CellsHistoGradient(seurat_obj, group.by = "annotation", breaks = seq(0, 1000, by = 100))
#' }
#' 
#' @import dplyr
#' @import ggplot2
#' @import Seurat
#' 
#' @export

CellsHistoGradient <- function(SeuratObject, group.by = NULL, scale.colors = "viridis", breaks = scales::extended_breaks()) {

  if (is.null(group.by)) {

    group.by <- "SeuratObject@active.ident"
    cells_per_feature <- as.data.frame(SeuratObject@active.ident) %>%
      group_by(.data[[group.by]]) %>%
      summarise(Cells = n()) %>%
      arrange(desc(Cells))

  } else {

    cells_per_feature <- SeuratObject@meta.data %>%
      group_by(.data[[group.by]]) %>%
      summarise(Cells = n()) %>%
      arrange(desc(Cells))
  }

  p <- ggplot(
    cells_per_feature,
    aes(
      x = factor(.data[[group.by]], levels = .data[[group.by]]),
      y = Cells,
      fill = Cells
    )
  ) +
    xlab(group.by) +
    scale_y_continuous(breaks = breaks) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_viridis_c(name = group.by, option = scale.colors) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    )
    
  return(p)
}
