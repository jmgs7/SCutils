#' @title CalculateQC
#' @description
#' Estimate common single-cell QC metrics and append them to a Seurat object's metadata.
#'
#' This function calculates percent-based metrics using Seurat::PercentageFeatureSet for
#' mitochondrial, ribosomal, hemoglobin, immunoglobulin, platelet-associated genes and
#' several individual marker genes, and also computes log10-transformed feature/count
#' values and a complexity ratio.
#'
#' If the Seurat object contains log-normalized data layers, the function will also calculate
#' the total number of log-normalized counts and features per cell, as well as cell cycle scoring
#' using updated S and G2/M phase gene sets.
#'
#' #' @details The following metadata columns are added to SeuObj:
#' \itemize{
#'   \item percent.mt: percentage of counts matching '^MT-'
#'   \item percent.ribo: percentage of counts matching '^RP[SL]'
#'   \item percent.hb: percentage of counts matching '^HB[^(P)]'
#'   \item percent.ig: percentage of counts matching '^IG'
#'   \item percent.plat: percentage of counts matching 'PECAM1|PF4'
#'   \item percent.MALAT1, percent.S100A9, percent.S100A8, percent.FCGR3B: individual gene percentages
#'   \item log10_nFeature_RNA, log10_nCount_RNA: base-10 logarithms of nFeature_RNA and nCount_RNA
#'   \item complexity: ratio log10_nFeature_RNA / log10_nCount_RNA
#' }
#'
#' @param SeuObj A Seurat object. The function reads counts from the active assay and adds metadata columns.
#'   Patterns are interpreted as regular expressions by Seurat::PercentageFeatureSet.
#' @param data.layers Default: `NULL`. If normalization has been applied to the SeuObj, you can provide
#'   the layers where the data is stored. If NULL, the function will attempt to find all log-normalized data layers
#'   in the Seurat object. If not found, the function will skip log-normalized QC calculations.
#' @param perform.cell.cycle.scoring Default: `TRUE`. If `TRUE`, the function will perform cell cycle scoring
#'   using the updated S and G2/M phase gene sets. If `FALSE`, cell cycle scoring will be skipped.
#' @return The input Seurat object with the new metadata columns added.
#' @examples
#' # CalculateQC(SeuObj)
#' # CalculateQC(SeuObj, data.layers = c("data", "data_layer2"))
#' # CalculateQC(SeuObj, perform.cell.cycle.scoring = FALSE)
#' @import Seurat
#' @import SeuratObject
#' @import BPCells
#' @import data.table
#' @export

CalculateQC <- function(
  SeuObj,
  data.layers = NULL,
  perform.cell.cycle.scoring = TRUE
) {
  # Estimation of metrics
  SeuObj@meta.data$percent.mt <- Seurat::PercentageFeatureSet(
    SeuObj,
    pattern = "^MT-"
  ) # Percentage of counts corresponding to mitochondrial genes.
  SeuObj@meta.data$percent.ribo <- Seurat::PercentageFeatureSet(
    SeuObj,
    "^RP[SL]"
  ) # Percentage of counts corresponding to ribosomal genes.
  SeuObj@meta.data$percent.hb <- Seurat::PercentageFeatureSet(
    SeuObj,
    "^HB[^(P)]"
  ) # Percentage of counts corresponding to hemoglobin.
  SeuObj@meta.data$percent.ig <- Seurat::PercentageFeatureSet(
    SeuObj,
    "^IG"
  ) # Percentage of counts corresponding to immunoglobulins.
  SeuObj@meta.data$percent.plat <- Seurat::PercentageFeatureSet(
    SeuObj,
    "PECAM1|PF4"
  ) # Percentage of counts corresponding to genes associated with platelets.
  SeuObj@meta.data$percent.MALAT1 <- Seurat::PercentageFeatureSet(
    SeuObj,
    pattern = "MALAT1"
  ) # Percentage of counts corresponding to MALAT1.
  SeuObj@meta.data$percent.S100A9 <- Seurat::PercentageFeatureSet(
    SeuObj,
    pattern = "S100A9"
  ) # Percentage of counts corresponding to S100A9.
  SeuObj@meta.data$percent.S100A8 <- Seurat::PercentageFeatureSet(
    SeuObj,
    pattern = "S100A8"
  ) # Percentage of counts corresponding to S100A8.
  SeuObj@meta.data$percent.FCGR3B <- Seurat::PercentageFeatureSet(
    SeuObj,
    pattern = "FCGR3B"
  ) # Percentage of counts corresponding to FCGR3B.
  SeuObj@meta.data$log10_nFeature_RNA <- log10(
    SeuObj@meta.data$nFeature_RNA
  )
  SeuObj@meta.data$log10_nCount_RNA <- log10(
    SeuObj@meta.data$nCount_RNA
  )
  SeuObj@meta.data$complexity <- SeuObj@meta.data$log10_nFeature_RNA /
    SeuObj@meta.data$log10_nCount_RNA # Complexity, corresponding to the amount of genes that are covered by the counts of each cell.

  # This steps are only if the Seurat object has log-normalized data layers, which are not always present.
  seurat.data.layers <- SeuratObject::Layers(
    SeuObj,
    assay = "RNA",
    search = "data"
  )
  if (is.null(data.layers)) {
    # Check if the Seurat object has log-normalized data layers, and if not, skip the log-normalized QC calculations.
    if (!any(grepl("data", seurat.data.layers))) {
      message(
        "No log-normalized data layers found in the Seurat object. Skipping log-normalized QC calculations."
      )
      return(SeuObj)
    }
    # Also check if the user provided data.layers are present in the Seurat object.
  } else if (!all(data.layers %in% seurat.data.layers)) {
    stop(
      "Some specified data layers are not present in the Seurat object. Please check the provided data.layers argument."
    )
  }

  DataLayersQC <- function(
    SeuObj,
    data.layers = NULL,
    perform.cell.cycle.scoring = TRUE
  ) {
    # Helper function to calculate QC metrics for each log-normalized data layer in the Seurat object.
    # We calculate the total number of log-normalized counts per cell, and the total
    # number of features detected per cell in log counts, and cell cycle scoring.

    # 1. Fetch all log-normalized 'data' layers names.
    if (is.null(data.layers)) {
      data.layers <- SeuratObject::Layers(
        SeuObj,
        assay = "RNA",
        search = "data"
      )
    }

    # 2. Cycle per layer and calculate the total number of log-normalized counts per cell,
    # and the total number of features detected per cell in log counts, and cell cycle scoring.
    # Data will be added to the Seurat object's metadata per layer.
    log.metadata <- lapply(data.layers, function(layer) {
      # Extract the layer data to compute QC metrics.
      layer.data <- SeuratObject::LayerData(
        SeuObj,
        assay = "RNA",
        layer = layer
      )

      # Create an emtpy data.frame with the same dimensions as the number of cells in the layer data.
      # This will store the QC metrics for each cell in the current layer.
      log.metadata <- data.frame(
        cell.id = colnames(layer.data),
        nCount_logRNA = numeric(ncol(layer.data)),
        nFeature_logRNA = numeric(ncol(layer.data))
      )

      # Calculate the total number of log-normalized counts and features per cell.
      log.metadata$nCount_logRNA <- BPCells::colSums(layer.data)
      log.metadata$nFeature_logRNA <- BPCells::colSums(layer.data > 0)

      return(log.metadata)
    }) |> # Convert the list of data frames to a single data frame.
      data.table::rbindlist() |>
      as.data.frame()

    # Set the row names of the log.metadata data frame to the cell IDs for proper alignment with the Seurat object's metadata.
    row.names(log.metadata) <- log.metadata$cell.id
    # Add the calculated QC metrics to the Seurat object's metadata.
    SeuObj@meta.data$nCount_logRNA <- log.metadata$nCount_logRNA
    SeuObj@meta.data$nFeature_logRNA <- log.metadata$nFeature_logRNA

    # Add cell cycle scoring results to the Seurat object's metadata if requested.
    # Perform cell cycle scoring if requested.
    if (perform.cell.cycle.scoring) {
      SeuObj <- Seurat::CellCycleScoring(
        SeuObj,
        # Uses Seurat S and G2M genes dataset.
        s.features = cc.genes.updated.2019$s.genes,
        g2m.features = cc.genes.updated.2019$g2m.genes
      )
    }

    return(SeuObj)
  }

  # Apply the helper function to the Seurat object.
  SeuObj <- DataLayersQC(
    SeuObj,
    data.layers = data.layers,
    perform.cell.cycle.scoring = perform.cell.cycle.scoring
  )

  ##TODO: Add MALAT1 test.
  return(SeuObj)
}
