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
#' #' @details The following metadata columns are added to SeuratObject:
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
#' @param SeuratObject A Seurat object. The function reads counts from the active assay and adds metadata columns.
#'   Patterns are interpreted as regular expressions by Seurat::PercentageFeatureSet.
#' @param assay Default: "RNA". The assay in SeuratObject to use for QC calculations.
#' @param layers Default: `NULL`. If normalization has been applied to the SeuratObject, you can provide
#'   the layers where the data is stored. If NULL, the function will attempt to find all log-normalized data layers
#'   in the Seurat object. If not found, the function will skip log-normalized QC calculations.
#' @param species Default: "human". The species of the dataset. Currently supports "human" and "mouse".
#' @param perform.cell.cycle.scoring Default: `TRUE`. If `TRUE`, the function will perform cell cycle scoring
#'   using the updated S and G2/M phase gene sets. If `FALSE`, cell cycle scoring will be skipped.
#' @param perform.MALAT1.test Default: `TRUE`. If `TRUE`, the function will apply the MALAT1 thresholding
#'   function to the Seurat object. If `FALSE`, the MALAT1 thresholding will be skipped.
#' @param ... Additional arguments passed to the underlying `.CalculateFeatureThresholdSeurat()`
#'   function, such as:
#'  \itemize{
#'    \item \code{bw.bandwidth}: Bandwidth for kernel density estimation (default: 0.01).
#'    \item \code{chosen.min}: Chosen minimum which a peak should be considered the dataset
#'      peak (default: 2).
#'    \item \code{smooth.spar}: Smoothing parameter for density estimation (default: 2).
#'    \item \code{abs.min}: Absolute minimum threshold (default: 1).
#'    \item \code{rough.max}: Rough expected position of the MALAT1 expression peak (default: 6).
#'    \item \code{conservative.threshold}: Conservative threshold to apply when impossible to find local
#'      minimum (default: 2).
#'  }
#'
#' @return SeuratObject. The input Seurat object with the new metadata columns added.
#'
#' @examples
#'   \dontrun{
#'     SeuratObject <- CalculateQC(SeuratObject)
#'     SeuratObject <- CalculateQC(SeuratObject, layers = c("data", "data_layer2"))
#'     SeuratObject <- CalculateQC(SeuratObject, species = "mouse", perform.cell.cycle.scoring = FALSE)
#'     SeuratObject <- CalculateQC(SeuratObject, assay = "RNA", perform.MALAT1.test = FALSE)
#' }
#'
#' @import Seurat
#' @import SeuratObject
#' @import BPCells
#' @importFrom data.table rbindlist
#'
#' @export
CalculateQC <- function(
  SeuratObject,
  assay = "RNA",
  layers = NULL,
  species = "human",
  perform.cell.cycle.scoring = TRUE,
  perform.MALAT1.test = TRUE,
  ...
) {
  # Input validation
  species <- tolower(species)
  if (!species %in% c("human", "mouse")) {
    stop(
      "Invalid species provided. Please use 'human' or 'mouse'."
    )
  }

  if (!is.character(assay) || length(assay) != 1) {
    stop(
      "Invalid assay provided. Please provide a single assay name as a character string."
    )
  }

  if (!is.null(layers) && (!is.character(layers) || length(layers) < 1)) {
    stop(
      "Invalid layers provided. Please provide a character vector of layer names or NULL."
    )
  }

  if (
    !is.logical(perform.cell.cycle.scoring) ||
      length(perform.cell.cycle.scoring) != 1
  ) {
    stop(
      "Invalid perform.cell.cycle.scoring provided. Please provide a single logical value (TRUE or FALSE)."
    )
  }

  if (!is.logical(perform.MALAT1.test) || length(perform.MALAT1.test) != 1) {
    stop(
      "Invalid perform.MALAT1.test provided. Please provide a single logical value (TRUE or FALSE)."
    )
  }

  # Generate pattern dataframe for select the appropriate gene patterns for each species.
  patterns <- data.frame(
    species = c("human", "mouse"),
    mt.pattern = c("^MT-", "^mt-"),
    ribo.pattern = c("^RP[SL]", "^Rp[sl]"),
    hb.pattern = c("^HB[^(P)]", "^Hb[^(p)]"),
    ig.pattern = c("^IG", "^Ig"),
    plat.pattern = c("PECAM1|PF4", "Pecam1|Pf4"),
    S100A9.pattern = c("S100A9", "S100a9"),
    S100A8.pattern = c("S100A8", "S100a8"),
    FCGR3B.pattern = c("FCGR3B", "Fcgr3b"),
    MALAT1.pattern = c("MALAT1", "Malat1")
  )
  rownames(patterns) <- patterns$species

  # Estimation of metrics
  SeuratObject@meta.data$percent.mt <- Seurat::PercentageFeatureSet(
    SeuratObject,
    pattern = patterns[species, "mt.pattern"]
  ) # Percentage of counts corresponding to mitochondrial genes.
  SeuratObject@meta.data$percent.ribo <- Seurat::PercentageFeatureSet(
    SeuratObject,
    pattern = patterns[species, "ribo.pattern"]
  ) # Percentage of counts corresponding to ribosomal genes.
  SeuratObject@meta.data$percent.hb <- Seurat::PercentageFeatureSet(
    SeuratObject,
    pattern = patterns[species, "hb.pattern"]
  ) # Percentage of counts corresponding to hemoglobin.
  SeuratObject@meta.data$percent.ig <- Seurat::PercentageFeatureSet(
    SeuratObject,
    pattern = patterns[species, "ig.pattern"]
  ) # Percentage of counts corresponding to immunoglobulins.
  SeuratObject@meta.data$percent.plat <- Seurat::PercentageFeatureSet(
    SeuratObject,
    pattern = patterns[species, "plat.pattern"]
  ) # Percentage of counts corresponding to genes associated with platelets.
  SeuratObject@meta.data$percent.MALAT1 <- Seurat::PercentageFeatureSet(
    SeuratObject,
    pattern = patterns[species, "MALAT1.pattern"]
  ) # Percentage of counts corresponding to MALAT1.
  SeuratObject@meta.data$percent.S100A9 <- Seurat::PercentageFeatureSet(
    SeuratObject,
    pattern = patterns[species, "S100A9.pattern"]
  ) # Percentage of counts corresponding to S100A9.
  SeuratObject@meta.data$percent.S100A8 <- Seurat::PercentageFeatureSet(
    SeuratObject,
    pattern = patterns[species, "S100A8.pattern"]
  ) # Percentage of counts corresponding to S100A8.
  SeuratObject@meta.data$percent.FCGR3B <- Seurat::PercentageFeatureSet(
    SeuratObject,
    pattern = patterns[species, "FCGR3B.pattern"]
  ) # Percentage of counts corresponding to FCGR3B.
  SeuratObject@meta.data$log10_nFeature_RNA <- log10(
    SeuratObject@meta.data$nFeature_RNA
  )
  SeuratObject@meta.data$log10_nCount_RNA <- log10(
    SeuratObject@meta.data$nCount_RNA
  )
  SeuratObject@meta.data$complexity <- SeuratObject@meta.data$log10_nFeature_RNA /
    SeuratObject@meta.data$log10_nCount_RNA # Complexity, corresponding to the amount of genes that are covered by the counts of each cell.

  # This steps are only if the Seurat object has log-normalized data layers, which are not always present.
  object.layers <- SeuratObject::Layers(
    SeuratObject,
    assay = assay
  )
  if (!is.null(layers)) {
    # Check if the Seurat object has log-normalized data layers, and if not, skip the log-normalized QC calculations.
    if (!all(layers %in% object.layers)) {
      stop(
        "Some specified layers are not present in the Seurat object. Please check the provided layers argument."
      )
    }
    # Also check if the user provided layers are present in the Seurat object.
  } else {
    if (!any(grepl("data", object.layers))) {
      warning(
        "User did not specified any layers and no log-normalized data layers were found in the Seurat object. Skipping log-normalized QC calculations."
      )
      return(SeuratObject)
    }
    layers <- object.layers[grepl("data", object.layers)]
  }

  DataLayersQC <- function(
    SeuratObject,
    assay = "RNA",
    layers = NULL,
    perform.cell.cycle.scoring = TRUE,
    perform.MALAT1.test = TRUE,
    ...
  ) {
    # Helper function to calculate QC metrics for each log-normalized data layer in the Seurat object.
    # We calculate the total number of log-normalized counts per cell, and the total
    # number of features detected per cell in log counts, and cell cycle scoring.

    # 1. Fetch all log-normalized 'data' layers names.
    if (is.null(layers)) {
      layers <- SeuratObject::Layers(
        SeuratObject,
        assay = assay,
        search = "data"
      )
    }

    # 2. Cycle per layer and calculate the total number of log-normalized counts per cell,
    # and the total number of features detected per cell in log counts, and cell cycle scoring.
    # Data will be added to the Seurat object's metadata per layer.
    log.metadata <- lapply(layers, function(layer) {
      # Extract the layer data to compute QC metrics.
      layer.data <- SeuratObject::LayerData(
        SeuratObject,
        assay = assay,
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
    # Remve the cell.id column from log.metadata as it is now redundant with the row names.
    log.metadata$cell.id <- NULL
    # Add the calculated QC metrics to the Seurat object's metadata.
    SeuratObject <- Seurat::AddMetaData(
      object = SeuratObject,
      metadata = log.metadata
    )

    # Add cell cycle scoring results to the Seurat object's metadata if requested.
    # Perform cell cycle scoring if requested.
    if (perform.cell.cycle.scoring) {
      SeuratObject <- suppressWarnings(Seurat::CellCycleScoring(
        SeuratObject,
        # Uses Seurat S and G2M genes dataset.
        s.features = cc.genes.updated.2019$s.genes,
        g2m.features = cc.genes.updated.2019$g2m.genes
      ))

      # Apply the MALAT1 thresholding function to the Seurat object upon user request
      if (perform.MALAT1.test) {
        SeuratObject <- .CalculateFeatureThresholdSeurat(
          SeuratObject = SeuratObject,
          assay = assay,
          layers = layers,
          feature = patterns[species, "MALAT1.pattern"],
          ...
        )
      }
    }

    return(SeuratObject)
  }

  # Apply the helper function to the Seurat object.
  SeuratObject <- DataLayersQC(
    SeuratObject,
    assay = assay,
    layers = layers,
    perform.cell.cycle.scoring = perform.cell.cycle.scoring,
    perform.MALAT1.test = perform.MALAT1.test,
    ...
  )

  return(SeuratObject)
}
