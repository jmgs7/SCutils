#' CalculateQC
#'
#' Estimate common single-cell QC metrics and append them to a Seurat object's metadata.
#'
#' This function calculates percent-based metrics using Seurat::PercentageFeatureSet for
#' mitochondrial, ribosomal, hemoglobin, immunoglobulin, platelet-associated genes and
#' several individual marker genes, and also computes log10-transformed feature/count
#' values and a complexity ratio.
#'
#' @param SeuratObject A Seurat object (v3/v4). The function reads counts from the active assay and adds metadata columns.
#' @details The following metadata columns are added to SeuratObject:
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
#' Patterns are interpreted as regular expressions by Seurat::PercentageFeatureSet.
#' @return The input Seurat object with the new metadata columns added.
#' @examples
#' # CalculateQC(SeuratObject)
#' @import Seurat
#' @export
CalculateQC <- function(SeuratObject) {
  # Estimation of metrics
  SeuratObject[["percent.mt"]] <- PercentageFeatureSet(
    SeuratObject,
    pattern = "^MT-"
  ) # Percentage of counts corresponding to mitochondrial genes.
  SeuratObject[["percent.ribo"]] <- PercentageFeatureSet(
    SeuratObject,
    "^RP[SL]"
  ) # Percentage of counts corresponding to ribosomal genes.
  SeuratObject[["percent.hb"]] <- PercentageFeatureSet(
    SeuratObject,
    "^HB[^(P)]"
  ) # Percentage of counts corresponding to hemoglobin.
  SeuratObject[["percent.ig"]] <- PercentageFeatureSet(SeuratObject, "^IG") # Percentage of counts corresponding to immunoglobulins.
  SeuratObject[["percent.plat"]] <- PercentageFeatureSet(
    SeuratObject,
    "PECAM1|PF4"
  ) # Percentage of counts corresponding to genes associated with platelets.
  SeuratObject[["percent.MALAT1"]] <- PercentageFeatureSet(
    SeuratObject,
    pattern = "MALAT1"
  ) # Percentage of counts corresponding to MALAT1.
  SeuratObject[["percent.S100A9"]] <- PercentageFeatureSet(
    SeuratObject,
    pattern = "S100A9"
  ) # Percentage of counts corresponding to S100A9.
  SeuratObject[["percent.S100A8"]] <- PercentageFeatureSet(
    SeuratObject,
    pattern = "S100A8"
  ) # Percentage of counts corresponding to S100A8.
  SeuratObject[["percent.FCGR3B"]] <- PercentageFeatureSet(
    SeuratObject,
    pattern = "FCGR3B"
  ) # Percentage of counts corresponding to FCGR3B.
  SeuratObject[["log10_nFeature_RNA"]] <- log10(SeuratObject[["nFeature_RNA"]])
  SeuratObject[["log10_nCount_RNA"]] <- log10(SeuratObject[["nCount_RNA"]])
  SeuratObject[["complexity"]] <- SeuratObject[["log10_nFeature_RNA"]] /
    SeuratObject[["log10_nCount_RNA"]] # Complexity, corresponding to the amount of genes that are covered by the counts of each cell.
  return(SeuratObject)
}
