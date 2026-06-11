#' Calculate the cellular detection rate
#'
#' Given a feature X cell table of raw counts, calculate the scaled cellular detection
#' rate for each cell. The cellular detection rate of a cell is a number between
#' 0 and 1 representing the fraction of the total number of features measured
#' that have a count > 0 in that cell. It is useful as a covariate in
#' differential expression analyses. This function returns a vector of scaled CDRs.
#'
#' @param SeuratObject A Seurat object.
#' @return The Seurat object with an additional "CDR" metadata slot with the calculated cdr.
#' @examples
#' \dontrun{
#' SeuratObject <- CalculateCDR(SeuratObject)
#' }
#' @export
#'
CalculateCDR <- function(SeuratObject) {
  mtry <- try(
    SeuratObject$CDR <- scale(colMeans(
      as.matrix(SeuratObject@assays$RNA@layers$counts) > 0
    )),
    silent = TRUE
  )

  if (inherits(mtry, "try-error")) {
    SeuratObject$CDR <- scale(colMeans(
      as.matrix(SeuratObject@assays$RNA@counts) > 0
    ))
  }

  return(SeuratObject)
}

# Credit:
# https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_MASTcpmDetRate.R
