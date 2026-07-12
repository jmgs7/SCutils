#' @title scGSEAmarkers
#' @description Calculates GSEA for all the clusters in a SeruatObject. Given the output from Seurat::FindAllMarkers(), calculate the GSEA for all clusters given a GSEA database. For example, use MSigDB C8 dataset to calculate enrichment for cell identity markers. This is intended to be a helper for manual annotation of the clusters.
#' @param cluster_markers A data.frame containing the DE markers for multiple clusters
#' from a SC experiment (i.e. a SeuratObject). Expects the output from Seurat::FindAllMarkers.
#' @param reference_markers A GSEA database loaded for fgsea input (see fgsea documentation).
#' @param padj.threshold Threshold for p-adj. Default is 1e-6.
#' @param only.pos Returns only positively enriched markers (NES > 0). Default is true.
#' @param workers Number of cores for parallel computing. Default is 4.
#' @return A list containing the GSEA for each cluster.
#' @examples
#' \dontrun{
#' results <- scGSEAmarkers(cluster_markers, reference_markers)
#' }
#' @import dplyr
#' @import fgsea
#' @export

scGSEAmarkers <- function(
  cluster_markers,
  reference_markers,
  padj.threshold = 1e-6,
  only.pos = TRUE,
  workers = 4
) {
  result <- list()
  cluster_list <- as.character(unique(cluster_markers$cluster))

  for (current_cluster in cluster_list) {
    cluster_stats <- cluster_markers |>
      filter(cluster == current_cluster) |>
      arrange(-avg_log2FC) |>
      select(c(gene, avg_log2FC)) |>
      tibble::deframe()

    res <- fgsea(
      pathways = reference_markers,
      stats = cluster_stats,
      BPPARAM = BiocParallel::MulticoreParam(workers = workers)
    ) |>
      filter(padj < padj.threshold) |>
      arrange(-NES)
    if (only.pos) {
      res <- res |> filter(NES > 0) |> arrange(-NES)
    }

    result[[current_cluster]] <- res
  }

  return(result)
}
