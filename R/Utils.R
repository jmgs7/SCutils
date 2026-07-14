ExtractFeatureTestResults <- function(SeuratObject, batch.col, feature = "MALAT1") {
  threshold.col <- paste0(feature, ".threshold")
  pass.col <- paste0(feature, ".pass")
  meta.data.cols <- colnames(SeuratObject@meta.data)

  if (!threshold.col %in% meta.data.cols) {
    stop(threshold.col, " column not found in meta.data")
  } else if (!pass.col %in% meta.data.cols) {
    stop(pass.col, " column not found in meta.data")
  } else if (!batch.col %in% meta.data.cols) {
    stop(batch.col, " column not found in meta.data")
  }

  # Create a data frame with all the information of the test.
  feature.test <- data.frame(
    batch = SeuratObject@meta.data[[batch.col]],
    feature.threshold = SeuratObject@meta.data[[threshold.col]],
    feature.pass = SeuratObject@meta.data[[pass.col]]
  ) |> # Get the threshold per library and compute the pass rate.
    dplyr::group_by(batch) |>
    dplyr::summarise(
      feature.threshold = dplyr::first(feature.threshold),
      nCells.total = dplyr::n(),
      nCells.pass = sum(feature.pass, na.rm = TRUE),
      pass.rate = mean(feature.pass, na.rm = TRUE),
      .groups = "drop"
    ) |>
    as.data.frame()

  # Keep the same order as the SeuratObject.
  row.names(feature.test) <- feature.test$batch
  feature.test <- feature.test[unique(SeuratObject@meta.data[[batch.col]]), ]

  # Keep the same names as the SeuratObject.
  names(feature.test) <- c(
    batch.col,
    threshold.col,
    "nCells.total",
    "nCells.pass",
    "pass.rate"
  )

  return(feature.test)
}
