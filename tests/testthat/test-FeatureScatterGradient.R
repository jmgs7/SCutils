library(testthat)
library(Seurat)

make_test_seurat <- function() {
  counts <- Matrix::Matrix(
    c(
      10,
      8,
      6,
      4,
      2,
      1,
      1,
      2,
      3,
      4,
      5,
      6,
      0,
      1,
      0,
      1,
      0,
      1,
      6,
      5,
      4,
      3,
      2,
      1
    ),
    nrow = 4,
    byrow = TRUE,
    sparse = TRUE,
    dimnames = list(
      c("GeneA", "GeneB", "GeneC", "GeneD"),
      paste0("cell", 1:6)
    )
  )

  object <- CreateSeuratObject(counts = counts)
  object <- NormalizeData(object, verbose = FALSE)
  object$percent.mt <- c(1, 2, 3, 4, 5, 6)
  object$batch <- c("A", "A", "B", "B", "C", NA)
  object$singleton_group <- c("s1", "s2", "s2", "s3", "s4", "s5")
  Idents(object) <- c("g1", "g1", "g2", "g2", "g3", "g3")

  object
}

extract_patchwork_text <- function(plot) {
  grob <- patchwork::patchworkGrob(plot)
  labels <- character()

  walk_grob <- function(x) {
    if (inherits(x, "text")) {
      labels <<- c(labels, x$label)
    }

    if (!is.null(x$grobs)) {
      for (i in seq_along(x$grobs)) {
        walk_grob(x$grobs[[i]])
      }
    }

    if (!is.null(x$children)) {
      for (child_name in names(x$children)) {
        walk_grob(x$children[[child_name]])
      }
    }
  }

  walk_grob(grob)
  unique(labels)
}

test_that("FeatureScatterGradient validates scalar inputs", {
  object <- make_test_seurat()

  expect_error(
    FeatureScatterGradient(
      "not_seurat",
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mt"
    ),
    regexp = "'SeuratObject' must be a Seurat object.",
    fixed = TRUE
  )

  expect_error(
    FeatureScatterGradient(object, character(), "nFeature_RNA", "percent.mt"),
    regexp = "'feature1' must be a single non-NA character value.",
    fixed = TRUE
  )

  expect_error(
    FeatureScatterGradient(object, "nCount_RNA", NA_character_, "percent.mt"),
    regexp = "'feature2' must be a single non-NA character value.",
    fixed = TRUE
  )

  expect_error(
    FeatureScatterGradient(object, "nCount_RNA", "nFeature_RNA", NA_character_),
    regexp = "'gradient' must be a single non-NA character value.",
    fixed = TRUE
  )

  expect_error(
    FeatureScatterGradient(
      object,
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mt",
      group.by = NA_character_
    ),
    regexp = "'group.by' must be NULL or a single non-NA character value.",
    fixed = TRUE
  )

  expect_error(
    FeatureScatterGradient(
      object,
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mt",
      pt.size = -1
    ),
    regexp = "'pt.size' must be a single non-negative numeric value.",
    fixed = TRUE
  )

  expect_error(
    FeatureScatterGradient(
      object,
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mt",
      lower.limit = "a"
    ),
    regexp = "'lower.limit' must be NULL or a single non-NA numeric value.",
    fixed = TRUE
  )

  expect_error(
    FeatureScatterGradient(
      object,
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mt",
      upper.limit = 2
    ),
    regexp = "'lower.limit' must be specified when 'upper.limit' is non-NULL.",
    fixed = TRUE
  )

  expect_error(
    FeatureScatterGradient(
      object,
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mt",
      lower.limit = 2,
      upper.limit = 1
    ),
    regexp = "'upper.limit' must be strictly greater than 'lower.limit'.",
    fixed = TRUE
  )

  expect_error(
    FeatureScatterGradient(
      object,
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mt",
      corr.method = "bad"
    ),
    regexp = "'corr.method' must be one of 'pearson', 'spearman', 'kendall'.",
    fixed = TRUE
  )

  expect_error(
    FeatureScatterGradient(
      object,
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mt",
      scale.colors = "bad"
    ),
    regexp = "'scale.colors' must be one of",
    fixed = TRUE
  )

  expect_error(
    FeatureScatterGradient(
      object,
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mt",
      layer1 = NA_character_
    ),
    regexp = "'layer1' must be NULL or a single non-NA character value.",
    fixed = TRUE
  )

  expect_error(
    FeatureScatterGradient(
      object,
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mt",
      layer2 = character()
    ),
    regexp = "'layer2' must be NULL or a single non-NA character value.",
    fixed = TRUE
  )

  expect_error(
    FeatureScatterGradient(
      object,
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mt",
      plot.title = c("a", "b")
    ),
    regexp = "'plot.title' must be NULL or a single non-NA character value.",
    fixed = TRUE
  )
})

test_that("FeatureScatterGradient validates data access failures clearly", {
  object <- make_test_seurat()

  expect_error(
    FeatureScatterGradient(
      object,
      "missing_feature",
      "nFeature_RNA",
      "percent.mt",
      group.by = NULL
    ),
    regexp = "Cannot fetch feature 'missing_feature'",
    fixed = TRUE
  )

  expect_error(
    FeatureScatterGradient(
      object,
      "nCount_RNA",
      "nFeature_RNA",
      "missing_gradient",
      group.by = NULL
    ),
    regexp = "Cannot fetch feature 'missing_gradient'",
    fixed = TRUE
  )

  expect_error(
    FeatureScatterGradient(
      object,
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mt",
      group.by = "missing_group"
    ),
    regexp = "Cannot fetch grouping variable 'missing_group'",
    fixed = TRUE
  )
})

test_that("FeatureScatterGradient returns a ggplot for ungrouped metadata features", {
  object <- make_test_seurat()

  plot <- FeatureScatterGradient(
    SeuratObject = object,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    gradient = "percent.mt",
    group.by = NULL
  )

  built <- ggplot2::ggplot_build(plot)

  expect_s3_class(plot, "ggplot")
  expect_identical(plot$labels$x, "nCount_RNA")
  expect_identical(plot$labels$y, "nFeature_RNA")
  expect_null(plot$labels$color)
  expect_match(plot$labels$title, "^Pearson =")
  expect_equal(nrow(built$data[[1]]), ncol(object))
})

test_that("FeatureScatterGradient respects layer-aware axis labels and sentinel layers", {
  object <- make_test_seurat()

  plot <- FeatureScatterGradient(
    SeuratObject = object,
    feature1 = "GeneA",
    feature2 = "GeneB",
    gradient = "percent.mt",
    group.by = NULL,
    layer1 = "counts",
    layer2 = "null"
  )

  built <- ggplot2::ggplot_build(plot)

  expect_identical(plot$labels$x, "GeneA_counts")
  expect_identical(plot$labels$y, "GeneB")
  expect_equal(nrow(built$data[[1]]), ncol(object))
})

test_that("FeatureScatterGradient supports alternative correlation methods", {
  object <- make_test_seurat()

  spearman_plot <- FeatureScatterGradient(
    SeuratObject = object,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    gradient = "percent.mt",
    group.by = NULL,
    corr.method = "spearman"
  )

  kendall_plot <- FeatureScatterGradient(
    SeuratObject = object,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    gradient = "percent.mt",
    group.by = NULL,
    corr.method = "kendall"
  )

  expect_match(spearman_plot$labels$title, "^Spearman =")
  expect_match(kendall_plot$labels$title, "^Kendall =")
})

test_that("FeatureScatterGradient applies explicit gradient limits", {
  object <- make_test_seurat()

  plot <- FeatureScatterGradient(
    SeuratObject = object,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    gradient = "percent.mt",
    group.by = NULL,
    lower.limit = 0,
    upper.limit = 10
  )

  colour_scale <- plot$scales$get_scales("colour")
  expect_equal(colour_scale$limits, c(0, 10))
})

test_that("FeatureScatterGradient returns a patchwork for grouped ident plots", {
  object <- make_test_seurat()

  plot <- FeatureScatterGradient(
    SeuratObject = object,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    gradient = "percent.mt",
    group.by = "ident"
  )

  expect_s3_class(plot, "patchwork")
  expect_equal(length(plot$patches$plots) + 1L, length(levels(Idents(object))))
  expect_true("nCount_RNA VS nFeature_RNA" %in% extract_patchwork_text(plot))
  expect_match(plot$labels$title, "^Pearson =")

  built <- ggplot2::ggplot_build(plot$patches$plots[[1]])
  expect_gt(nrow(built$data[[1]]), 0)
})

test_that("FeatureScatterGradient uses custom grouped main title", {
  object <- make_test_seurat()

  plot <- FeatureScatterGradient(
    SeuratObject = object,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    gradient = "percent.mt",
    group.by = "ident",
    plot.title = "Custom grouped title"
  )

  expect_true("Custom grouped title" %in% extract_patchwork_text(plot))
  expect_identical(plot$patches$annotation$theme$plot.title$face, "bold")
})

test_that("FeatureScatterGradient drops missing grouping values before plotting", {
  object <- make_test_seurat()

  plot <- FeatureScatterGradient(
    SeuratObject = object,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    gradient = "percent.mt",
    group.by = "batch"
  )

  expect_equal(length(plot$patches$plots) + 1L, 3)
  expect_true("nCount_RNA VS nFeature_RNA" %in% extract_patchwork_text(plot))
})

test_that("FeatureScatterGradient handles groups with insufficient correlation data", {
  object <- make_test_seurat()

  plot <- FeatureScatterGradient(
    SeuratObject = object,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    gradient = "percent.mt",
    group.by = "singleton_group"
  )

  panel_titles <- c(
    plot$labels$title,
    vapply(
      plot$patches$plots,
      function(p) p$labels$title,
      character(1)
    )
  )
  expect_true(any(grepl("NA$", panel_titles)))
})

test_that("FeatureScatterGradient errors when all grouping values are missing", {
  object <- make_test_seurat()
  object$all_missing <- rep(NA_character_, ncol(object))

  expect_error(
    FeatureScatterGradient(
      object,
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mt",
      group.by = "all_missing"
    ),
    regexp = "No groups available to plot after filtering missing grouping values.",
    fixed = TRUE
  )
})
