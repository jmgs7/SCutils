#!/usr/bin/env Rscript

# Sync feature threshold functions source files from the vendored/submodule repository into
# the top-level SCutils R/ directory so they are part of the package source.

# Define the project-relative root paths.
project.root <- "."
vendor.root <- file.path(project.root, "vendor", "feature_threshold")
target.root <- file.path(project.root, "R")

# Define the source -> target file mapping explicitly.
# Left side: file in the MALAT1 repo.
# Right side: destination file in SCutils/R/.

source.files <- c(
  file.path(vendor.root, "R", "ComputeFeatureThreshold.R"),
  file.path(vendor.root, "R", "ComputeFeatureThresholdSeurat.R")
)
file.map <- c(
  file.path(target.root, "ComputeFeatureThreshold.R"),
  file.path(target.root, "ComputeFeatureThresholdSeurat.R")
)
names(file.map) <- source.files

# Check that the vendored repo exists.
if (!dir.exists(vendor.root)) {
  stop(
    "Vendor directory not found: ",
    normalizePath(vendor.root, mustWork = FALSE),
    "\n",
    "Did you initialize the git submodule?"
  )
}

# Check that all source files exist before copying anything.
missing.files <- source.files[!file.exists(source.files)]

if (length(missing.files) > 0L) {
  stop(
    "The following source files are missing in the MALAT1 repository:\n- ",
    paste(missing.files, collapse = "\n- ")
  )
}

# Create target directories if they do not already exist.
target.dirs <- unique(dirname(unname(file.map)))

for (target.dir in target.dirs) {
  if (!dir.exists(target.dir)) {
    ok <- dir.create(target.dir, recursive = TRUE, showWarnings = FALSE)
    if (!ok) {
      stop("Failed to create target directory: ", target.dir)
    }
  }
}

# Copy the files one by one so failures are easy to diagnose.
copy.ok <- vapply(
  X = seq_along(file.map),
  FUN = function(i) {
    from.file <- names(file.map)[[i]]
    to.file <- unname(file.map)[[i]]

    file.copy(
      from = from.file,
      to = to.file,
      overwrite = TRUE
    )
  },
  FUN.VALUE = logical(1)
)

# Abort if any copy failed.
if (!all(copy.ok)) {
  failed.files <- names(file.map)[!copy.ok]
  stop(
    "Failed to copy the following files:\n- ",
    paste(failed.files, collapse = "\n- ")
  )
}

# Print a compact summary for the terminal/log.
message("Synced MALAT1 files into SCutils:")
for (i in seq_along(file.map)) {
  message("- ", names(file.map)[[i]], " -> ", unname(file.map)[[i]])
}
