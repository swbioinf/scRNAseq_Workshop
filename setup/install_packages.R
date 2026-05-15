#!/usr/bin/env Rscript
# R package installation for workshop dependencies

# Change library directory to ${HOME}/.library so root is not needed
user_lib <- "~/.library"
.libPaths(user_lib)

# Install BiocManager first, may prevent conflicting dependencies if
# CRAN packages are installed first
install.packages("BiocManager", repos = "https://cran.csiro.au/", lib = user_lib)

bioc_pkgs <- c(
  "BiocGenerics",
  "SingleR",
  "celldex",
  "BiocGenerics",
  "DelayedArray",
  "DelayedMatrixStats",
  "limma",
  "S4Vectors",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "edgeR"
)

BiocManager::install(bioc_pkgs, lib = user_lib, ask = FALSE, update = FALSE)

# Install CRAN packges
cran_pkgs <- c(
  "Seurat",
  "dplyr",
  "remotes",
  "R.utils",
  "harmony",
  "hdf5r",
  "clustree",
  "RColorBrewer",
  "tidyverse",
  "pander",
  "qs2",
  "plotly",
  "ggplotly"
)

install.packages(cran_pkgs, lib = user_lib, repos = "https://cran.csiro.au/")

# Lastly github pkgs
remotes::install_github("immunogenomics/presto")

# Validate installation

cat("========== Installation validation ==========\n")
cat("Target library:", user_lib, "\n\n")

expected_pkgs <- c("BiocManager", bioc_pkgs, cran_pkgs, "presto")

check_package <- function(pkg, lib) {
  installed <- dir.exists(file.path(lib, pkg))
  loadable  <- FALSE
  version   <- NA_character_
  message   <- NULL

  if (installed) {
    loadable <- suppressWarnings(
      requireNamespace(pkg, quietly = TRUE, lib.loc = lib)
    )
    if (loadable) {
      version <- tryCatch(
        as.character(packageVersion(pkg, lib.loc = lib)),
        error = function(e) NA_character_
      )
    } else {
      message <- "installed but fails to load"
    }
  } else {
    message <- "not installed"
  }

  list(
    package   = pkg,
    installed = installed,
    loadable  = loadable,
    version   = version,
    status    = if (installed && loadable) "OK" else message
  )
}

results <- lapply(expected_pkgs, check_package, lib = user_lib)

# Print results table
cat(sprintf("%-35s %-10s %s\n", "Package", "Version", "Status"))
cat(strrep("-", 60), "\n")
for (r in results) {
  cat(sprintf(
    "%-35s %-10s %s\n",
    r$package,
    if (is.na(r$version)) "-" else r$version,
    r$status
  ))
}

cat("\n")

failed <- Filter(function(r) !identical(r$status, "OK"), results)

if (length(failed) == 0) {
  cat("All", length(expected_pkgs), "packages installed and loadable.\n")
} else {
  cat(sprintf("FAILED: %d package(s) did not pass validation:\n", length(failed)))
  for (r in failed) {
    cat(sprintf("  - %s: %s\n", r$package, r$status))
  }
  quit(status = 1)
}
