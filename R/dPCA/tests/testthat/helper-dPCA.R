# Helper file for testthat tests
# This file is automatically sourced before tests run

# Source the package files directly for development testing
# This allows tests to run even before the package is formally installed

# Find the package root directory
find_pkg_root <- function() {
  # Start from the current test directory
  test_dir <- getwd()

  # Go up until we find DESCRIPTION
  check_dir <- test_dir
  for (i in 1:10) {
    if (file.exists(file.path(check_dir, "DESCRIPTION"))) {
      return(check_dir)
    }
    check_dir <- dirname(check_dir)
  }

  # Try relative to testthat directory
  possible_roots <- c(
    file.path(test_dir, "..", ".."),
    file.path(test_dir, "..")
  )

  for (root in possible_roots) {
    if (file.exists(file.path(root, "DESCRIPTION"))) {
      return(normalizePath(root))
    }
  }

  NULL
}

# Source R files if package is not installed
if (!requireNamespace("dPCA", quietly = TRUE)) {
  pkg_root <- find_pkg_root()
  if (!is.null(pkg_root)) {
    r_dir <- file.path(pkg_root, "R")
    if (dir.exists(r_dir)) {
      source(file.path(r_dir, "utils.R"))
      source(file.path(r_dir, "dPCA.R"))
    }
  }
}
