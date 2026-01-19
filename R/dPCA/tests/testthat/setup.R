# Setup file for testthat tests
# This file is automatically sourced before tests run

# Source the package files directly for development testing
# In installed package, these would come from the NAMESPACE

pkg_dir <- system.file(package = "dPCA")
if (pkg_dir == "" || !dir.exists(file.path(pkg_dir, "R"))) {
  # Package not installed, source files directly
  pkg_root <- normalizePath(file.path(dirname(getwd()), ".."))
  if (file.exists(file.path(pkg_root, "R", "utils.R"))) {
    source(file.path(pkg_root, "R", "utils.R"))
    source(file.path(pkg_root, "R", "dPCA.R"))
  }
}
