# Tests comparing R implementation to Python implementation
# These tests require reticulate and a Python environment with dPCA installed

skip_if_no_python <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    skip("reticulate not available")
  }

  tryCatch({
    reticulate::py_config()
  }, error = function(e) {
    skip("Python not available")
  })

  # Check if dPCA is installed
  tryCatch({
    reticulate::py_run_string("import dPCA")
  }, error = function(e) {
    skip("Python dPCA package not installed")
  })
}


test_that("R marginalization matches Python marginalization", {
  skip_if_no_python()

  reticulate::py_run_string("
import numpy as np
from dPCA import dPCA

np.random.seed(42)
n_neurons = 10
n_time = 5
n_stim = 3

X_py = np.random.randn(n_neurons, n_time, n_stim)

dpca_py = dPCA.dPCA(labels='ts', n_components=3)
mXs_py = dpca_py._marginalize(X_py)
  ")

  # Get Python data
  X_py <- reticulate::py$X_py
  mXs_py <- reticulate::py$mXs_py

  # Convert to R array
  X_r <- array(X_py, dim = dim(X_py))

  # Compute R marginalizations
  model <- dPCA(labels = "ts", n_components = 3)
  mXs_r <- dpca_marginalize(model, X_r)

  # Compare marginalizations
  for (key in names(mXs_r)) {
    py_mX <- mXs_py[[key]]
    r_mX <- mXs_r[[key]]

    # Allow for small numerical differences
    expect_equal(as.vector(r_mX), as.vector(py_mX), tolerance = 1e-10,
                 label = paste("marginalization", key))
  }
})


test_that("R fit matches Python fit", {
  skip_if_no_python()

  reticulate::py_run_string("
import numpy as np
from dPCA import dPCA

np.random.seed(123)
n_neurons = 15
n_time = 8
n_stim = 4

X_py = np.random.randn(n_neurons, n_time, n_stim)

dpca_py = dPCA.dPCA(labels='ts', n_components=5, regularizer=0)
dpca_py.fit(X_py)
  ")

  # Get Python results
  X_py <- reticulate::py$X_py
  dpca_py <- reticulate::py$dpca_py
  P_py <- dpca_py$P
  D_py <- dpca_py$D

  # R fit
  X_r <- array(X_py, dim = dim(X_py))
  model <- dPCA(labels = "ts", n_components = 5, regularizer = 0)

  # Need to set same random seed behavior - use deterministic approach
  # Since randomized SVD uses random initialization, we compare structure not exact values
  model <- dpca_fit(model, X_r)

  # Check shapes match
  for (key in names(model$D)) {
    expect_equal(dim(model$D[[key]]), dim(D_py[[key]]),
                 label = paste("D shape", key))
    expect_equal(dim(model$P[[key]]), dim(P_py[[key]]),
                 label = paste("P shape", key))
  }

  # The columns of D and P are eigenvectors - they may differ by sign
  # But D^T @ P should be approximately identity-like for orthonormal bases
  for (key in names(model$D)) {
    r_prod <- t(model$D[[key]]) %*% model$P[[key]]
    py_prod <- t(D_py[[key]]) %*% P_py[[key]]

    # Both should be approximately identity (or close to it)
    expect_equal(abs(diag(r_prod)), abs(diag(py_prod)), tolerance = 0.1,
                 label = paste("D^T @ P diagonal", key))
  }
})


test_that("R transform output has same shape as Python", {
  skip_if_no_python()

  reticulate::py_run_string("
import numpy as np
from dPCA import dPCA

np.random.seed(456)
n_neurons = 12
n_time = 6
n_stim = 3

X_py = np.random.randn(n_neurons, n_time, n_stim)

dpca_py = dPCA.dPCA(labels='ts', n_components=4)
Z_py = dpca_py.fit_transform(X_py)
  ")

  X_py <- reticulate::py$X_py
  Z_py <- reticulate::py$Z_py

  X_r <- array(X_py, dim = dim(X_py))
  model <- dPCA(labels = "ts", n_components = 4)
  result <- dpca_fit_transform(model, X_r)
  Z_r <- result$Z

  # Check shapes match
  for (key in names(Z_r)) {
    expect_equal(dim(Z_r[[key]]), dim(Z_py[[key]]),
                 label = paste("Z shape", key))
  }
})


test_that("R denoise_mask matches Python denoise_mask", {
  skip_if_no_python()

  reticulate::py_run_string("
import numpy as np
from dPCA.utils import denoise_mask

mask = np.array([0, 1, 0, 0, 1, 1, 0, 1, 0], dtype=np.int32)
result_py = denoise_mask(mask.copy(), 2)
  ")

  mask_py <- reticulate::py$mask
  result_py <- reticulate::py$result_py

  mask_r <- as.integer(mask_py)
  result_r <- denoise_mask(mask_r, 2)

  expect_equal(as.integer(result_r), as.integer(result_py))
})


test_that("R classification matches Python classification", {
  skip_if_no_python()

  reticulate::py_run_string("
import numpy as np
from dPCA.utils import classification

np.random.seed(789)
class_means = np.array([[0., 0.], [1., 1.], [2., 2.]])
test_data = np.array([[0.1, 0.1], [0.9, 0.9], [2.1, 2.1]])
result_py = classification(class_means, test_data)
  ")

  class_means_py <- reticulate::py$class_means
  test_data_py <- reticulate::py$test_data
  result_py <- reticulate::py$result_py

  class_means_r <- matrix(class_means_py, nrow = 3, ncol = 2)
  test_data_r <- matrix(test_data_py, nrow = 3, ncol = 2)
  result_r <- classification(class_means_r, test_data_r)

  expect_equal(result_r, as.vector(result_py), tolerance = 1e-10)
})


test_that("R get_parameter_combinations matches Python", {
  skip_if_no_python()

  reticulate::py_run_string("
from dPCA import dPCA

dpca_py = dPCA.dPCA(labels='tsd', n_components=3)
pcombs_py = dpca_py._get_parameter_combinations(join=False)
  ")

  pcombs_py <- reticulate::py$pcombs_py

  model <- dPCA(labels = "tsd", n_components = 3)
  pcombs_r <- get_parameter_combinations(model, join = FALSE)

  # Check same keys
  expect_equal(sort(names(pcombs_r)), sort(names(pcombs_py)))

  # Check each combination (accounting for Python sets vs R vectors)
  for (key in names(pcombs_r)) {
    py_set <- sort(as.integer(pcombs_py[[key]]))
    r_vec <- sort(pcombs_r[[key]])
    expect_equal(r_vec, py_set, label = paste("combination", key))
  }
})


test_that("R and Python produce consistent explained variance ratios", {
  skip_if_no_python()

  reticulate::py_run_string("
import numpy as np
from dPCA import dPCA

np.random.seed(999)
n_neurons = 20
n_time = 10
n_stim = 5

X_py = np.random.randn(n_neurons, n_time, n_stim)

dpca_py = dPCA.dPCA(labels='ts', n_components=5)
dpca_py.fit(X_py)
Z_py = dpca_py.transform(X_py)
evr_py = dpca_py.explained_variance_ratio_
  ")

  X_py <- reticulate::py$X_py
  evr_py <- reticulate::py$evr_py

  X_r <- array(X_py, dim = dim(X_py))
  model <- dPCA(labels = "ts", n_components = 5)
  model <- dpca_fit(model, X_r)
  Z_r <- dpca_transform(model, X_r)

  # Explained variance ratios should be similar (not exact due to different SVD random states)
  # Check that they're in the same ballpark (same order of magnitude)
  for (key in names(model$explained_variance_ratio_)) {
    r_evr <- model$explained_variance_ratio_[[key]]
    py_evr <- evr_py[[key]]

    # Check that the sum is similar (total explained variance for this marginalization)
    expect_equal(sum(r_evr), sum(py_evr), tolerance = 0.1,
                 label = paste("total explained variance", key))
  }
})


test_that("R reconstruction quality is similar to Python", {
  skip_if_no_python()

  reticulate::py_run_string("
import numpy as np
from dPCA import dPCA

np.random.seed(111)
n_neurons = 15
n_time = 8
n_stim = 4

X_py = np.random.randn(n_neurons, n_time, n_stim)

dpca_py = dPCA.dPCA(labels='ts', n_components=10)
dpca_py.fit(X_py)

# Compute reconstruction error for 't' marginalization
Z_t = dpca_py.transform(X_py, marginalization='t')
X_recon = dpca_py.inverse_transform(Z_t, 't')
recon_error_py = np.sum((X_py - X_recon)**2) / np.sum(X_py**2)
  ")

  X_py <- reticulate::py$X_py
  recon_error_py <- reticulate::py$recon_error_py

  X_r <- array(X_py, dim = dim(X_py))
  model <- dPCA(labels = "ts", n_components = 10)
  model <- dpca_fit(model, X_r)

  Z_t <- dpca_transform(model, X_r, marginalization = "t")
  X_recon <- dpca_inverse_transform(model, Z_t, "t")
  recon_error_r <- sum((X_r - X_recon)^2) / sum(X_r^2)

  # Reconstruction errors should be similar
  expect_equal(recon_error_r, recon_error_py, tolerance = 0.05)
})
