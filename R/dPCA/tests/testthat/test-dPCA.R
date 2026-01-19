# Tests for dPCA core functionality

test_that("dPCA object is created correctly", {
  # Test with string labels
  model <- dPCA(labels = "ts", n_components = 5)
  expect_s3_class(model, "dPCA")
  expect_equal(model$labels, "ts")
  expect_equal(model$n_components, 5)
  expect_equal(names(model$marginalizations), c("t", "s", "ts"))

  # Test with integer labels
  model <- dPCA(labels = 3, n_components = 10)
  expect_equal(model$labels, "abc")
  expect_equal(length(model$marginalizations), 7)  # a, b, c, ab, ac, bc, abc

  # Test with join
  model <- dPCA(labels = "ts", join = list(ts = c("s", "ts")), n_components = 5)
  expect_true("ts" %in% names(model$marginalizations))
  expect_false("s" %in% names(model$marginalizations))
})


test_that("get_parameter_combinations generates correct subsets", {
  model <- dPCA(labels = "xyz", n_components = 5)

  # Without join
  pcombs <- get_parameter_combinations(model, join = FALSE)
  expect_equal(names(pcombs), c("x", "y", "z", "xy", "xz", "yz", "xyz"))

  # Check indices (0-indexed like Python)
  expect_equal(pcombs[["x"]], 0)
  expect_equal(pcombs[["y"]], 1)
  expect_equal(pcombs[["xy"]], c(0, 1))
  expect_equal(pcombs[["xyz"]], c(0, 1, 2))
})


test_that("dpca_marginalize computes correct marginalizations", {
  # Create simple test data
  set.seed(123)
  n_neurons <- 10
  n_time <- 5
  n_stim <- 3

  # X has shape (n_neurons, n_time, n_stim)
  X <- array(rnorm(n_neurons * n_time * n_stim), dim = c(n_neurons, n_time, n_stim))

  model <- dPCA(labels = "ts", n_components = 3)
  mXs <- dpca_marginalize(model, X)

  # Check that marginalizations were computed
  expect_true(all(c("t", "s", "ts") %in% names(mXs)))

  # Check dimensions (all should be n_neurons x (n_time * n_stim))
  expect_equal(dim(mXs[["t"]]), c(n_neurons, n_time * n_stim))
  expect_equal(dim(mXs[["s"]]), c(n_neurons, n_time * n_stim))
  expect_equal(dim(mXs[["ts"]]), c(n_neurons, n_time * n_stim))

  # Sum of marginalizations should approximately equal centered data
  X_centered <- X - array(rowMeans(matrix(X, nrow = n_neurons)), dim = dim(X))
  X_flat <- matrix(X_centered, nrow = n_neurons)
  mX_sum <- mXs[["t"]] + mXs[["s"]] + mXs[["ts"]]
  expect_equal(mX_sum, X_flat, tolerance = 1e-10)
})


test_that("dpca_fit and dpca_transform work correctly", {
  set.seed(456)
  n_neurons <- 20
  n_time <- 10
  n_stim <- 4

  X <- array(rnorm(n_neurons * n_time * n_stim), dim = c(n_neurons, n_time, n_stim))

  model <- dPCA(labels = "ts", n_components = 5)
  model <- dpca_fit(model, X)

  # Check that P and D matrices are computed
  expect_true(!is.null(model$P))
  expect_true(!is.null(model$D))
  expect_true(all(c("t", "s", "ts") %in% names(model$P)))
  expect_true(all(c("t", "s", "ts") %in% names(model$D)))

  # Check dimensions
  expect_equal(nrow(model$P[["t"]]), n_neurons)
  expect_equal(ncol(model$P[["t"]]), 5)
  expect_equal(nrow(model$D[["t"]]), n_neurons)
  expect_equal(ncol(model$D[["t"]]), 5)

  # Transform
  Z <- dpca_transform(model, X)
  expect_true(all(c("t", "s", "ts") %in% names(Z)))
  expect_equal(dim(Z[["t"]]), c(5, n_time, n_stim))

  # Transform with specific marginalization
  Z_t <- dpca_transform(model, X, marginalization = "t")
  expect_equal(dim(Z_t), c(5, n_time, n_stim))
})


test_that("dpca_fit_transform works correctly", {
  set.seed(789)
  n_neurons <- 15
  n_time <- 8
  n_stim <- 3

  X <- array(rnorm(n_neurons * n_time * n_stim), dim = c(n_neurons, n_time, n_stim))

  model <- dPCA(labels = "ts", n_components = 4)
  result <- dpca_fit_transform(model, X)

  expect_true(!is.null(result$model$P))
  expect_true(!is.null(result$Z))
  expect_true(all(c("t", "s", "ts") %in% names(result$Z)))
})


test_that("dpca_inverse_transform and dpca_reconstruct work", {
  set.seed(101)
  n_neurons <- 12
  n_time <- 6
  n_stim <- 4

  X <- array(rnorm(n_neurons * n_time * n_stim), dim = c(n_neurons, n_time, n_stim))

  model <- dPCA(labels = "ts", n_components = 5)
  model <- dpca_fit(model, X)

  # Transform and inverse transform
  Z_t <- dpca_transform(model, X, marginalization = "t")
  X_recon <- dpca_inverse_transform(model, Z_t, "t")

  expect_equal(dim(X_recon), dim(X))

  # Reconstruct should be equivalent
  X_recon2 <- dpca_reconstruct(model, X, "t")
  expect_equal(dim(X_recon2), dim(X))
})


test_that("regularization is applied correctly", {
  set.seed(202)
  n_neurons <- 15
  n_time <- 8
  n_stim <- 3

  X <- array(rnorm(n_neurons * n_time * n_stim), dim = c(n_neurons, n_time, n_stim))

  # Without regularization
  model1 <- dPCA(labels = "ts", n_components = 4, regularizer = 0)
  model1 <- dpca_fit(model1, X)

  # With regularization
  model2 <- dPCA(labels = "ts", n_components = 4, regularizer = 0.01)
  model2 <- dpca_fit(model2, X)

  # Results should be different
  expect_false(all(model1$D[["t"]] == model2$D[["t"]]))
})


test_that("dPCA works with different number of components per marginalization", {
  set.seed(303)
  n_neurons <- 20
  n_time <- 10
  n_stim <- 5

  X <- array(rnorm(n_neurons * n_time * n_stim), dim = c(n_neurons, n_time, n_stim))

  model <- dPCA(labels = "ts", n_components = list(t = 3, s = 2, ts = 4))
  model <- dpca_fit(model, X)

  expect_equal(ncol(model$D[["t"]]), 3)
  expect_equal(ncol(model$D[["s"]]), 2)
  expect_equal(ncol(model$D[["ts"]]), 4)

  Z <- dpca_transform(model, X)
  expect_equal(dim(Z[["t"]])[1], 3)
  expect_equal(dim(Z[["s"]])[1], 2)
  expect_equal(dim(Z[["ts"]])[1], 4)
})


test_that("dPCA works with 3 parameter dimensions", {
  set.seed(404)
  n_neurons <- 15
  n_time <- 5
  n_stim <- 3
  n_dec <- 2

  # X has shape (n_neurons, n_time, n_stim, n_dec)
  X <- array(rnorm(n_neurons * n_time * n_stim * n_dec),
             dim = c(n_neurons, n_time, n_stim, n_dec))

  model <- dPCA(labels = "tsd", n_components = 3)
  model <- dpca_fit(model, X)

  # Should have 7 marginalizations: t, s, d, ts, td, sd, tsd
  expect_equal(length(model$D), 7)
  expect_true(all(c("t", "s", "d", "ts", "td", "sd", "tsd") %in% names(model$D)))

  Z <- dpca_transform(model, X)
  expect_equal(dim(Z[["t"]]), c(3, n_time, n_stim, n_dec))
  expect_equal(dim(Z[["tsd"]]), c(3, n_time, n_stim, n_dec))
})
