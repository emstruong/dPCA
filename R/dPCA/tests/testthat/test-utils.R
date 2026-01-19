# Tests for utility functions

test_that("denoise_mask removes isolated significant points", {
  # Single isolated point should be removed with n_consecutive = 2
  mask <- c(0, 1, 0, 0, 1, 1, 0, 1, 0)
  result <- denoise_mask(mask, n_consecutive = 2)
  expected <- c(0, 0, 0, 0, 1, 1, 0, 0, 0)
  expect_equal(result, expected)

  # With n_consecutive = 1, nothing should change
  mask <- c(0, 1, 0, 0, 1, 1, 0, 1, 0)
  result <- denoise_mask(mask, n_consecutive = 1)
  expect_equal(result, mask)

  # All ones should stay
  mask <- c(1, 1, 1, 1, 1)
  result <- denoise_mask(mask, n_consecutive = 3)
  expect_equal(result, mask)

  # Two consecutive with n_consecutive = 3 should be removed

  mask <- c(0, 1, 1, 0, 1, 1, 1, 0)
  result <- denoise_mask(mask, n_consecutive = 3)
  expected <- c(0, 0, 0, 0, 1, 1, 1, 0)
  expect_equal(result, expected)
})


test_that("classification computes nearest neighbor accuracy", {
  # Perfect classification
  class_means <- matrix(c(0, 1, 2, 0, 1, 2), nrow = 3, ncol = 2)
  test_data <- matrix(c(0.1, 0.9, 2.1, 0.1, 0.9, 2.1), nrow = 3, ncol = 2)
  result <- classification(class_means, test_data)
  expect_equal(result, c(1, 1))

  # 2/3 correct classification
  class_means <- matrix(c(0, 1, 2), nrow = 3, ncol = 1)
  test_data <- matrix(c(0.1, 1.9, 2.1), nrow = 3, ncol = 1)  # Second one wrong (closer to 2)
  result <- classification(class_means, test_data)
  expect_equal(result, 2/3)
})


test_that("shuffle2D shuffles rows respecting NaN", {
  set.seed(42)
  X <- matrix(1:20, nrow = 5, ncol = 4)
  X_original <- X

  # Shuffle should change the order
  X_shuffled <- shuffle2D(X)

  # Check that all original values are still present (just reordered by rows)
  expect_true(all(sort(X_shuffled[, 1]) == sort(X_original[, 1])))

  # With NA values, only non-NA rows should shuffle
  set.seed(42)
  X <- matrix(c(1, 2, NA, 4, 5,
                 6, 7, NA, 9, 10,
                 11, 12, NA, 14, 15,
                 16, 17, NA, 19, 20), nrow = 5, ncol = 4)
  X_shuffled <- shuffle2D(X)

  # Row 3 (with NA in first column) should be unchanged
  expect_true(is.na(X_shuffled[3, 1]))
  expect_true(is.na(X_shuffled[3, 2]))
})
