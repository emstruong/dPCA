#' Shuffle 2D Matrix In-Place
#'
#' Shuffles rows of a 2D matrix while respecting NaN values.
#' Only shuffles rows that do not have NaN in the first column.
#'
#' @param X A numeric matrix to shuffle (modified in place)
#' @return The shuffled matrix (invisibly, since modification is in-place)
#' @export
shuffle2D <- function(X) {
  # Find rows where the first column is not NA
  idx <- which(!is.na(X[, 1]))
  K <- ncol(X)
  TT <- length(idx)

  if (TT <= 1) {
    return(invisible(X))
  }

  # Fisher-Yates shuffle
  for (i in seq(TT, 2, by = -1)) {
    j <- sample.int(i, 1)
    n <- idx[i]
    m <- idx[j]
    # Swap rows
    temp <- X[n, ]
    X[n, ] <- X[m, ]
    X[m, ] <- temp
  }

  invisible(X)
}


#' Classification by Nearest Class Mean
#'
#' Classify test points according to the nearest class mean.
#' Used for significance testing in dPCA.
#'
#' @param class_means A matrix of class means (Q classes x T time points)
#' @param test A matrix of test data (Q test points x T time points)
#' @return A numeric vector of classification performance for each time point
#' @export
classification <- function(class_means, test) {
  Q <- nrow(class_means)
  TT <- ncol(class_means)
  performance <- numeric(TT)

  # Classify every data point in test according to nearest class mean
  for (t in seq_len(TT)) {
    for (p in seq_len(Q)) {
      argmin <- 1
      distance <- abs(class_means[argmin, t] - test[p, t])

      # Find closest class mean
      if (Q > 1) {
        for (q in 2:Q) {
          new_dist <- abs(class_means[q, t] - test[p, t])
          if (new_dist < distance) {
            distance <- new_dist
            argmin <- q
          }
        }
      }

      # Add 1 if class is correct
      if (argmin == p) {
        performance[t] <- performance[t] + 1
      }
    }
    performance[t] <- performance[t] / Q
  }

  performance
}


#' Denoise Significance Mask
#'
#' Remove isolated significant points that don't have enough consecutive neighbors.
#'
#' @param mask A binary vector (0/1) indicating significance
#' @param n_consecutive Minimum number of consecutive significant points required
#' @return A denoised binary vector
#' @export
denoise_mask <- function(mask, n_consecutive) {
  subseq <- 0
  N <- length(mask)
  result <- mask

  for (n in seq_len(N)) {
    if (result[n] == 1) {
      subseq <- subseq + 1
    } else {
      if (subseq < n_consecutive && subseq > 0) {
        # Set the previous subsequence to 0
        for (k in seq(n - subseq, n - 1)) {
          result[k] <- 0
        }
      }
      subseq <- 0
    }
  }

  # Handle trailing sequence
  if (subseq > 0 && subseq < n_consecutive) {
    for (k in seq(N - subseq + 1, N)) {
      result[k] <- 0
    }
  }

  result
}
