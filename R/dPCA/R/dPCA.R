#' Demixed Principal Component Analysis (dPCA)
#'
#' dPCA is a linear dimensionality reduction technique that automatically discovers
#' and highlights the essential features of complex population activities. The
#' population activity is decomposed into a few demixed components that capture most
#' of the variance in the data and that highlight the dynamic tuning of the population
#' to various task parameters, such as stimuli, decisions, rewards, etc.
#'
#' @param labels Character string of labels for feature axes (e.g., "ts" for time and stimulus)
#'               or an integer specifying the number of labels to generate from alphabet
#' @param join Optional named list specifying parameter combinations to join.
#'             For example, list(ts = c("s", "ts")) joins stimulus and time-stimulus interactions.
#' @param n_components Number of components to keep. Can be a single integer (same for all
#'                     marginalizations) or a named list (e.g., list(t = 10, ts = 5)).
#' @param regularizer Regularization parameter. NULL or 0 for no regularization, numeric for
#'                    specific value (as multiple of data variance), or "auto" to optimize.
#' @param n_iter Number of iterations for randomized SVD solver.
#'
#' @return A dPCA object (list) containing model parameters and methods
#'
#' @examples
#' \dontrun{
#' # Create a dPCA model for time and stimulus
#' model <- dPCA(labels = "ts", n_components = 10)
#'
#' # Fit the model
#' model <- dpca_fit(model, X)
#'
#' # Transform data
#' Z <- dpca_transform(model, X)
#' }
#'
#' @export
dPCA <- function(labels = NULL, join = NULL, n_components = 10, regularizer = NULL, n_iter = 0) {
  # Create labels from alphabet if not provided as string
  if (is.numeric(labels) && length(labels) == 1) {
    alphabet <- letters
    labels <- paste0(alphabet[1:labels], collapse = "")
  } else if (!is.character(labels) || nchar(labels) == 0) {
    stop("Wrong type for labels. Please either set labels to the number of variables or provide the axis labels as a single string of characters (like 'ts' for time and stimulus)")
  }

  # Handle regularizer
  opt_regularizer_flag <- identical(regularizer, "auto")
  if (is.null(regularizer)) {
    regularizer <- 0
  } else if (identical(regularizer, "auto")) {
    regularizer <- 0
    message("You chose to determine the regularization parameter automatically. This can ",
            "take substantial time and grows linearly with the number of crossvalidation ",
            "folds. Use $n_trials to set the number of folds (default = 3). Similarly, ",
            "use $protect to set the list of axes not to shuffle during splitting.")
  }

  # Create model object
  model <- list(
    labels = labels,
    join = join,
    regularizer = regularizer,
    opt_regularizer_flag = opt_regularizer_flag,
    n_components = n_components,
    n_iter = n_iter,
    debug = 2,
    n_trials = 3,
    protect = NULL,
    P = NULL,  # Encoding matrices
    D = NULL,  # Decoding matrices
    explained_variance_ratio_ = NULL
  )

  # Compute marginalizations
  model$marginalizations <- get_parameter_combinations(model, join = TRUE)

  class(model) <- c("dPCA", "list")
  model
}


#' Get Parameter Combinations
#'
#' Returns all parameter combinations for the given labels.
#' For labels = "xyz", returns list(x = c(0), y = c(1), z = c(2),
#' xy = c(0,1), xz = c(0,2), yz = c(1,2), xyz = c(0,1,2))
#'
#' @param model A dPCA model object
#' @param join Logical, whether to condense according to model$join
#' @return A named list of parameter combinations (0-indexed like Python)
#' @keywords internal
get_parameter_combinations <- function(model, join = TRUE) {
  labels <- model$labels
  n_labels <- nchar(labels)
  label_chars <- strsplit(labels, "")[[1]]

  # Generate all subsets (0-indexed to match Python)
  indices <- 0:(n_labels - 1)

  # Get all non-empty subsets
  subsets <- list()
  for (r in 1:n_labels) {
    combs <- combn(indices, r, simplify = FALSE)
    subsets <- c(subsets, combs)
  }

  # Create named list (using set-like representation)
  pcombs <- list()
  for (subset in subsets) {
    key <- paste0(label_chars[subset + 1], collapse = "")
    pcombs[[key]] <- subset  # Store as numeric vector (0-indexed)
  }

  # Condense dict if join is specified
  if (!is.null(model$join) && join) {
    for (key in names(model$join)) {
      combs <- model$join[[key]]
      tmp <- lapply(combs, function(comb) pcombs[[comb]])

      # Remove the combined keys
      for (comb in combs) {
        pcombs[[comb]] <- NULL
      }

      # Add the joined key with list of subsets
      pcombs[[key]] <- tmp
    }
  }

  pcombs
}


#' Marginalize Data Matrix
#'
#' Compute marginalized versions of the data for each parameter combination.
#'
#' @param model A dPCA model object
#' @param X Array of shape (n_samples, n_features_1, n_features_2, ...)
#' @param save_memory Logical, use memory-saving mode (slower)
#' @return Named list of marginalized data matrices
#' @keywords internal
dpca_marginalize <- function(model, X, save_memory = FALSE) {
  labels <- model$labels
  n_labels <- nchar(labels)
  label_chars <- strsplit(labels, "")[[1]]

  # Helper function: mean over multiple axes with optional expansion

  mmean <- function(Z, axes, expand = FALSE) {
    result <- Z
    # Sort axes in descending order for safe sequential application
    for (ax in sort(axes, decreasing = TRUE)) {
      result <- apply(result, setdiff(seq_along(dim(result)), ax), mean)
      if (expand) {
        # Re-expand the dimension
        new_dims <- dim(Z)
        new_dims[ax] <- 1
        if (is.null(dim(result))) {
          result <- array(result, dim = 1)
        }
        result <- array(result, dim = new_dims[-ax])
        # Insert dimension back
        result <- aperm(
          array(result, dim = c(dim(result)[1:(ax-1)], 1, if(ax <= length(dim(result))) dim(result)[ax:length(dim(result))] else NULL)),
          c(seq_len(ax-1), ax, if(ax <= length(dim(result))) (ax):length(dim(result)) + 1 else NULL)
        )
      }
    }
    result
  }

  # Copy and center data
  Xres <- X
  n_features <- dim(X)[1]
  flat_X <- array(X, dim = c(n_features, prod(dim(X)[-1])))
  row_means <- rowMeans(flat_X)
  Xres <- X - array(row_means, dim = dim(X))

  # Get parameter combinations without join
  pcombs <- get_parameter_combinations(model, join = FALSE)

  # Full set of indices (0-indexed)
  S <- pcombs[[names(pcombs)[length(pcombs)]]]

  # Initialize marginalizations
  Xmargs <- list()

  if (save_memory) {
    for (key in names(pcombs)) {
      phi <- pcombs[[key]]
      S_without_phi <- setdiff(S, phi)

      if (length(S_without_phi) > 0) {
        # Mean over axes not in phi (convert to 1-indexed for R, +2 because first axis is features)
        axes_to_mean <- S_without_phi + 2  # +1 for R indexing, +1 for feature axis
        Xmargs[[key]] <- apply(Xres, setdiff(seq_along(dim(Xres)), axes_to_mean), mean)
        # Expand back to original dimensions
        target_dims <- dim(Xres)
        target_dims[axes_to_mean] <- 1
        Xmargs[[key]] <- array(Xmargs[[key]], dim = target_dims)
      } else {
        Xmargs[[key]] <- Xres
      }
      Xres <- Xres - Xmargs[[key]]
    }
  } else {
    # Efficient precomputation of means
    pre_mean <- list()

    for (key in names(pcombs)) {
      phi <- pcombs[[key]]
      if (length(key) == 1) {
        # Single character key - compute mean directly
        axis_to_mean <- phi[1] + 2  # +1 for R indexing, +1 for feature axis
        other_axes <- setdiff(seq_along(dim(Xres)), axis_to_mean)
        mean_result <- apply(Xres, other_axes, mean)
        target_dims <- dim(Xres)
        target_dims[axis_to_mean] <- 1
        pre_mean[[key]] <- array(mean_result, dim = target_dims)
      } else {
        # Multiple characters - compute iteratively from previous
        prev_key <- substr(key, 1, nchar(key) - 1)
        last_phi <- phi[length(phi)]
        axis_to_mean <- last_phi + 2

        other_axes <- setdiff(seq_along(dim(pre_mean[[prev_key]])), axis_to_mean)
        if (length(other_axes) > 0) {
          mean_result <- apply(pre_mean[[prev_key]], other_axes, mean)
        } else {
          mean_result <- mean(pre_mean[[prev_key]])
        }
        target_dims <- dim(pre_mean[[prev_key]])
        target_dims[axis_to_mean] <- 1
        pre_mean[[key]] <- array(mean_result, dim = target_dims)
      }
    }

    # Compute marginalizations
    for (key in names(pcombs)) {
      phi <- pcombs[[key]]
      # Get characters not in key
      key_without_phi <- paste0(setdiff(label_chars, strsplit(key, "")[[1]]), collapse = "")

      # Get the appropriate pre-mean
      if (nchar(key_without_phi) > 0) {
        X_base <- pre_mean[[key_without_phi]]
      } else {
        X_base <- Xres
      }

      if (nchar(key) > 1) {
        # Subtract all proper subsets
        key_chars <- strsplit(key, "")[[1]]
        # Get all proper subsets of key
        proper_subsets <- character(0)
        for (r in 1:(length(key_chars) - 1)) {
          subset_combs <- combn(key_chars, r, simplify = FALSE)
          proper_subsets <- c(proper_subsets, sapply(subset_combs, paste0, collapse = ""))
        }

        result <- X_base
        for (subset in proper_subsets) {
          result <- result - Xmargs[[subset]]
        }
        Xmargs[[key]] <- result
      } else {
        Xmargs[[key]] <- X_base
      }
    }
  }

  # Condense according to join
  if (!is.null(model$join)) {
    for (key in names(model$join)) {
      combs <- model$join[[key]]

      # Determine target shape
      Xshape <- rep(1, n_labels + 1)
      for (comb in combs) {
        sh <- dim(Xmargs[[comb]])
        non_one <- which(sh > 1)
        Xshape[non_one] <- sh[non_one]
      }

      tmp <- array(0, dim = Xshape)
      for (comb in combs) {
        tmp <- tmp + Xmargs[[comb]]
        Xmargs[[comb]] <- NULL
      }
      Xmargs[[key]] <- tmp
    }
  }

  # Convert to dense 2D format (n_features x n_conditions)
  for (key in names(Xmargs)) {
    # Expand to full size if needed
    if (!identical(dim(Xmargs[[key]]), dim(X))) {
      Xmargs[[key]] <- array(Xmargs[[key]], dim = dim(X))
    }
    Xmargs[[key]] <- matrix(Xmargs[[key]], nrow = n_features)
  }

  Xmargs
}


#' Add Regularization to Data
#'
#' Prepares data matrices for the randomized dPCA solver by adding regularization.
#'
#' @param Y Data array
#' @param mYs Marginalized data list
#' @param lam Regularization parameter
#' @param SVD Optional SVD decomposition
#' @return List with regularized data, marginalizations, and pseudo-inverse
#' @keywords internal
add_regularization <- function(Y, mYs, lam, SVD = NULL) {
  n_features <- dim(Y)[1]
  flat_Y <- matrix(Y, nrow = n_features)

  # Add regularization to Y: [Y, lam * I]
  regY <- cbind(flat_Y, lam * diag(n_features))

  # Add zeros to marginalizations
  regmYs <- list()
  for (key in names(mYs)) {
    regmYs[[key]] <- cbind(mYs[[key]], matrix(0, nrow = n_features, ncol = n_features))
  }

  # Compute pseudo-inverse
  if (!is.null(SVD)) {
    U <- SVD$u
    s <- SVD$d
    V <- SVD$v

    M <- t(U) / (s^2 + lam^2)
    pregY <- rbind(t(V) * s, lam * U) %*% M
  } else {
    YYt <- flat_Y %*% t(flat_Y)
    pregY <- t(regY) %*% solve(YYt + lam^2 * diag(n_features))
  }

  list(regY = regY, regmYs = regmYs, pregY = pregY)
}


#' Randomized dPCA Solver
#'
#' Solves the dPCA minimization problem using randomized SVD.
#'
#' @param model A dPCA model object
#' @param X Data array
#' @param mXs Marginalized data list
#' @param pinvX Pseudo-inverse of X (optional)
#' @return List with P (encoding) and D (decoding) matrices
#' @keywords internal
randomized_dpca <- function(model, X, mXs, pinvX = NULL) {
  n_features <- dim(X)[1]
  rX <- matrix(X, nrow = n_features)

  if (is.null(pinvX)) {
    pinvX <- MASS::ginv(rX)
  }

  P <- list()
  D <- list()

  for (key in names(mXs)) {
    mX <- matrix(mXs[[key]], nrow = n_features)
    C <- mX %*% pinvX

    # Determine number of components
    if (is.list(model$n_components)) {
      n_comp <- model$n_components[[key]]
    } else {
      n_comp <- model$n_components
    }

    # Randomized SVD using irlba
    CX <- C %*% rX
    svd_result <- irlba::irlba(CX, nv = n_comp, nu = n_comp)

    U <- svd_result$u
    P[[key]] <- U
    D[[key]] <- t(t(U) %*% C)
  }

  list(P = P, D = D)
}


#' Fit dPCA Model
#'
#' Fit the dPCA model to training data.
#'
#' @param model A dPCA model object
#' @param X Array of shape (n_samples, n_features_1, n_features_2, ...)
#' @param trialX Optional trial-by-trial data for regularization optimization
#' @param mXs Optional pre-computed marginalizations
#' @param center Logical, whether to center the data
#' @param optimize Logical, whether to optimize regularization
#'
#' @return The fitted dPCA model object
#' @export
dpca_fit <- function(model, X, trialX = NULL, mXs = NULL, center = TRUE, optimize = TRUE) {
  n_features <- dim(X)[1]
  labels <- model$labels

  # Center data
  if (center) {
    flat_X <- matrix(X, nrow = n_features)
    row_means <- rowMeans(flat_X)
    X <- X - array(row_means, dim = dim(X))
  }

  # Marginalize data
  if (is.null(mXs)) {
    mXs <- dpca_marginalize(model, X)
  }

  # Optimize regularization if requested
  if (model$opt_regularizer_flag && optimize) {
    if (model$debug > 0) {
      message("Start optimizing regularization.")
    }
    if (is.null(trialX)) {
      stop("To optimize the regularization parameter, the trial-by-trial data trialX needs to be provided.")
    }
    model <- optimize_regularization(model, X, trialX, center = FALSE)
  }

  # Add regularization
  if (model$regularizer > 0) {
    reg_result <- add_regularization(X, mXs, model$regularizer * sum(X^2))
    regX <- reg_result$regY
    regmXs <- reg_result$regmYs
    pregX <- reg_result$pregY
  } else {
    regX <- X
    regmXs <- mXs
    pregX <- MASS::ginv(matrix(X, nrow = n_features))
  }

  # Compute closed-form solution
  result <- randomized_dpca(model, regX, regmXs, pinvX = pregX)
  model$P <- result$P
  model$D <- result$D

  model
}


#' Transform Data Using dPCA Model
#'
#' Apply dimensionality reduction to data using fitted dPCA model.
#'
#' @param model A fitted dPCA model object
#' @param X Array of shape (n_samples, n_features_1, n_features_2, ...)
#' @param marginalization Optional specific marginalization to return
#'
#' @return Named list of transformed data for each marginalization, or single array
#' @export
dpca_transform <- function(model, X, marginalization = NULL) {
  # Zero-mean the data
  n_features <- dim(X)[1]
  flat_X <- matrix(X, nrow = n_features)
  row_means <- rowMeans(flat_X)
  X <- X - array(row_means, dim = dim(X))

  total_variance <- sum((X - mean(X))^2)

  # Helper to compute variance explained
  marginal_variances <- function(marginal) {
    D <- model$D[[marginal]]
    Xr <- matrix(X, nrow = n_features)
    sapply(1:ncol(D), function(k) {
      sum((D[, k] %*% Xr)^2) / total_variance
    })
  }

  if (!is.null(marginalization)) {
    D <- model$D[[marginalization]]
    Xr <- matrix(X, nrow = n_features)
    X_transformed <- t(D) %*% Xr
    X_transformed <- array(X_transformed, dim = c(ncol(D), dim(X)[-1]))
    model$explained_variance_ratio_ <- list()
    model$explained_variance_ratio_[[marginalization]] <- marginal_variances(marginalization)
    return(X_transformed)
  } else {
    X_transformed <- list()
    model$explained_variance_ratio_ <- list()

    for (key in names(model$marginalizations)) {
      D <- model$D[[key]]
      Xr <- matrix(X, nrow = n_features)
      transformed <- t(D) %*% Xr
      X_transformed[[key]] <- array(transformed, dim = c(ncol(D), dim(X)[-1]))
      model$explained_variance_ratio_[[key]] <- marginal_variances(key)
    }
    return(X_transformed)
  }
}


#' Fit and Transform dPCA Model
#'
#' Fit the dPCA model and apply dimensionality reduction in one step.
#'
#' @param model A dPCA model object
#' @param X Array of shape (n_samples, n_features_1, n_features_2, ...)
#' @param trialX Optional trial-by-trial data
#'
#' @return List with transformed data and updated model
#' @export
dpca_fit_transform <- function(model, X, trialX = NULL) {
  model <- dpca_fit(model, X, trialX = trialX)
  Z <- dpca_transform(model, X)
  list(model = model, Z = Z)
}


#' Inverse Transform dPCA Data
#'
#' Transform data back to original space.
#'
#' @param model A fitted dPCA model object
#' @param X Transformed data array
#' @param marginalization The marginalization used for transformation
#'
#' @return Array in original data space
#' @export
dpca_inverse_transform <- function(model, X, marginalization) {
  # Zero-mean
  n_components <- dim(X)[1]
  flat_X <- matrix(X, nrow = n_components)
  row_means <- rowMeans(flat_X)
  X <- X - array(row_means, dim = dim(X))

  P <- model$P[[marginalization]]
  Xr <- matrix(X, nrow = n_components)
  X_transformed <- P %*% Xr
  array(X_transformed, dim = c(nrow(P), dim(X)[-1]))
}


#' Reconstruct Data Through dPCA
#'
#' Transform data to reduced space and back. Equivalent to
#' inverse_transform(transform(X)).
#'
#' @param model A fitted dPCA model object
#' @param X Original data array
#' @param marginalization The marginalization to use
#'
#' @return Reconstructed data array
#' @export
dpca_reconstruct <- function(model, X, marginalization) {
  Z <- dpca_transform(model, X, marginalization = marginalization)
  dpca_inverse_transform(model, Z, marginalization)
}


#' Get Number of Samples per Condition
#'
#' Compute number of valid (non-NA) samples for each parameter combination.
#'
#' @param model A dPCA model object
#' @param trialX Trial-by-trial data array
#' @param protect Axes to protect during splitting
#' @return Array of sample counts
#' @keywords internal
get_n_samples <- function(model, trialX, protect = NULL) {
  n_unprotect <- length(dim(trialX)) - length(protect) - 1
  n_protect <- if (!is.null(protect)) length(protect) else 0

  if (n_protect > 0) {
    # Count non-NA samples
    # Index to select first element along protected axes
    idx <- c(rep(TRUE, 1 + n_unprotect), rep(1, n_protect))
    subset_X <- do.call(`[`, c(list(trialX), as.list(c(TRUE, rep(TRUE, n_unprotect), rep(1, n_protect)))))
    N_samples <- dim(trialX)[1] - apply(subset_X, seq(2, n_unprotect + 1), function(x) sum(is.na(x)))
  } else {
    N_samples <- dim(trialX)[1] - apply(trialX, seq(2, length(dim(trialX))), function(x) sum(is.na(x)))
  }

  N_samples
}


#' Train-Test Split for dPCA
#'
#' Split data into training and validation sets for cross-validation.
#'
#' @param model A dPCA model object
#' @param X Mean data array
#' @param trialX Trial-by-trial data array
#' @param N_samples Optional pre-computed sample counts
#'
#' @return List with trainX and validX arrays
#' @export
dpca_train_test_split <- function(model, X, trialX, N_samples = NULL) {
  protect <- model$protect
  n_unprotect <- length(dim(X)) - if (!is.null(protect)) length(protect) else 0
  n_protect <- if (!is.null(protect)) length(protect) else 0

  # Compute number of samples if not provided
  if (is.null(N_samples)) {
    N_samples <- get_n_samples(model, trialX, protect = protect)
  }

  # Random indices for validation set
  idx <- array(floor(runif(length(N_samples)) * as.vector(N_samples)) + 1, dim = dim(N_samples))

  # Select validation data
  validX <- array(NA, dim = dim(trialX)[-1])

  # Iterate over conditions
  it <- arrayInd(seq_along(N_samples), dim(N_samples))
  for (i in seq_len(nrow(it))) {
    cond_idx <- it[i, ]
    trial_idx <- idx[matrix(cond_idx, nrow = 1)]

    # Build index for trialX
    full_idx <- c(trial_idx, cond_idx)
    if (n_protect > 0) {
      full_idx <- c(full_idx, rep(TRUE, n_protect))
    }

    # Build index for validX
    valid_idx <- cond_idx
    if (n_protect > 0) {
      valid_idx <- c(valid_idx, rep(TRUE, n_protect))
    }

    # Extract and assign
    validX[matrix(valid_idx, nrow = 1)] <- trialX[matrix(full_idx, nrow = 1)]
  }

  # Compute training data
  N_samples_expanded <- array(N_samples, dim = dim(X))
  trainX <- X * (N_samples_expanded / (N_samples_expanded - 1)) - validX / (N_samples_expanded - 1)

  # Re-center both datasets
  n_features <- dim(X)[1]
  train_means <- rowMeans(matrix(trainX, nrow = n_features))
  valid_means <- rowMeans(matrix(validX, nrow = n_features))
  trainX <- trainX - array(train_means, dim = dim(X))
  validX <- validX - array(valid_means, dim = dim(X))

  list(trainX = trainX, validX = validX)
}


#' Shuffle Labels In-Place
#'
#' Shuffle labels between conditions in trial-by-trial data.
#'
#' @param model A dPCA model object
#' @param trialX Trial-by-trial data array (will be modified)
#'
#' @return The shuffled trial data
#' @export
dpca_shuffle_labels <- function(model, trialX) {
  protect <- model$protect
  n_protect <- if (!is.null(protect)) length(protect) else 0

  # Reshape to 2D for shuffling
  original_shape <- dim(trialX)

  if (n_protect > 0) {
    # Reshape: (n_trials * n_unprotected, n_protected)
    n_unprotected_dims <- length(original_shape) - n_protect
    new_shape <- c(prod(original_shape[1:n_unprotected_dims]), prod(original_shape[(n_unprotected_dims + 1):length(original_shape)]))
    trialX <- matrix(trialX, nrow = new_shape[1], ncol = new_shape[2])
  } else {
    trialX <- matrix(trialX, nrow = prod(original_shape), ncol = 1)
  }

  # Shuffle
  trialX <- shuffle2D(trialX)

  # Reshape back
  array(trialX, dim = original_shape)
}


#' Scoring Function for Cross-Validation
#'
#' Compute reconstruction error for cross-validation.
#'
#' @param model A fitted dPCA model object
#' @param X Validation data
#' @param mXs Marginalized validation data
#' @param mean Logical, return mean over marginalizations
#' @return Reconstruction error (single value or named list)
#' @keywords internal
dpca_score <- function(model, X, mXs, mean = TRUE) {
  n_features <- dim(X)[1]
  flat_X <- matrix(X, nrow = n_features)

  error <- list()
  PDY <- list()
  trPD <- list()

  for (key in names(mXs)) {
    PDY[[key]] <- model$P[[key]] %*% (t(model$D[[key]]) %*% flat_X)
    trPD[[key]] <- rowSums(model$P[[key]] * model$D[[key]])
    error[[key]] <- sum((mXs[[key]] - PDY[[key]] + outer(trPD[[key]], rep(1, ncol(flat_X))) * flat_X)^2)
  }

  if (mean) {
    sum(unlist(error))
  } else {
    error
  }
}


#' Cross-Validation Score
#'
#' Calculate cross-validation scores for different regularization parameters.
#'
#' @param model A dPCA model object
#' @param lams Vector of regularization parameters to test
#' @param X Mean data array
#' @param trialX Trial-by-trial data array
#' @param mean Logical, average over marginalizations
#'
#' @return Matrix or list of cross-validation scores
#' @export
dpca_crossval_score <- function(model, lams, X, trialX, mean = TRUE) {
  n_trials <- model$n_trials

  if (mean) {
    scores <- matrix(0, nrow = n_trials, ncol = length(lams))
  } else {
    scores <- lapply(names(model$marginalizations), function(key) {
      matrix(0, nrow = n_trials, ncol = length(lams))
    })
    names(scores) <- names(model$marginalizations)
  }

  # Get sample counts
  N_samples <- get_n_samples(model, trialX, protect = model$protect)

  for (trial in seq_len(n_trials)) {
    message("Starting trial ", trial, "/", n_trials)

    # Train-test split
    split_result <- dpca_train_test_split(model, X, trialX, N_samples = N_samples)
    trainX <- split_result$trainX
    validX <- split_result$validX

    # Marginalize
    trainmXs <- dpca_marginalize(model, trainX)
    validmXs <- dpca_marginalize(model, validX)

    for (k in seq_along(lams)) {
      # Fit with this lambda
      model$regularizer <- lams[k]
      model <- dpca_fit(model, trainX, mXs = trainmXs, optimize = FALSE)

      # Score
      if (mean) {
        scores[trial, k] <- dpca_score(model, validX, validmXs, mean = TRUE)
      } else {
        tmp <- dpca_score(model, validX, validmXs, mean = FALSE)
        for (key in names(model$marginalizations)) {
          scores[[key]][trial, k] <- tmp[[key]]
        }
      }
    }
  }

  scores
}


#' Optimize Regularization Parameter
#'
#' Find optimal regularization parameter via cross-validation.
#'
#' @param model A dPCA model object
#' @param X Mean data array
#' @param trialX Trial-by-trial data array
#' @param center Logical, center data first
#' @param lams Vector of lambda values to test, or "auto"
#' @return Updated model with optimal regularizer
#' @keywords internal
optimize_regularization <- function(model, X, trialX, center = TRUE, lams = "auto") {
  # Center data
  if (center) {
    n_features <- dim(X)[1]
    flat_X <- matrix(X, nrow = n_features)
    row_means <- rowMeans(flat_X)
    X <- X - array(row_means, dim = dim(X))
  }

  # Generate lambda values if auto
  if (identical(lams, "auto")) {
    N <- 45
    lams <- 1.4^(0:(N-1)) * 1e-7
  }

  # Cross-validation
  scores <- dpca_crossval_score(model, lams, X, trialX, mean = FALSE)

  # Total score across marginalizations
  total_score <- Reduce(`+`, lapply(scores, function(s) rowMeans(s)))
  total_score <- colMeans(matrix(total_score, ncol = length(lams)))

  # Warning if at boundary
  opt_idx <- which.min(total_score)
  if (opt_idx == 1 || opt_idx == length(total_score)) {
    if (model$debug > 0) {
      warning("Optimal regularization parameter lies at the boundary of the search interval. ",
              "Please provide different search list (key: lams).")
    }
  }

  model$regularizer <- lams[opt_idx]

  if (model$debug > 1) {
    message("Optimized regularization, optimal lambda = ", model$regularizer)
    message("Regularization will be fixed; to compute the optimal parameter again, ",
            "set opt_regularizer_flag to TRUE.")
  }

  model$opt_regularizer_flag <- FALSE
  model
}


#' Significance Analysis
#'
#' Cross-validated significance analysis of dPCA model using classification.
#'
#' @param model A fitted dPCA model object
#' @param X Mean data array
#' @param trialX Trial-by-trial data array
#' @param n_shuffles Number of label shuffles (default 100)
#' @param n_splits Number of train-test splits per shuffle (default 100)
#' @param n_consecutive Minimum consecutive significant points required (default 1)
#' @param axis NULL or TRUE for time-resolved analysis
#' @param full Logical, return all scores
#'
#' @return List of significance masks (and optionally scores)
#' @export
dpca_significance_analysis <- function(model, X, trialX, n_shuffles = 100, n_splits = 100,
                                       n_consecutive = 1, axis = NULL, full = FALSE) {
  stopifnot(is.null(axis) || isTRUE(axis))

  # Helper to compute mean score
  compute_mean_score <- function(X, trialX, n_splits) {
    K <- if (is.null(axis)) 1 else dim(X)[length(dim(X))]

    if (is.numeric(model$n_components)) {
      scores <- lapply(keys, function(key) {
        array(NA, dim = c(model$n_components, n_splits, K))
      })
    } else {
      scores <- lapply(keys, function(key) {
        array(NA, dim = c(model$n_components[[key]], n_splits, K))
      })
    }
    names(scores) <- keys

    for (shuffle in seq_len(n_splits)) {
      cat(".")

      # Train-validation split
      split_result <- dpca_train_test_split(model, X, trialX)
      trainX <- split_result$trainX
      validX <- split_result$validX

      # Fit and transform
      fit_result <- dpca_fit_transform(model, trainX)
      model <- fit_result$model
      trainZ <- fit_result$Z
      validZ <- dpca_transform(model, validX)

      for (key in keys) {
        ncomps <- if (is.numeric(model$n_components)) model$n_components else model$n_components[[key]]

        # Get axes to average over (not in key)
        axset <- model$marginalizations[[key]]
        if (is.list(axset)) {
          axset <- unique(unlist(axset))
        }
        axes_to_avg <- setdiff(0:(length(dim(X)) - 2), axset)

        # Average over non-key axes
        for (ax in sort(axes_to_avg + 2, decreasing = TRUE)) {
          trainZ[[key]] <- apply(trainZ[[key]], setdiff(seq_along(dim(trainZ[[key]])), ax), mean)
          validZ[[key]] <- apply(validZ[[key]], setdiff(seq_along(dim(validZ[[key]])), ax), mean)
        }

        # Reshape for classification
        if ((length(dim(X)) - 2) %in% axset && !is.null(axis)) {
          trainZ[[key]] <- array(trainZ[[key]], dim = c(ncomps, prod(dim(trainZ[[key]])[-1]) / K, K))
          validZ[[key]] <- array(validZ[[key]], dim = c(ncomps, prod(dim(validZ[[key]])[-1]) / K, K))
        } else {
          trainZ[[key]] <- array(trainZ[[key]], dim = c(ncomps, prod(dim(trainZ[[key]])[-1]), 1))
          validZ[[key]] <- array(validZ[[key]], dim = c(ncomps, prod(dim(validZ[[key]])[-1]), 1))
        }

        # Classification
        for (comp in seq_len(ncomps)) {
          scores[[key]][comp, shuffle, ] <- classification(trainZ[[key]][comp, , , drop = FALSE],
                                                            validZ[[key]][comp, , , drop = FALSE])
        }
      }
    }

    # Average over splits
    for (key in keys) {
      scores[[key]] <- apply(scores[[key]], c(1, 3), mean, na.rm = TRUE)
    }

    scores
  }

  # Optimize regularization if needed
  if (model$opt_regularizer_flag) {
    message("Regularization not optimized yet; start optimization now.")
    model <- optimize_regularization(model, X, trialX)
  }

  keys <- setdiff(names(model$marginalizations), model$labels[nchar(model$labels)])

  # Copy trial data for shuffling
  trialX <- trialX + 0  # Force copy

  # Compute score of original data
  message("Compute score of data: ", appendLF = FALSE)
  true_score <- compute_mean_score(X, trialX, n_splits)
  message(" Finished.")

  # Shuffled scores
  scores <- lapply(keys, function(key) list())
  names(scores) <- keys

  for (it in seq_len(n_shuffles)) {
    message("\rCompute score of shuffled data: ", it, "/", n_shuffles, appendLF = FALSE)

    # Shuffle labels
    trialX <- dpca_shuffle_labels(model, trialX)

    # Recompute mean
    X <- apply(trialX, seq(2, length(dim(trialX))), mean, na.rm = TRUE)

    score <- compute_mean_score(X, trialX, n_splits)

    for (key in keys) {
      scores[[key]] <- c(scores[[key]], list(score[[key]]))
    }
  }
  message("")

  # Compute masks
  masks <- list()
  for (key in keys) {
    max_score <- apply(simplify2array(scores[[key]]), 1:2, max)
    masks[[key]] <- true_score[[key]] > max_score
  }

  # Denoise masks
  if (n_consecutive > 1) {
    for (key in keys) {
      for (k in seq_len(nrow(masks[[key]]))) {
        masks[[key]][k, ] <- denoise_mask(as.integer(masks[[key]][k, ]), n_consecutive)
      }
    }
  }

  if (full) {
    list(masks = masks, true_score = true_score, scores = scores)
  } else {
    masks
  }
}


#' Print dPCA Object
#'
#' @param x A dPCA object
#' @param ... Additional arguments (ignored)
#' @export
print.dPCA <- function(x, ...) {
  cat("dPCA Model\n")
  cat("  Labels:", x$labels, "\n")
  cat("  Marginalizations:", paste(names(x$marginalizations), collapse = ", "), "\n")
  cat("  n_components:", if (is.list(x$n_components)) paste(names(x$n_components), unlist(x$n_components), sep = "=", collapse = ", ") else x$n_components, "\n")
  cat("  Regularizer:", x$regularizer, "\n")
  cat("  Fitted:", !is.null(x$P), "\n")
  invisible(x)
}
