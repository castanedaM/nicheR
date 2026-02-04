#' Find extreme covariance limits for a specific pair
#'
#' @description
#' Calculates the maximum and minimum covariance values between two 
#' variables that maintain a positive definite (PD) matrix.
#'
#' @usage
#' covariance_pairs(varcov_matrix, i = 1, j = 2, tol = 1e-6)
#'
#' @param varcov_matrix A square numerical matrix.
#' @param i Integer. Row index of the covariance.
#' @param j Integer. Column index of the covariance.
#' @param tol Small value to subtract from limits to ensure strict PD.
#' 
#' @return A numeric vector with min and max covariance limits.
#'
#' @noRd

covariance_pairs <- function(varcov_matrix, i = 1, j = 2, tol = 1e-6) {
  n <- nrow(varcov_matrix)
  
  # 2D Logic: Deterministic Solution
  if (n == 2) {
    limit <- sqrt(varcov_matrix[1, 1] * varcov_matrix[2, 2])
    return(c(-limit + tol, limit - tol))
  }
  
  # N-Dimensional Case: Numerical Solution
  det_func <- function(x, m, r, c) {
    m[r, c] <- x
    m[c, r] <- x
    return(det(m))
  }
  
  ## Initial search limit based on 2D variances
  search_l <- sqrt(varcov_matrix[i, i] * varcov_matrix[j, j])
  
  ## Find roots with interval extension "yes" to ensure sign change
  ### Max covariance (positive root)
  up_r <- uniroot(det_func, lower = 0, upper = search_l, 
                  m = varcov_matrix, r = i, c = j, 
                  extendInt = "yes")$root
  
  ### Min covariance (negative root)
  lo_r <- uniroot(det_func, lower = -search_l, upper = 0, 
                  m = varcov_matrix, r = i, c = j, 
                  extendInt = "yes")$root
  
  return(c(lo_r + tol, up_r - tol))
}



#' Calculate safe covariance ranges for positive definiteness
#'
#' @description
#' Identifies the maximum and minimum covariance values that maintain a 
#' positive definite (PD) variance-covariance matrix. For niches with more 
#' than two dimensions, it independently shrinks the positive and negative 
#' theoretical limits until the "all-maximum" and "all-minimum" covariance 
#' scenarios are globally valid.
#'
#' @usage
#' covariance_limits(varcov_matrix, tol = 1e-6)
#'
#' @param varcov_matrix A square numerical matrix with variances on diagonal.
#' @param tol Small value to subtract from limits to ensure strict PD.
#' 
#' @return A data frame with columns 'min' and 'max', where each row 
#' represents a pair of dimensions named using the column names of 
#' the input matrix.
#' 
#' @details 
#' The function begins by calculating the deterministic 2D limits for each 
#' pair of variables ($|cov_{ij}| < \sqrt{\sigma_i^2 \sigma_j^2}$). 
#'
#' In higher dimensions (> 2), satisfying all pairwise limits is necessary 
#' but not sufficient to ensure the entire matrix is PD. Specifically, 
#' contradictory negative correlations can lead to non-transitivity errors. 
#' To provide a "safe zone" for users, the function independently shrinks 
#' the maximum and minimum bounds by 1% increments. 
#' 
#' It tests the "all-maximum" and "all-minimum" matrices using Cholesky 
#' decomposition. If a test fails, the respective bounds are reduced until 
#' a globally valid state is found. This independent approach ensures that 
#' geometric constraints in the negative correlation space do not 
#' unnecessarily penalize the limits in the positive correlation space.
#' 
#' @export
#'
#' @examples
#' # Example 1: 2D Matrix (Standard Deterministic Limits)
#' v_mat_2d <- matrix(c(10, 0, 0, 5), nrow = 2)
#' colnames(v_mat_2d) <- c("temp", "precip")
#' covariance_limits(v_mat_2d)
#'
#' # Example 2: 3D Matrix (Independent Shrinkage Exploration)
#' v_mat_3d <- matrix(c(10, 0, 0, 0, 5, 0, 0, 0, 7), nrow = 3)
#' colnames(v_mat_3d) <- c("temp", "precip", "hum")
#' covariance_limits(v_mat_3d)

covariance_limits <- function(varcov_matrix, tol = 1e-6) {
  n <- nrow(varcov_matrix)
  d_names <- colnames(varcov_matrix)
  if (is.null(d_names)) d_names <- paste0("dim", 1:n)
  
  # 1. Calculate theoretical 2D limits
  mins <- numeric()
  maxs <- numeric()
  row_labels <- character()
  
  for (j in 1:(n - 1)) {
    for (i in (j + 1):n) {
      row_labels <- c(row_labels, paste0(d_names[j], "-", d_names[i]))
      limit <- sqrt(varcov_matrix[i, i] * varcov_matrix[j, j]) - tol
      mins <- c(mins, -limit)
      maxs <- c(maxs, limit)
    }
  }
  
  # 2. Independent Global Safety Check for n > 2
  if (n > 2) {
    shrink_max <- 1.0
    shrink_min <- 1.0
    pd_max_ok <- FALSE
    pd_min_ok <- FALSE
    
    ## Shrink Maximums (Positive Space)
    while (!pd_max_ok && shrink_max > 0) {
      test_max <- varcov_matrix
      counter <- 1
      for (j in 1:(n - 1)) {
        for (i in (j + 1):n) {
          test_max[i, j] <- test_max[j, i] <- maxs[counter] * shrink_max
          counter <- counter + 1
        }
      }
      if (!inherits(try(chol(test_max), silent = TRUE), "try-error")) {
        pd_max_ok <- TRUE
        maxs <- maxs * shrink_max
      } else {
        shrink_max <- shrink_max - 0.01
      }
    }
    
    ## Shrink Minimums (Negative Space)
    while (!pd_min_ok && shrink_min > 0) {
      test_min <- varcov_matrix
      counter <- 1
      for (j in 1:(n - 1)) {
        for (i in (j + 1):n) {
          test_min[i, j] <- test_min[j, i] <- mins[counter] * shrink_min
          counter <- counter + 1
        }
      }
      if (!inherits(try(chol(test_min), silent = TRUE), "try-error")) {
        pd_min_ok <- TRUE
        mins <- mins * shrink_min
      } else {
        shrink_min <- shrink_min - 0.01
      }
    }
  }
  
  res_df <- data.frame(min = mins, max = maxs)
  rownames(res_df) <- row_labels
  return(res_df)
}



#' Update covariance values and calculate remaining safe limits
#'
#' @description
#' Updates a variance-covariance matrix with specific values and identifies 
#' the safe limits for all remaining zero-covariance combinations.
#'
#' @usage
#' update_covariance(varcov_matrix, covariance = 0, tol = 1e-6)
#'
#' @param varcov_matrix A square numerical matrix with variances on diagonal.
#' @param covariance A single value to apply to all off-diagonals, or a 
#' named vector where names follow the "dim1-dim2" format.
#' @param tol Small value to subtract from limits to ensure strict PD.
#' 
#' @return A list containing \code{up_mat} (updated matrix) and \code{limits} 
#' (data frame of safe limits for remaining zero-cov pairs).
#' 
#' @export
#'
#' @examples
#' # Setup
#' v_mat_3d <- matrix(c(10, 0, 0, 0, 5, 0, 0, 0, 7), 3, 3)
#' colnames(v_mat_3d) <- c("temp", "precip", "hum")
#'
#' # 1. 2D Update: Change only the single available covariance
#' update_covariance(v_mat_3d[1:2, 1:2], covariance = c("temp-precip" = 3))
#'
#' # 2. 3D Full Update: Set all combinations to 2.0
#' update_covariance(v_mat_3d, covariance = 2.0)
#'
#' # 3. 3D Partial: Define two, find limits for the remaining precip-hum
#' update_covariance(v_mat_3d, covariance = c("temp-precip" = -4, "temp-hum" = 5))

update_covariance <- function(varcov_matrix, covariance = 0, tol = 1e-6) {
  n <- nrow(varcov_matrix)
  d_names <- colnames(varcov_matrix)
  if (is.null(d_names)) d_names <- paste0("dim", 1:n)
  
  new_mat <- varcov_matrix
  
  # 1. Update the matrix with user-provided covariances
  if (length(covariance) == 1 && is.null(names(covariance))) {
    ## Case: Single value for all off-diagonals
    new_mat[lower.tri(new_mat)] <- covariance
    new_mat[upper.tri(new_mat)] <- t(new_mat)[upper.tri(new_mat)]
  } else if (!is.null(names(covariance))) {
    ## Case: Named vector for specific pairs
    for (name in names(covariance)) {
      parts <- strsplit(name, "-")[[1]]
      i <- which(d_names == parts[2])
      j <- which(d_names == parts[1])
      
      if (length(i) == 1 && length(j) == 1) {
        new_mat[i, j] <- new_mat[j, i] <- covariance[name]
      } else {
        warning(paste("Combination", name, "not found in matrix names."))
      }
    }
  }

  # 2. Immediate PD Check
  if (inherits(try(chol(new_mat), silent = TRUE), "try-error")) {
    stop("The provided covariance values result in a non-positive definite matrix.")
  }

  # 3. Calculate limits for remaining combinations (where cov == 0)
  mins <- numeric()
  maxs <- numeric()
  row_labels <- character()
  
  for (j in 1:(n - 1)) {
    for (i in (j + 1):n) {
      if (new_mat[i, j] == 0) {
        row_labels <- c(row_labels, paste0(d_names[j], "-", d_names[i]))
        
        ## Use numerical logic to find current safe range
        lims <- covariance_pairs(new_mat, i, j, tol)
        mins <- c(mins, lims[1])
        maxs <- c(maxs, lims[2])
      }
    }
  }
  
  res_df <- if (length(row_labels) > 0) {
    data.frame(min = mins, max = maxs, row.names = row_labels)
  } else {
    NULL # All combinations are already non-zero
  }
  
  return(list(updated_matrix = new_mat, remaining_limits = res_df))
}