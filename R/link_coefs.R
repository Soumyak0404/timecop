#' Compute Hermite coefficients and link function constants
#'
#' @param param_list List. A list of marginal parameters
#' @param k Numeric. The value at which Hermite coefficient infinite sums terminate. Default is 100
#' @param family_list List. List of marginal distributions
#' @param corr Logical. Correlations or covariances
#' @return A numeric vector of k link function coefficients.
#' @keywords internal

link_coefs <- function(param_list, k, family_list, corr){

  # get SDs of marginal distributions
  sd_marg <- lapply(seq_along(family_list), function(i) {

    param <- param_list[[i]]

    if (family_list[[i]] == "Bernoulli") {

      sd_marg <- sqrt(param * (1 - param))

    } else if (family_list[[i]] == "Poisson") {

      sd_marg <- sqrt(param)

    } else if (family_list[[i]] == "Ordinal") {

      # For X in {0,1,2}, param = c(p0, p1, p2)
      # E[X]   = p1 + 2 p2
      # E[X^2] = p1 + 4 p2
      # Var(X) = E[X^2] - E[X]^2

      if (length(param) != 3) {
        stop("For Ordinal family, param must be c(p0, p1, p2).", call. = FALSE)
      }

      if (any(param < 0) || abs(sum(param) - 1) > 1e-8) {
        stop("For Ordinal family, probabilities must be nonnegative and sum to 1.", call. = FALSE)
      }

      EX <- param[2] + 2 * param[3]
      EX2 <- param[2] + 4 * param[3]

      sd_marg <- sqrt(EX2 - EX^2)

    } else if (family_list[[i]] == "Gaussian") {

      sd_marg <- 1
    }

    return(sd_marg)
  })

  # get Hermite coefficients
  g_i <- hermite_coefs(param_list[[1]], k, family_list[[1]])
  g_j <- hermite_coefs(param_list[[2]], k, family_list[[2]])

  # Polynomial of link function
  if (corr == FALSE) {

    l_ij <- g_i * g_j * factorial(1:k)

  } else {

    l_ij <- g_i * g_j * factorial(1:k) / (sd_marg[[1]] * sd_marg[[2]])
  }

  return(l_ij)
}

