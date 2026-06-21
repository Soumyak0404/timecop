#' Estimate marginal parameters from data
#'
#' @param data Matrix. A d by n multivariate time series matrix
#' @param family List. A list of length d with the names of each distribution
#' @param d Numeric. The number of variables
#' @return A list of length d of estimated marginal parameters.
#' @keywords internal

estimate_marginal_params <- function(data, family, d) {
  param_hat <- vector("list", d)
  for (i in seq_len(d)) {
    if (family[[i]] %in% c("Bernoulli", "Poisson")) {
      param_hat[[i]] <- mean(data[i, ], na.rm = TRUE)
    } else if (family[[i]] == "Ordinal") {

      x_i <- data[i, ]

      if (any(!x_i %in% c(0, 1, 2))) {
        stop("Ordinal data must only contain values 0, 1, and 2.", call. = FALSE)
      }

      counts <- tabulate(match(x_i, c(0, 1, 2)), nbins = 3)

      param_hat[[i]] <- counts / sum(counts)

    } else if (family[[i]] == "Gaussian") {
      param_hat[[i]] <- c(mean(data[i, ], na.rm = TRUE),
                          var(data[i, ],  na.rm = TRUE))
    }
  }
  param_hat
}
