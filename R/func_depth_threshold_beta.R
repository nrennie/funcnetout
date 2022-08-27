#' Functional Depth Threshold Beta
#'
#' This function calculates the threshold for the functional depths in order to
#' classify outliers by bootstrapping and fitting a beta distribution. This approach is
#' faster than func_depth_threshold() but makes a Beta distribution assumption for the functional
#' depths which may not always be appropriate.
#'
#' @param data A matrix where each row is an observation of a time series.
#' @param times A vector of length equal to the number of columns as data specifying the time points recorded at.
#' @param perc The percentile of the distribution to use. Default 0.01.
#' @param B The number of bootstraps to use. If 0, distribution is fitted to empirical depths. Default 10.
#' @return A numeric giving the threshold.
#' @seealso [func_depth_threshold()]
#' @export
#'

func_depth_threshold_beta <- function(data,
                                      times,
                                      perc = 0.01,
                                      B = 10) {

  #calculate depths and weights
  diffs <- array(t(data), dim = c(ncol(data), nrow(data), 1))
  weights <- mrfDepth::mfd(diffs, time = times, type = "projdepth")$MFDdepthZ
  if (B == 0) {
    m <- mean(weights)
    s <- stats::var(weights)
    alpha <- (m ^ 2) * (((1 - m) / s) - (1 / m))
    beta <- alpha * ((1 / m) - 1)
    threshold <- stats::qbeta(perc, shape1 = alpha, shape2 = beta)
    threshold
  } else {

    #create an empty list to store bootstrap samples
    w <- numeric(length = B)
    diff1 <- as.matrix(diffs[, , 1])

    for (i in 1:B) {
      #create bootstrap samples
      b <- diff1[, sample(seq_len(ncol(diff1)),
                          size = ncol(diff1),
                          replace = TRUE,
                          prob = weights)]
      #calculate smoothing matrix
      s <- matrix(MASS::mvrnorm(nrow(data), rep(0, ncol(data)), Sigma = 0.05 * stats::cov(data)),
                  nrow = ncol(data),
                  ncol = nrow(data),
                  byrow = TRUE)
      #smooth samples
      y <- array(b + s, dim = c(ncol(data), nrow(data), 1))
      #calculate percentile of each set of weights
      k <- mrfDepth::mfd(y, time = times, type = "projdepth")$MFDdepthZ
      m <- mean(k)
      s <- stats::var(k)
      alpha <- (m ^ 2) * (((1 - m) / s) - (1 / m))
      beta <- alpha * ((1 / m) - 1)
      #return percentile of distribution
      w[[i]] <- stats::qbeta(perc, shape1 = alpha, shape2 = beta)
    }

    #return median
    return(stats::median(w))
  }
}
