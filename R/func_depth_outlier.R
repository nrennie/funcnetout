#' Functional Depth Outliers
#'
#' This function calculates the functional depth and returns the outliers.
#' @param data A matrix where each row is an observation of a time series. Row names should be a unique ID.
#' @param times A vector of length equal to the number of columns as data specifying the time points recorded at.
#' @param threshold The choice of threshold calculation. One of c("beta", "bootstrap"). Default "beta".
#' @param perc The percentile of the distribution to use. Default 0.01.
#' @param B The number of bootstraps to use. If 0, distribution is fitted to empirical depths. Default 10.
#' @param maxiter Maximum number of iterations for outlier detection. Default 10.
#' @seealso [func_depth_threshold()]
#' @seealso [func_depth_threshold_beta()]
#' @return A vector of outliers based on rownames.
#' @export

func_depth_outlier <- function(data,
                               times,
                               threshold = "beta",
                               perc = 0.01,
                               B = 10,
                               maxiter = 10) {

  #check data format
  if (length(unique(rownames(data))) != nrow(data)) {
    stop("Row names must be a unique datetime")
  }

  #calculate threshold
  if (threshold == "beta") {
    C <- func_depth_threshold_beta(data, times = times, perc = perc, B = B)
  } else if (threshold == "bootstrap") {
    C <- func_depth_threshold_beta(data, times = times, perc = perc, B = B)
  }

  #format and calculate depths
  d <- array(t(data), dim = c(ncol(data), nrow(data), 1))
  rownames(d) <- seq_len(nrow(data))
  fit <- mrfDepth::mfd(d, time = times, type = "projdepth")
  depths <- fit$MFDdepthX

  #find outliers
  outliers_int <- which(depths <= C)
  outliers <- rownames(data)[outliers_int]
  i <- 0

  #iteratively remove outliers
  while (length(which(depths <= C)) > 0 && i < maxiter && length(outliers) < nrow(data)) {
    i <- i + 1
    d1 <- data[-outliers_int, ]
    d <- array(t(d1), dim = c(ncol(d1), nrow(d1), 1))
    fit <- mrfDepth::mfd(d, time = times, type = "projdepth")
    depths <- fit$MFDdepthX
    outliers <- c(outliers, rownames(d1[which(depths <= C), ]))
  }

  #return outliers
  return(outliers)
}
