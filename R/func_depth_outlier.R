#' Functional Depth Outliers
#'
#' This function negates the in function and returns true if the item before notin is not in the item after.
#'
#' @param data A matrix where each row is an observation of a time series. Row names should be a unique ID.
#' @param times A vector of length equal to the number of columns as data specifying the time points recorded at.
#' @param threshold The choice of threshold calculation. One of c("beta", "bootstrap"). Default "beta".
#' @param perc The percentile of the distribution to use. Default 0.01.
#' @param B The number of bootstraps to use. If 0, distribution is fitted to empirical depths. Default 10.
#' @param maxiter Maximum number of iterations for outlier detection. Default 10.
#' @return A numeric giving the threshold.
#' @seealso [func_depth_threshold()]
#' @seealso [func_depth_threshold_beta()]
#' @return A boolean vector.
#' @export
#'
#'

func_depth_outlier <- function(data, times, threshold="beta", perc=0.01, B=10, maxiter=10){

  #check data format
  if (length(unique(rownames(data))) != nrow(data)) {
    stop('Row names must be a unique ID')
  }

  #calculate threshold
  if (threshold == "beta"){
    C <- func_depth_threshold_beta(data, times=times, perc=perc, B=B)
  } else if (threshold == "bootstrap") {
    C <- func_depth_threshold_beta(data, times=times, perc=perc, B=B)
  }

  #format and calculate depths
  d <- array(t(data), dim=c(ncol(data),nrow(data),1))
  fit <- mfd(d)
  depths <- fit$MFDdepthX

  #find outliers
  outliers <- which(depths <= C)
  i = 0

  #iteratively remove outliers
  while (length(which(depths <= C))>0 && i < maxiter){
    i = i + 1
    d1 <- data[-outliers,]
    d <- array(t(d1), dim=c(ncol(d1),nrow(d1),1))
    fit <- mfd(d)
    depths <- fit$MFDdepthX
    outliers <- c(outliers, as.integer(rownames(d1[which(depths <= C),])))
  }

  #return outliers
  outliers
}
