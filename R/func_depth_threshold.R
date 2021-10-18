#' Functional Depth Threshold
#'
#' This function calculates the threshold for the functional depths in order to
#' classify outliers using a bootstrap approach.
#'
#' @param data A matrix where each row is an observation of a time series.
#' @param times A vector of length equal to the number of columns as data specifying the time points recorded at.
#' @param perc The percentile of the distribution to use. Default 0.01.
#' @param B The number of bootstraps to use. If 0, distribution is fitted to empirical depths. Default 1000.
#' @return A numeric giving the threshold.
#' @seealso [func_depth_threshold_beta()]
#' @importFrom mrfDepth mfd
#' @export
#'
#'

func_depth_threshold <- function(data,times,perc=0.01,B=1000){
  if (B <= 0){
    stop('Number of bootstraps must be positive')
  }
  #calculate depths and weights
  diffs <- array(t(data), dim=c(ncol(data),nrow(data),1))
  weights <- mfd(diffs, time=times, type="projdepth")$MFDdepthZ
  #create an empty list to store bootstrap samples
  w <- numeric(length = B)
  diff1 <- as.matrix(diffs[,,1])
  for (i in 1:B){
    #create bootstrap samples
    b <- diff1[,sample(1:ncol(diff1),size=ncol(diff1),replace=TRUE,prob=weights)]
    #calculate smoothing matrix
    s <- matrix(mvrnorm(nrow(data),rep(0,ncol(data)),Sigma=0.05*cov(data)),nrow=ncol(data),ncol=n,byrow=TRUE)
    #smooth samples
    y <- array(b+s,dim=c(ncol(data),nrow(data),1))
    #calculate percentile of each set of weights
    k <- mfd(y, time=times, type="projdepth")$MFDdepthZ
    w[[i]] <- sort(k)[ceiling(perc*B)]
  }
  median(w)
}
