#' Functional Depth
#'
#' This function calculates the functional depths.
#' @param data A matrix where each row is an observation of a time series. Row names should be a unique ID.
#' @param times A vector of length equal to the number of columns as data specifying the time points recorded at.
#' @return A numeric vector of depths.
#' @export

func_depth <- function(data, times) {
  #check data format
  if (length(unique(rownames(data))) != nrow(data)) {
    stop("Row names must be a unique datetime")
  }

  d <- array(t(data), dim = c(ncol(data), nrow(data), 1))
  fit <- mrfDepth::mfd(d, time = times, type = "projdepth")
  depths <- fit$MFDdepthX
  rownames(depths) <- rownames(data)
  return(depths)
}
