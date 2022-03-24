#' Functional extrapolation
#'
#' This function extrapolates partially observed data using historic and arima processes.
#'
#' @param data A matrix where each row is an observation of a time series. Row names should be a unique datetime.
#' @param alpha A weighting function to weight between historic and arima forecasts
#' @return A boolean vector.
#' @export

func_extrapolation <- function(data,
                               alpha = 0) {

  # alpha should be same length as ncol(data)
  if (length(alpha) != ncol(data)) {
    stop("Length of weighting vector must be equal to number of time points")
    }

  # prep data
  d <- data

  #remove all rows with no observations
  d1 <- d[apply(d, 1, function(x) !all(is.na(x))), ]

  #check which cases are completely observed i.e. are historic data
  d_obs <- d[which(stats::complete.cases(d1)), ]
  if (nrow(d_obs) == nrow(data)) {
    return(extrapolation = data)
  }

  #run extrapolation on remaining values
  d_ext <- d[which(!stats::complete.cases(d1)), ]

  #obtain general forecast based on historic data (completely observed)
  g_mat <- historic_forecast(d_obs)

  #for each partially observed run extrap
  extrap <- matrix(NA, nrow = nrow(d_ext), ncol = ncol(d_ext))
  colnames(extrap) <- colnames(d_obs)

  for (i in seq_len(nrow(d_ext))) {
    x <- d_ext[i, ]
    n_obs <- length(x[!is.na(x)])

    #historic forecast
    hist_g <- which(weekdays(as.Date(rownames(g_mat))) == weekdays(as.Date(rownames(d_ext)[i])))
    g <- colMeans(g_mat[hist_g, ])[(n_obs + 1):(length(x))]

    #arima forecast
    f <- as.vector(forecast::forecast(forecast::auto.arima(as.numeric(x[1:(n_obs)])),
                            h = ncol(data) - (n_obs))$mean)

    #weighted forecast
    fc <- as.numeric((1 - (alpha[n_obs])) * g + (alpha[n_obs]) * f)
    extrap[i, ] <- as.numeric(c(x[1:n_obs], fc))
  }

  #join observed and extrapolated together
  extrapolation <- rbind(d_obs, extrap)
  return(extrapolation = extrapolation)
}
