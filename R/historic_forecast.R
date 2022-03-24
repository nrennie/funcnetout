#' Historic forecast
#'
#' This function computes a historic forecast using linear regression, with days of the week as variables.
#'
#' @param d_obs A matrix where each row is an observation of a time series. Row names should be a unique datetime.
#' @return A numeric matrix
#' @export

historic_forecast <- function(d_obs) {
  dat <- as.matrix(d_obs)
  dates <- as.Date(rownames(d_obs))

  # define factors based on day of week
  day <- factor(weekdays(dates)) #nolint

  #apply a regression to each time point
  fit_matrix <- apply(dat, 2, function(y) {
      mod <- stats::lm(y ~ day)
      mod_fit <- mod$fitted.values
      return(mod_fit)
    }
  )

  # prediction matrix
  f_mat <- cbind(d_obs, fit_matrix)
  pred_mat <- f_mat[!duplicated(round(f_mat, 10)), ]

  return(pred_mat)
}
