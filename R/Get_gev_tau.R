#' Returns the probability for extreme values
#'
#'@description
#' Fits GEV to the upper or lower tail and returns the probability for extreme values out of a historical range
#'
#'@details
#' By integrating a GEV tail fit, quantiles beyond the range of a historical data can be estimated.
#' The GEV distribution is defined by three parameters: location, scale, and shape.
#' The shape parameter is particularly important as it determines the behavior of the tail.
#'
#' @param y_val Single or Array (or matrix) of extreme values beyond a historical range
#' @param sample_data Array(or matrix) of sample data (observed or modelled)
#' @param upper logical value to indicate the location of extreme values (tails on the right side or both sides). A variable with zero bound (e.g., precipitation) is upper=T.
#'
#' @return tau_extrapolated: Estimated probability for the given extreme values
#' @export

get_gev_tau <- function(y_val, sample_data, upper=T) {
  library(extRemes)
  # Fit to the top 10% of the data to characterize the tail behavior
  if (upper) {
    threshold <- quantile(sample_data, 0.9, na.rm = TRUE)
    tail_data <- sample_data[sample_data > threshold]
  } else {
    threshold <- quantile(sample_data, 0.1, na.rm = TRUE)
    tail_data <- sample_data[sample_data < threshold]
  }

  # Fit GEV using Maximum Likelihood Estimation (MLE) [cite: 558, 568]
  # type = "GEV" specifies the Generalized Extreme Value distribution [cite: 814]
  fit <- fevd(tail_data, type = "GEV", method = "MLE")
  params <- distill(fit)

  # Use pevd (Probability function for EVD) to calculate F(x)
  # This corresponds to G(x; mu, sigma, xi) in the paper
  tau_extrapolated <- pevd(y_val,
                           loc = params["location"],
                           scale = params["scale"],
                           shape = params["shape"],
                           type = "GEV")
  return(tau_extrapolated)
}

