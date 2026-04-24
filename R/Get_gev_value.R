#' Inverse GEV
#'
#'@description
#' Returns value for a given probability
#'
#'@details
#' By integrating a GEV tail fit, quantiles beyond the range of a historical data can be estimated.
#' For a given probability, the GEV distribution is defined with a sample data and return inverse values.
#'
#' @param tau Single or Array (or matrix) of probability for extreme values beyond a historical range
#' @param sample_data Array(or matrix) of sample data (observed or modelled)
#' @param upper logical value to indicate the location of extreme values (tails on the right side or both sides). A variable with zero bound (e.g., precipitation) is upper=T.
#'
#' @return Inverse GEV values for given probability
#' @export

get_gev_value <- function(tau, sample_data,upper=T) {
  library(extRemes)
  if (upper) {
    threshold <- quantile(sample_data, 0.9, na.rm = TRUE)
    tail_data <- sample_data[sample_data > threshold]
  } else {
    threshold <- quantile(sample_data, 0.1, na.rm = TRUE)
    tail_data <- sample_data[sample_data < threshold]
  }

  fit <- fevd(tail_data, type = "GEV", method = "MLE")
  params <- distill(fit)

  # Use qevd (Quantile function for EVD) to calculate the value for tau
  return(qevd(tau,
              loc = params["location"],
              scale = params["scale"],
              shape = params["shape"],
              type = "GEV"))
}
