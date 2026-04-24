#' Bias-Correction for operational forecasts with extreme tails
#'
#'@description
#' Quantile Mapping (QM) based bias-Correction considering sub-seasonal or seasonal forecasts with extremes out of historical range

#'@details
#' While QM is a non-parametric approach, when forecasted extreme events are out of historical range, those extremes are fitted to an extreme probabilistic distribution (GEV) to estimate non-exceedance probability.
#' Thus, this is a hybrid approach with non-parametric and parametric approaches.This approach is to bias-correct for subseasonal or seasonal forecasts.
#'
#'
#' @references
#'
#' @param o.h Vector of observed values during a historical period
#' @param m.h Vector of modelled values during a hindcast period
#' @param m.f Vector of forecast values during a forecast period (1-month to 6-month ahead forecasts)
#' @param ratio logical value indicating if samples are of a ratio quantity. TRUE: Preserving relative ratio (e.g., precipitation), FALSE: Preserving relative absolute value (e.g., temperature)
#' @param trace Replace values less than trace with exact zeros (defaults to 0.05)
#' @param trace.calc Treat values below trace.calc as censored (defaults to 0.5*trace)
#' @param jitter.factor jitter to accommodate ties (defaults to 0.01). If jitter.factor > 0, apply a small amount of jitter to accommodate ties due to limited measurement precision.
#' @param n.tau NULL --> number of empirical quantiles (NULL=sample length)
#' @param ratio.max numeric value indicating the maximum proportional change allowed for ratio variables (defaults to 2)
#' @param ratio.max.trace Values below which ratio.max is applied (default 10*trace)
#' @param ECBC logical value indicating whether mhat.p outputs should be ordered according to o.h ranks using the empirical copula-bias correction (ECBC) algorithm (defaults to FALSE)
#' @param ties Method used to handle ties when calculating ordinal ranks (defaults to 'first')
#' @param subsample use subsample draws of size n.tau to calculate empirical quantiles; if NULL, calculate normally
#' @param pp.type type of plotting position used in quantile (defaults to 7)
#' @param upper.tail logical value to identify extreme values in the upper tail (e.g., precipitation and wind (minimum bound is zero) = T, temperature=F)
#'
#' @return mhat.h Vector of bias corrected m.h values for the historical period
#' @return mhat.p Vector of bias corrected m.f values for the forecast period
#' @export
#'
#' @examples
#'
BC_Forecasts_Tails <-function(o.h, m.h,m.f, ratio=FALSE, trace=0.05, trace.calc=0.5*trace,
           jitter.factor=0, n.tau=NULL, ratio.max=2, ratio.max.trace=10*trace,
           ECBC=FALSE, ties='first', subsample=NULL, pp.type=7,upper.tail=T){

    Not.Missing<-which(!is.na(o.h))
    o.h<-o.h[Not.Missing]

    if(jitter.factor==0 &&
       (length(unique(o.h))==1 ||
        length(unique(m.h))==1 ||
        length(unique(m.f))==1)){
      jitter.factor <- sqrt(.Machine$double.eps)
    }
    if(jitter.factor > 0){
      o.h <- jitter(o.h, jitter.factor)
      m.h <- jitter(m.h, jitter.factor)
      m.f <- jitter(m.f, jitter.factor)
    }
    # For ratio data, treat exact zeros as left censored values less than
    # trace.calc
    if(ratio){
      epsilon <- .Machine$double.eps
      o.h[o.h < trace.calc] <- runif(sum(o.h < trace.calc), min=epsilon,
                                     max=trace.calc)
      m.h[m.h < trace.calc] <- runif(sum(m.h < trace.calc), min=epsilon,
                                     max=trace.calc)
      m.f[m.f < trace.calc] <- runif(sum(m.f < trace.calc), min=epsilon,
                                     max=trace.calc)
    }
    # Calculate empirical quantiles
    n <- length(m.h)
    if(is.null(n.tau)) n.tau <- n
    tau <- seq(0, 1, length=n.tau)
    if(!is.null(subsample)){
      quant.o.h <- rowMeans(apply(replicate(subsample,
                                            sample(o.h, size=length(tau))),
                                  2, quantile, probs=tau, type=pp.type))
      quant.m.h <- rowMeans(apply(replicate(subsample,
                                            sample(m.h, size=length(tau))),
                                  2, quantile, probs=tau, type=pp.type))
      # quant.m.f <- rowMeans(apply(replicate(subsample,
      #                                       sample(m.h, size=length(tau))),
      #                             2, quantile, probs=tau, type=pp.type))
    } else{
      quant.o.h <- quantile(o.h, tau, type=pp.type)
      quant.m.h <- quantile(m.h, tau, type=pp.type)
      #quant.m.f <- quantile(m.f, tau, type=pp.type)
    }
    # Apply quantile mapping bias correction

    tau.m.f <- approx(quant.m.h, tau, m.f, rule=2)$y

    if(ratio){
      mhat.p <- approx(tau, quant.o.h, tau.m.f, rule=2)$y
    } else{
      mhat.p <- approx(tau, quant.o.h, tau.m.f, rule=2)$y
    }
    mhat.h <- approx(quant.m.h, quant.o.h, m.h, rule=2)$y

    #---------------------------------------------------
    # GEV tail extrapolation
    if (upper.tail) {
      Index.Extremes<-which (m.f > max(m.h))
      if (length(Index.Extremes) > 0) {
        tau.f.extreme<-get_gev_tau(m.f[Index.Extremes], m.h)
        mhat.p.extreme<-get_gev_value(tau.f.extreme, o.h)
        mhat.p[Index.Extremes]<-mhat.p.extreme
      }
    } else { # for temperature
      Index.Extremes.upper<-which (m.f > max(m.h))
      Index.Extremes.lower<-which (m.f < min(m.h))
      if (length(Index.Extremes.upper) > 0) {
        tau.f.extreme<-get_gev_tau(m.f[Index.Extremes.upper], m.h,upper=T)
        mhat.p.extreme<-get_gev_value(tau.f.extreme, o.h,upper=T)
        mhat.p[Index.Extremes.upper]<-mhat.p.extreme
      }
      if (length(Index.Extremes.lower) > 0) {
        tau.f.extreme<-get_gev_tau(m.f[Index.Extremes.lower], m.h,upper=F)
        mhat.p.extreme<-get_gev_value(tau.f.extreme, o.h,upper=F)
        mhat.p[Index.Extremes.lower]<-mhat.p.extreme
      }
    }

    #-------------------------------------------
    # For ratio data, set values less than trace to zero
    if(ratio){
      mhat.h[mhat.h < trace] <- 0
      mhat.p[mhat.p < trace] <- 0
    }
    if(ECBC){
      # empirical copula coupling/Schaake shuffle
      if(length(mhat.p)==length(o.h)){
        mhat.p <- sort(mhat.p)[rank(o.h, ties.method=ties)]
      } else{
        stop('Schaake shuffle failed due to incompatible lengths')
      }
    }
    list(mhat.h=mhat.h, mhat.p=mhat.p)
  }
