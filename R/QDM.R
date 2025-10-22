#' Quantile Delta Mapping (QDM)
#'
#'@description
#' QDM bias correction for preserving relative changes in quantiles

#'@details
#' QDM is equivalent to the equidistant and equiratio forms of quantile mapping (Cannon et al., 2015)
#'
#' @references
#' Cannon, A.J., Sobie, S.R., and Murdock, T.Q. 2015. Bias correction of
#' simulated precipitation by quantile mapping: How well do methods preserve
#' relative changes in quantiles and extremes? Journal of Climate, 28: 6938-6959. doi:10.1175/JCLI-D-14-00754.1
#'
#' @param o.h Vector of observed values during a historical period
#' @param m.h Vector of modelled values during a historical period
#' @param m.p Vector of modelled values during a projection period
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
#'
#' @return mhat.h Vector of bias corrected m.h values for the historical period
#' @return mhat.p Vector of bias corrected m.p values for the projection period
#' @export
#'
#' @examples
#'
QDM <-function(o.h, m.h, m.p, ratio=FALSE, trace=0.05, trace.calc=0.5*trace,
           jitter.factor=0, n.tau=NULL, ratio.max=2, ratio.max.trace=10*trace,
           ECBC=FALSE, ties='first', subsample=NULL, pp.type=7){

    Not.Missing<-which(!is.na(o.h))
    o.h<-o.h[Not.Missing]

    if(jitter.factor==0 &&
       (length(unique(o.h))==1 ||
        length(unique(m.h))==1 ||
        length(unique(m.p))==1)){
      jitter.factor <- sqrt(.Machine$double.eps)
    }
    if(jitter.factor > 0){
      o.h <- jitter(o.h, jitter.factor)
      m.h <- jitter(m.h, jitter.factor)
      m.p <- jitter(m.p, jitter.factor)
    }
    # For ratio data, treat exact zeros as left censored values less than
    # trace.calc
    if(ratio){
      epsilon <- .Machine$double.eps
      o.h[o.h < trace.calc] <- runif(sum(o.h < trace.calc), min=epsilon,
                                     max=trace.calc)
      m.h[m.h < trace.calc] <- runif(sum(m.h < trace.calc), min=epsilon,
                                     max=trace.calc)
      m.p[m.p < trace.calc] <- runif(sum(m.p < trace.calc), min=epsilon,
                                     max=trace.calc)
    }
    # Calculate empirical quantiles
    n <- length(m.p)
    if(is.null(n.tau)) n.tau <- n
    tau <- seq(0, 1, length=n.tau)
    if(!is.null(subsample)){
      quant.o.h <- rowMeans(apply(replicate(subsample,
                                            sample(o.h, size=length(tau))),
                                  2, quantile, probs=tau, type=pp.type))
      quant.m.h <- rowMeans(apply(replicate(subsample,
                                            sample(m.h, size=length(tau))),
                                  2, quantile, probs=tau, type=pp.type))
      quant.m.p <- rowMeans(apply(replicate(subsample,
                                            sample(m.p, size=length(tau))),
                                  2, quantile, probs=tau, type=pp.type))
    } else{
      quant.o.h <- quantile(o.h, tau, type=pp.type)
      quant.m.h <- quantile(m.h, tau, type=pp.type)
      quant.m.p <- quantile(m.p, tau, type=pp.type)
    }
    # Apply quantile delta mapping bias correction
    tau.m.p <- approx(quant.m.p, tau, m.p, rule=2)$y
    if(ratio){
      approx.t.qmc.tmp <- approx(tau, quant.m.h, tau.m.p, rule=2)$y
      delta.m <- m.p/approx.t.qmc.tmp
      delta.m[(delta.m > ratio.max) &
                (approx.t.qmc.tmp < ratio.max.trace)] <- ratio.max
      mhat.p <- approx(tau, quant.o.h, tau.m.p, rule=2)$y*delta.m
    } else{
      delta.m <- m.p - approx(tau, quant.m.h, tau.m.p, rule=2)$y
      mhat.p <- approx(tau, quant.o.h, tau.m.p, rule=2)$y + delta.m
    }
    mhat.h <- approx(quant.m.h, quant.o.h, m.h, rule=2)$y
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
