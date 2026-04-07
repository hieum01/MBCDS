#' Multivariate Bias Correction with Distribution-free Shuffle approach (MBCDS)
#'@description
#' MBCDS bias correction for preserving relative changes in quantiles and interdependence between climate variables
#'
#'@details
#' QDM is applied to bias-corrrect GCM outputs and Distribution-free Shuffle is used for preserving interdependence between climate variables.
#' The struture of arguments are exactly same as QDM except for multi-variables (array or matrix) of o.h, m.h, and m.p while they are vectors (univariate).
#'
#' @references
#' MBCDS:
#' Eum, H.-I., Gupta, A., Dibike, Y., 2020. Effects of univariate and multivariate statistical downscaling methods on climatic and hydrologic indicators for Alberta, Canada. Journal of Hydrology 588, 125065. https://doi.org/10.1016/j.jhydrol.2020.125065
#'
#' QDM:
#' Cannon, A.J., Sobie, S.R., and Murdock, T.Q. 2015. Bias correction of simulated precipitation by quantile mapping: How well do methods preserve relative changes in quantiles and extremes? Journal of Climate, 28: 6938-6959. doi:10.1175/JCLI-D-14-00754.1
#'
#' @param o.h Array(or matrix) of observed values during a historical period
#' @param m.h Array(or matrix) of modelled values during a historical period
#' @param m.f Array(or matrix) of modelled values during a future (projection) period
#' @param ratio.seq Array of logical values indicating if samples are of a ratio quantity. TRUE: Preserving relative ratio (e.g., precipitation), FALSE: Preserving relative absolute value (e.g., temperature)
#' @param trace Replace values less than trace with exact zeros (defaults to 0.05)
#' @param trace.calc Treat values below trace.calc as censored (defaults to 0.5*trace)
#' @param jitter.factor jitter to accommodate ties (defaults to 0.01). If jitter.factor > 0, apply a small amount of jitter to accommodate ties due to limited measurement precision.
#' @param n.tau NULL --> number of empirical quantiles (NULL=sample length)
#' @param ratio.max numeric value indicating the maximum proportional change allowed for ratio variables (defaults to 2)
#' @param ratio.max.trace Values below which ratio.max is applied (default 10*trace)
#' @param ties Method used to handle ties when calculating ordinal ranks (defaults to 'first')
#' @param subsample use subsample draws of size n.tau to calculate empirical quantiles; if NULL, calculate normally
#' @param pp.type type of plotting position used in quantile (defaults to 7)
#' @param pr.i Index of precipitation variable
#'
#' @return mhat.h Array of bias corrected m.h values for the historical period
#' @return mhat.p Array of bias corrected m.p values for the projection period
#'
#' @export
#'
#' @examples
MBCDS <-
  function(o.h, m.h, m.f,ratio.seq=rep(FALSE, ncol(o.h)),trace=0.05,
           trace.calc=0.01, jitter.factor=0, n.tau=NULL, ratio.max=2,
           ratio.max.trace=10*trace, ties='first',subsample=NULL, pp.type=7,pr.i=1) {

    ratio.col<-which(ratio.seq)
    Threshold_Pr<-trace.calc

    if(length(trace.calc)==1)
      trace.calc <- rep(trace.calc, ncol(o.h))
    if(length(trace)==1)
      trace <- rep(trace, ncol(o.h))
    if(length(jitter.factor)==1)
      jitter.factor <- rep(jitter.factor, ncol(o.h))
    if(length(ratio.max))
      ratio.max <- rep(ratio.max, ncol(o.h))
    if(length(ratio.max.trace)==1)
      ratio.max.trace <- rep(ratio.max.trace, ncol(o.h))

    m.h.qmap <- m.h
    m.f.qmap <- m.f
    # Quantile delta mapping bias correction
    for(i in seq(ncol(o.h))){
      fit.qmap <- QDM(o.h=o.h[,i], m.h=m.h[,i], m.p=m.f[,i],
                      ratio=ratio.seq[i], trace.calc=trace.calc[i],
                      trace=trace[i], jitter.factor=jitter.factor[i],
                      n.tau=n.tau, ratio.max=ratio.max[i],
                      ratio.max.trace=ratio.max.trace[i],
                      subsample=subsample, pp.type=pp.type)
      m.h.qmap[,i] <- fit.qmap$mhat.h
      m.f.qmap[,i] <- fit.qmap$mhat.p
    }

    m.h <- m.h.qmap
    m.f <- m.f.qmap
    #=========================================================================================
    #Jitter to numbers of zero for ratio variables (e.g., precipitation)
    epsilon <- .Machine$double.eps
    o.c[o.c[,pr.i] < Threshold_Pr,pr.i] <- runif(sum(o.c[,pr.i] < Threshold_Pr), min=epsilon,
                                                 max=Threshold_Pr)
    m.c[m.c[,pr.i] < Threshold_Pr,pr.i] <- runif(sum(m.c[,pr.i] < Threshold_Pr), min=epsilon,
                                                 max=Threshold_Pr)
    m.f[m.f[,pr.i] < Threshold_Pr,pr.i] <- runif(sum(m.f[,pr.i] < Threshold_Pr), min=epsilon,
                                                 max=Threshold_Pr)
    #=========================================================================================
    Ncol<-ncol(o.h)
    Nrow<-nrow(o.h)
    #View(Data)
    Rank_o.h<-sapply(1:ncol(o.h),function(x) {rank(o.h[,x])})
    Rank_m.h<-sapply(1:ncol(m.h),function(x) {rank(m.h[,x])})
    Rank_m.f<-sapply(1:ncol(m.f),function(x) {rank(m.f[,x])})

    #----------------------------------
    # van der Waerden scores
    #----------------------------------
    p<-Rank_o.h/(Nrow+1)
    W_o.h<-qnorm(p,0,1)

    p<-Rank_m.h/(Nrow+1)
    W_m.h<-qnorm(p,0,1)

    p<-Rank_m.f/(Nrow+1)
    W_m.f<-qnorm(p,0,1)

    #==============================================
    # Confirm if dependence structure is preserved
    #==============================================
    #cor(W_o.h,method='spearman')
    #cor(W_m.h,method='spearman')
    #library(psych)
    #pairs.panels(W)
    #pairs.panels(W_m.f)
    #==============================================
    # Cholesky decomposition matrix
    #==============================================
    P_o.h<-as.matrix(chol(nearPD(cor(W_o.h,method='spearman'))$mat))  # need Matrix package for nearPD
    P_m.h<-as.matrix(chol(nearPD(cor(W_m.h,method='spearman'))$mat))
    P_m.f<-as.matrix(chol(nearPD(cor(W_m.f,method='spearman'))$mat))
    #==============================================
    # New score matrix
    #==============================================
    W_m.h.hat<-W_m.h%*% solve(P_m.h)%*%P_o.h
    W_m.f.hat<-W_m.f%*% solve(P_m.f)%*%P_o.h
    #========================================================================
    # Rearrangement of original data corresponding to the new score matrix
    #========================================================================
    m.h.hat<-m.h
    m.f.hat<-m.f
    for(i in seq(ncol(W_m.h.hat))){
      # Replace ordinal ranks with the original data
      m.h.hat[,i] <- sort(m.h[,i])[rank(W_m.h.hat[,i])]
      m.f.hat[,i] <- sort(m.f[,i])[rank(W_m.f.hat[,i])]
    }

    m.c.hat[m.c.hat[,pr.i] < Threshold_Pr,pr.i] <- 0.0
    m.f.hat[m.f.hat[,pr.i] < Threshold_Pr,pr.i] <- 0.0
    list(mhat.h=m.h.hat, mhat.p=m.f.hat)
  }

################################################################################
