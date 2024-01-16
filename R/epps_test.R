#' The asymptotic Epps and Pulley Test for normality.
#'
#' Performs the asymptotic Epps test of normality for univariate time series.
#' Computes the p-value using the asymptotic Gamma Distribution.
#'
#' @usage epps.test(y, lambda = c(1,2))
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a stationary
#' time series.
#' @param lambda a numeric vector for evaluating the characteristic function. This values
#' could be selected by the user for a better test performance. By default, the values
#' are `c(1,2)`, another plausible option is to select random values.
#'
#' @return A list with class \code{"h.test"} containing the following components:
#'  \item{statistic:}{the Epps statistic.}
#'  \item{parameter:}{the test degrees freedoms.}
#'  \item{p.value:}{the p value.}
#'  \item{alternative:}{a character string describing the alternative hypothesis.}
#'  \item{method:}{a character string \dQuote{Epps test}.}
#'  \item{data.name:}{a character string giving the name of the data.}
#'
#' @details
#' The Epps test minimize the process' empirical characteristic function using a
#' quadratic loss in terms of the process two first moments.  \emph{Nieto-Reyes, A.,
#' Cuesta-Albertos, J. & Gamboa, F. (2014)} upgrade the test implementation by
#' allowing the option of evaluating the characteristic function with random values.
#' The \code{amoebam()} function of \emph{Press, W.H., Teukolsky, S.A., Vetterling,
#' W.T. and  Flannery, B.P. (2007)}, performs the optimal search.
#'
#' @export
#'
#' @author Asael Alonzo Matamoros and Alicia Nieto-Reyes.
#'
#' @seealso \code{\link{lobato.test}}
#'
#' @references
#' Epps, T.W. (1987). Testing that a stationary time series is Gaussian. \emph{The
#' Annals of Statistic}. 15(4), 1683-1698.
#'
#' Nieto-Reyes, A., Cuesta-Albertos, J. & Gamboa, F. (2014). A random-projection
#' based test of Gaussianity for stationary processes. \emph{Computational
#' Statistics & Data Analysis, Elsevier}, vol. 75(C), pages 124-141.
#'
#' Press, W.H., Teukolsky, S.A., Vetterling, W.T. and  Flannery, B.P. (2007).
#' Numerical Recipes. The Art of Scientific Computing. \emph{Cambridge
#' University Press}.
#'
#' @examples
#' # Generating an stationary arma process
#' y = arima.sim(100,model = list(ar = 0.3))
#' epps.test(y)
#'
#'# Epps tests with random lambda values
#' y = arima.sim(100,model = list(ar = c(0.3,0.2)))
#' epps.test(y, lambda = rnorm(2,mean = 1,sd = 0.1))
#'
epps.test = function(y, lambda = c(1,2)){

  if( !is.numeric(y) & !is(y,class2 = "ts") )
    stop("y object must be numeric or a time series")

  if( anyNA(y) )
    stop("The time series contains missing values")

  # checking stationarity
  cc = uroot.test(y)
  if(!cc$stationary)
    warning("y has a unit root, epps.test requires stationary process")

  # checking seasonality
  if(frequency(y) > 1){
    cc = seasonal.test(y)
    if(cc$seasonal)
      warning("y has a seasonal unit root, epps.test requires stationary process")
  }

  dname = deparse(substitute(y))
  alt = paste(dname,"does not follow a Gaussian Process")
  tstat = epps.statistic(y, lambda = lambda)
  names(tstat) <- "epps"
  df = 2
  names(df) <- "df"
  pval = pchisq(q = tstat,df = 2,lower.tail = FALSE)

  rval <- list(statistic =tstat, parameter = df, p.value = pval,
                 alternative = alt,
                 method = "Epps test", data.name = dname)
  class(rval) <- "htest"
  return(rval)
}
#' Estimates the Epps statistic.
#'
#' Estimates the Epps statistic minimizing the quadratic loss of the process'
#' characteristic function in terms of the first two moments.
#'
#' @usage epps.statistic(y, lambda = c(1,2))
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a stationary time series.
#' @param lambda a numeric vector for evaluating the characteristic function. This values
#' could be selected by the user for a better test performance. By default, the values
#' are `c(1,2)`, another plausible option is to select random values.
#'
#' @details
#' The Epps test minimize the process' empirical characteristic function using a
#' quadratic loss in terms of the process two first moments. \emph{Nieto-Reyes, A.,
#' Cuesta-Albertos, J. & Gamboa, F. (2014)} upgrade the test implementation by
#' allowing the option of evaluating the characteristic function with random values.
#'
#' This function is the equivalent of \code{Sub} in \emph{Nieto-Reyes, A.,
#' Cuesta-Albertos, J. & Gamboa, F. (2014)}. This function uses a quadratic
#' optimization solver implemented by \emph{Press, W.H., Teukolsky, S.A.,
#' Vetterling, W.T. and  Flannery, B.P. (2007)}.
#'
#' @return a real value with the Epps test's statistic.
#'
#' @importFrom MASS ginv
#' @export
#'
#' @author Alicia Nieto-Reyes and Asael Alonzo Matamoros.
#'
#' @seealso \code{\link{lobato.statistic}}
#'
#' @references
#' Epps, T.W. (1987). Testing that a stationary time series is Gaussian. \emph{The
#' Annals of Statistic}. 15(4), 1683-1698.
#'
#' Nieto-Reyes, A., Cuesta-Albertos, J. & Gamboa, F. (2014). A random-projection
#' based test of Gaussianity for stationary processes. \emph{Computational
#' Statistics & Data Analysis, Elsevier}, vol. 75(C), pages 124-141.
#'
#' Press, W.H., Teukolsky, S.A., Vetterling, W.T. and  Flannery, B.P. (2007).
#' Numerical Recipes. The Art of Scientific Computing. \emph{Cambridge
#' University Press}.
#'
#' @examples
#' # Generating an stationary arma process
#' y = arima.sim(100,model = list(ar = 0.3))
#' epps.statistic(y)
#'
epps.statistic =  function(y, lambda = c(1, 2)){
  n = length(y);
  N = 2*length(lambda);
  rn = floor(n^.4)

  lambda_std = lambda/sd(y)

  gm = t(lambda_std %*% t(y))
  gm = rbind(cos(gm), sin(gm))
  gm = matrix(gm, nrow = n, ncol = N)

  gn = apply(gm, 2, mean)

  z = gm - matrix(gn, nrow = n, ncol = N, byrow = TRUE)
  Gm = MASS::ginv(cov(z));

  me = mean(y); sts = sd(y)/sqrt(n); ts2 = sqrt(2/n); Va = sd(y)^2;
  P  = matrix(c(me-sts, me+sts, me, Va*(1-ts2), Va*(1-ts2), Va*(1+ts2)),
              nrow = 3, ncol = 2)
  Y  = cbind(rep(0, 3))

  for (i in 1:3) Y[i] = Quadratic(P[i,], gn, lambda_std, Gm)

  return(amoebam(P, Y, n, gn, lambda_std, Gm))
}
#' The Sieve Bootstrap Epps and Pulley test for normality.
#'
#' Performs the approximated Epps and Pulley's test of normality for univariate time series.
#' Computes the p-value using Psaradakis and Vavra's (2020) sieve bootstrap procedure.
#'
#' @usage epps_bootstrap.test(y, lambda = c(1,2), reps = 500, h = 100, seed = NULL)
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a stationary
#' time series.
#' @param lambda a numeric vector for evaluating the characteristic function.
#' @param reps an integer with the total bootstrap repetitions.
#' @param reps an integer with the total bootstrap repetitions.
#' @param h an integer with the first \code{burn-in} sieve bootstrap replicates.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#'
#' @return A list with class \code{"h.test"} containing the following components:
#'  \item{statistic:}{the sieve bootstrap Epps and Pulley's statistic.}
#'  \item{p.value:}{the p value for the test.}
#'  \item{alternative:}{a character string describing the alternative hypothesis.}
#'  \item{method:}{a character string \dQuote{Sieve-Bootstrap Epps' test}.}
#'  \item{data.name:}{a character string giving the name of the data.}
#'
#' @details
#' The Epps test minimize the process' empirical characteristic function using a
#' quadratic loss in terms of the process two first moments, \emph{Epps, T.W. (1987)}.
#' Approximates the p-value using a sieve-bootstrap procedure \emph{Psaradakis, Z.
#' and Vávra, M. (2020)}.
#'
#' @export
#'
#' @author Asael Alonzo Matamoros and Alicia Nieto-Reyes.
#'
#' @seealso \code{\link{lobato.statistic}}, \code{\link{epps.test}}
#'
#' @references
#' Psaradakis, Z. and Vávra, M. (2020) Normality tests for dependent
#' data: large-sample and bootstrap approaches. Communications in
#' \emph{Statistics-Simulation and Computation 49 (2)}. ISSN 0361-0918.
#'
#' Nieto-Reyes, A., Cuesta-Albertos, J. & Gamboa, F. (2014). A random-projection
#' based test of Gaussianity for stationary processes. \emph{Computational
#' Statistics & Data Analysis, Elsevier}, vol. 75(C), pages 124-141.
#'
#' Epps, T.W. (1987). Testing that a stationary time series is Gaussian. \emph{The
#' Annals of Statistic}. 15(4), 1683-1698.
#'
#' @examples
#' # Generating an stationary arma process
#' y = arima.sim(300, model = list(ar = 0.3))
#' epps_bootstrap.test(y, reps = 1000)
#'
epps_bootstrap.test = function(y, lambda = c(1, 2), reps = 500, h = 100,
                               seed = NULL){

  if( !is.numeric(y) & !is(y,class2 = "ts") )
    stop("y object must be numeric or a time series")

  if( anyNA(y) )
    stop("The time series contains missing values")

  # checking stationarity
  cc = uroot.test(y)
  if(!cc$stationary)
    warning("y has a unit root, vavra.test requires stationary process")

  # checking seasonality
  if(frequency(y) > 1){
    cc = seasonal.test(y)
    if(cc$seasonal)
      warning("y has a seasonal unit root, vavra.test requires stationary process")
  }

  if (!is.null(seed))
    set.seed(seed)

  htest = vavra.test(y = y, lambda = lambda, normality = "epps",
                     reps = reps, h = h, seed = seed)

  return(htest)
}
#'
#' @noRd
#'
Quadratic = function(par, gn, lambda, Gm){
  N = 2*length(lambda);
  ml = par[1]*lambda
  e = exp(-par[2]*(lambda)^2/2)

  g0 = rbind(e*cos(ml),e*sin(ml))
  g0 = as.numeric(g0)

  target = (gn - g0) %*% Gm %*% cbind(gn - g0)
  as.numeric(target)
}
#'
#' @noRd
#'
amoebam = function(P, Y, n, gn, lambda, Gm){
  NDIM = 2; FTOL = .0001; NMAX = 20; ALPHA = 1;
  BETA = 0.5; GAMMA = 2; ITMAX = 500; PR = rep(0,NMAX);
  PRR = rep(0,NMAX); MPTS = NDIM+1; ITER = 0; IXXX = 1;
  while (IXXX==1) {
    ILO = 1; if(Y[1]>Y[2]){ IHI = 1; INHI = 2;} else { IHI = 2; INHI = 1; }
    for (I in 1:MPTS) {
      if(Y[I]<Y[ILO]){ ILO = I }
      if(Y[I]>Y[IHI]){ INHI = IHI; IHI = I; }
      else if(Y[I]>Y[INHI] && I != IHI){ INHI = I }
    }
    RTOL = 2*abs(Y[IHI]-Y[ILO])/(abs(Y[IHI])+abs(Y[ILO]))
    if(!is.na(RTOL) && RTOL<FTOL) { return(n*min(Y)) }
    if(ITER==ITMAX) { return(n*min(Y)) }
    ITER = ITER+1; PBAR = rep(0,NDIM);
    for (I in 1:MPTS) {
      if(I != IHI){ for (J in 1:NDIM) { PBAR[J] = PBAR[J]+P[I,J] } }
    }
    for (J in 1:NDIM) {
      PBAR[J] = PBAR[J]/NDIM
      PR[J] = (1.+ALPHA)*PBAR[J]-ALPHA*P[IHI,J]
    }
    YPR = Quadratic(PR,gn,lambda,Gm)
    if (YPR<=Y[ILO]){
      for (J in 1:NDIM) { PRR[J] = GAMMA*PR[J]+(1.-GAMMA)*PBAR[J] }
      YPRR = Quadratic(PRR,gn,lambda,Gm)
      if (YPRR < Y[ILO]){
        for (J in 1:NDIM) { P[IHI,J] = PRR[J] }
        Y[IHI] = YPRR
      }
      else {
        for (J in 1:NDIM) { P[IHI,J] = PR[J] }; Y[IHI] = YPR;
      }
    }
    else if(YPR >= Y[INHI]){
      if(YPR < Y[IHI]){
        for (J in 1:NDIM) { P[IHI,J] = PR[J] }
        Y[IHI] = YPR
      }
      for (J in 1:NDIM) { PRR[J] = BETA*P[IHI,J]+(1.-BETA)*PBAR[J] }
      YPRR = Quadratic(PRR,gn,lambda,Gm)
      if(YPRR<Y[IHI]){
        for (J in 1:NDIM) { P[IHI,J] = PRR[J] }
        Y[IHI] = YPRR
      }
      else{
        for (I in 1:MPTS) {
          if(I != ILO){
            for (J in 1:NDIM) {
              PR[J] = 0.5*(P[I,J]+P[ILO,J]); P[I,J] = PR[J];
            }
            Y[I] = Quadratic(PR,gn,lambda,Gm)
          }
        }
      }
    } else{
      for (J in 1:NDIM) { P[IHI,J] = PR[J] }; Y[IHI] = YPR;
    }
  }
}
