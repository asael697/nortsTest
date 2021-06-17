#' The Psaradakis and  Vavra test for normality
#'
#' Performs the Psaradakis and Vavra distance test for normality. The null hypothesis (H0),
#' is that the given data follows a Gaussian process.
#'
#' @usage  vavra.test(y,reps = 1000,h = 100,seed = NULL)
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a stationary time series.
#' @param reps an integer with the total bootstrap repetitions.
#' @param h an integer with the first \code{burn-in} sieve bootstrap replicates.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#'
#' @return a h.test class with the main results of the Epps hypothesis test. The
#' h.test class have the following values:
#' \itemize{
#'  \item{"bootstrap A"}{The sieve bootstrap A statistic}
#'  \item{"p.value"}{The p value}
#'  \item{"alternative"}{The alternative hypothesis}
#'  \item{"method"}{The used method}
#'  \item{"data.name"}{The data name.}
#' }
#'
#' @details
#' The Psaradakis and Vavra test approximates the empirical distribution
#' function of the Anderson Darling's statistic, using a sieve bootstrap
#' approximation. The test was proposed by \emph{Psaradakis, Z. & Vavra, M (20.17)}.
#'
#' @export
#'
#' @author Asael Alonzo Matamoros.
#'
#' @seealso \code{\link{lobato.test}},\code{\link{epps.test}}
#'
#' @references
#' Psaradakis, Z. & Vavra, M. (2017). A distance test of normality for a wide class
#' of stationary process. \emph{Journal of Econometrics and Statistics}. 2, 50-60.
#'
#' Bulmann, P. (1997). Sieve Bootstrap for time series. \emph{Bernoulli}.
#' 3(2), 123 -148.
#'
#' @examples
#' # Generating an stationary arma process
#' y = arima.sim(100,model = list(ar = 0.3))
#' vavra.test(y)
#'
vavra.test = function(y,reps = 1000,h = 100,seed = NULL){

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

  ad0 = ad.statistic(y)

  ad = vavra.sample(y = as.numeric(y),reps = reps,seed = seed,h = h)

  # Additional values
  dname = deparse(substitute(y))
  alt = paste(dname,"does not follow a Gaussian Process")
  # Bootstrap test statistic
  tstat = mean(ad)
  names(tstat) = "bootstrap A"

  # Bootstrap p.value
  pval = mean(ad > ad0)

  rval = list(statistic =tstat, p.value = pval,
              alternative = alt,
              method = "Psaradakis-Vavra test", data.name = dname)
  class(rval) = "htest"
  return(rval)
}
#' vavra test's sieve Bootstrap sample for Anderson Darling statistic
#'
#' Generates a sieve bootstrap sample of the Anderson-Darling
#' statistic test.
#'
#' @usage vavra.sample(y,reps = 1000,h = 100,seed = NULL)
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a stationary time series.
#' @param reps an integer with the total bootstrap repetitions.
#' @param h an integer with the first \code{burn-in} sieve bootstrap replicates.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#'
#' @details
#' The Vavra test approximates the empirical distribution function of the
#' Anderson-Darlings statistic, using a sieve bootstrap approximation.
#' The test was proposed by \emph{Psaradakis, Z. & Vavra, M (20.17)}.
#'
#' This function is the equivalent of \code{xarsieve} of
#' \emph{Psaradakis, Z. &  Vavra, M (20.17)}.
#'
#' @return A numeric array with the Anderson Darling sieve bootstrap sample
#'
#' @export
#'
#' @author Asael Alonzo Matamoros.
#'
#' @seealso \code{\link{epps.statistic}} \code{\link{lobato.statistic}}
#'
#' @references
#' Psaradakis, Z. & Vavra, M. (2017). A distance test of normality for a wide class
#' of stationary process. \emph{Journal of Econometrics and Statistics}. 2, 50-60.
#'
#' Bulmann, P. (1997). Sieve Bootstrap for time series. \emph{Bernoulli}.
#' 3(2), 123 -148.
#'
#' @examples
#' # Generating an stationary arma process
#' y = arima.sim(100,model = list(ar = 0.3))
#' adbs = vavra.sample(y)
#' mean(adbs)
#'
vavra.sample = function(y,reps = 1000,h = 100,seed = NULL){

  if( !is.numeric(y) & !is(y,class2 = "ts") )
    stop("y object must be numeric or a time series")

  if( anyNA(y) )
    stop("The time series contains missing values")

  if (!is.null(seed))
    set.seed(seed)

  n = length(y)

  yrep = sieve.bootstrap(y = as.numeric(y),reps = reps,seed = seed,h = h)
  adb = apply(yrep, 1,ad.statistic)
  return(adb)
}
#' Generates a sieve bootstrap sample
#'
#' The function generates a sieve bootstrap sample for a univariate stochastic process.
#'
#' @usage  sieve.bootstrap(y,reps = 1000,pmax = NULL,h = 100,seed = NULL)
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a stationary time series.
#' @param reps an integer with the total bootstrap repetitions.
#' @param pmax an integer with the max considered lags for the generated \code{ar(p)} process.
#' By default, \code{pmax = NULL}.
#' @param h an integer with the first \code{burn-in} sieve bootstrap replicates.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#'
#' @return A matrix or \code{reps} row and \code{n} columns, with the sieve bootstrap
#' sample and \code{n} the time series length.
#'
#' @details
#' simulates bootstrap samples for the stochastic process y, using a stationary
#' auto-regressive model of order \code{"pmax"}, \code{AR(pmax)}. If \code{pmax = NULL} (\emph{default}),
#' the function estimates the process maximum lags using an \code{AIC} as a model
#' selection criteria.
#'
#' @importFrom forecast auto.arima
#' @importFrom stats arima
#' @export
#'
#' @author Asael Alonzo Matamoros
#'
#' @seealso \code{\link{lobato.test}},\code{\link{epps.test}}
#'
#' @references
#' Bulmann, P. (1997). Sieve Bootstrap for time series. \emph{Bernoulli}.
#' 3(2), 123 -148.
#'
#' @examples
#' # Generating an stationary arma process
#' y = arima.sim(100,model = list(ar = 0.3))
#' M = sieve.bootstrap(y)
#'
sieve.bootstrap = function(y,reps = 1000,pmax = NULL,h = 100,seed = NULL){

  if( !is.numeric(y) & !is(y,class2 = "ts") )
    stop("y object must be numeric or a time series")

  if( anyNA(y) )
    stop("The time series contains missing values")

  if (!is.null(seed))
    set.seed(seed)

  n = length(y)
  if(is.null(pmax)){
    pmax = floor(log(n)^2)
    mod = forecast::auto.arima(y,d = 0,D = 0,max.p = pmax,
                               allowmean = TRUE,max.q = 0,
                               max.P = 0,max.Q = 0)
    p = length(mod$coef)
  }
  else{
    if(floor(log(n)^2)> pmax)
      p = pmax
    else
      p = 1
  }

  if(p <= 0) p = 1

  mod = stats::arima(y,order = c(p,0,0),include.mean = TRUE)
  phi = mod$coef
  p = length(phi)

  N = n + h
  sig = sd(mod$residuals)/sqrt(n-2*p-1)
  yrep = matrix(0,nrow = reps,ncol = n)

  for(i in 1:reps){
    yb = y[n:(n-p)]
    for(j in (p+1):N){
      e = rnorm(1,mean = 0,sd = sig)
      yb[j] = sum(phi[1:(p-1)]*yb[(j-p+1):(j-1)]) + phi[p] + e
    }
    yrep[i,] = yb[(h+1):N]
  }
  row.names(yrep) = paste0("repetition",1:reps)
  colnames(yrep)  = paste0("observation",1:n)

  return(yrep)
}
#' Anderson-Darling statistic
#'
#' @importFrom  nortest ad.test
#' @noRd
#'
ad.statistic = function(y){
  m = nortest::ad.test(y)$statistic
  names(m) = NULL
  return(m)
}
