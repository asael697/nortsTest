#' The Psaradakis and  Vávra test for normality.
#'
#' Performs the Psaradakis and Vávra distance test for normality. The null
#' hypothesis (H0), is that the given data follows a Gaussian process.
#'
#' @usage vavra.test(y, normality = c("ad","lobato","jb","cvm","epps"),
#'                   reps = 1000, h = 100, seed = NULL, c = 1, lambda = c(1,2))
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a
#' stationary time series.
#' @param normality A character string naming the desired test for checking
#' normality. Valid values are \code{"epps"} for the Epps, \code{"lobato"} for
#' Lobato and Velasco's, \code{"jb"} for the Jarque and Bera, \code{"ad"} for
#' Anderson Darling test, and \code{"cvm"} for the Cramer Von Mises' test.
#' The default value is \code{"ad"} test.
#' @param reps an integer with the total bootstrap repetitions.
#' @param h an integer with the first \code{burn-in} sieve bootstrap replicates.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#' @param c a positive real value used as argument for the Lobato's test.
#' @param lambda a numeric vector used as argument for the Epps's test.
#'
#' @return A list with class \code{"h.test"} containing the following components:
#'  \item{statistic:}{the sieve bootstrap A statistic.}
#'  \item{p.value:}{the p value for the test.}
#'  \item{alternative:}{a character string describing the alternative hypothesis.}
#'  \item{method:}{a character string \dQuote{Psaradakis and Vávra test}.}
#'  \item{data.name:}{a character string giving the name of the data.}
#'
#' @details
#' The Psaradakis and Vávra test approximates the empirical distribution
#' function of the Anderson Darling's statistic, using a sieve bootstrap
#' approximation. The test was proposed by \emph{Psaradakis, Z. & Vávra, M.
#' (20.17)}.
#'
#' @export
#'
#' @author Asael Alonzo Matamoros.
#'
#' @seealso \code{\link{lobato.test}}, \code{\link{epps.test}}
#'
#' @references
#' Psaradakis, Z. and Vávra, M. (2020) Normality tests for dependent
#' data: large-sample and bootstrap approaches. Communications in
#' \emph{Statistics-Simulation and Computation 49 (2)}. ISSN 0361-0918.
#'
#' Psaradakis, Z. & Vávra, M. (2017). A distance test of normality for a wide class
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
vavra.test = function(y, normality = c("ad","lobato","jb","cvm","epps"),
                      reps = 1000, h = 100, seed = NULL, c = 1, lambda = c(1,2)){

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

  if(!is.null(seed))
    set.seed(seed)

  normality = match.arg(normality)

  ad0 = norm.stat(y, normality = normality, c = c, lambda = lambda)

  ad = vavra.sample(y = as.numeric(y), normality = normality,
                    reps = reps,seed = seed,h = h)

  # Additional values
  dname = deparse(substitute(y))
  alt = paste(dname,"does not follow a Gaussian Process")
  # Bootstrap test statistic
  tstat = ad0
  names(tstat) = paste("bootstrap",normality,sep = "-")

  mtd =ifelse(normality == "ad",
              "Psaradakis-Vavra test",
              paste("Sieve-Bootstrap",normality,"test"))

  # Bootstrap p.value
  pval = mean(ad > ad0)

  rval = list(statistic = tstat,
              p.value = pval,
              alternative = alt,
              method = mtd,
              data.name = dname)
  class(rval) = "htest"
  return(rval)
}
#' Vávra test's sieve Bootstrap sample for Anderson Darling statistic
#'
#' Generates a sieve bootstrap sample of the Anderson-Darling
#' statistic test.
#'
#' @usage vavra.sample(y, normality = c("ad","lobato","jb","cvm","shapiro","epps"),
#'                     reps = 1000, h = 100, seed = NULL, c = 1, lambda = c(1,2))
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a stationary
#' time series.
#' @param normality A character string naming the desired test for checking normality.
#' Valid values are \code{"epps"} for the Epps, \code{"lobato"} for Lobato and Velasco's,
#' \code{"jb"} for the Jarque and Bera, \code{"ad"} for Anderson Darling test,\code{"cvm"}
#' for the Cramer Von Mises' test, and \code{"shapiro"} for the Shapiro's test.
#' The default value is \code{"ad"} test.
#' @param reps an integer with the total bootstrap repetitions.
#' @param h an integer with the first \code{burn-in} sieve bootstrap replicates.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#' @param c a positive real value used as argument for the Lobato's test.
#' @param lambda a numeric vector used as argument for the Epps's test.
#'
#' @details
#' The Vávra test approximates the empirical distribution function of the
#' Anderson-Darlings statistic, using a sieve bootstrap approximation.
#' The test was proposed by \emph{Psaradakis, Z. & Vávra, M (20.17)}.
#'
#' This function is the equivalent of \code{xarsieve} of
#' \emph{Psaradakis, Z. &  Vávra, M (20.17)}.
#'
#' @return A numeric array with the Anderson Darling sieve bootstrap sample
#'
#' @export
#'
#' @author Asael Alonzo Matamoros.
#'
#' @seealso \code{\link{epps.statistic}}, \code{\link{lobato.statistic}}
#'
#' @references
#' Psaradakis, Z. and Vávra, M. (2020) Normality tests for dependent
#' data: large-sample and bootstrap approaches. Communications in
#' \emph{Statistics-Simulation and Computation 49 (2)}. ISSN 0361-0918.
#'
#' Psaradakis, Z. & Vávra, M. (2017). A distance test of normality for a wide class
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
vavra.sample = function(y, normality = c("ad","lobato","jb","cvm","shapiro","epps"),
                        reps = 1000, h = 100, seed = NULL, c = 1, lambda = c(1,2)){

  if(!is.numeric(y) & !is(y,class2 = "ts"))
    stop("y object must be numeric or a time series")

  if(anyNA(y))
    stop("The time series contains missing values")

  if(!is.null(seed))
    set.seed(seed)

  normality = match.arg(normality)
  n = length(y)

  yrep = sieve.bootstrap(y = as.numeric(y),reps = reps,seed = seed,h = h)
  adb = apply(yrep, 1, norm.stat, normality = normality, c = c, lambda = lambda)
  return(adb)
}
#' Generates a sieve bootstrap sample.
#'
#' The function generates a sieve bootstrap sample for a univariate linear
#' stochastic process.
#'
#' @usage sieve.bootstrap(y,reps = 1000,pmax = NULL,h = 100,seed = NULL)
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a
#' stationary time series.
#' @param reps an integer with the total bootstrap repetitions.
#' @param pmax an integer with the max considered lags for the generated
#' \code{ar(p)} process. By default, \code{pmax = NULL}.
#' @param h an integer with the first \code{burn-in} sieve bootstrap replicates.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#'
#' @return A matrix or \code{reps} row and \code{n} columns, with the sieve
#' bootstrap sample and \code{n} the time series length.
#'
#' @details
#' simulates bootstrap samples for the stochastic process y, using a stationary
#' auto-regressive model of order \code{"pmax"}, \code{AR(pmax)}. If
#' \code{pmax = NULL} (\emph{default}), the function estimates the process maximum
#' lags using an \code{AIC} as a model selection criteria.
#'
#' @importFrom forecast auto.arima
#' @importFrom stats arima
#' @export
#'
#' @author Asael Alonzo Matamoros.
#'
#' @seealso \code{\link{lobato.test}}, \code{\link{epps.test}}.
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

  if(!is.numeric(y) & !is(y,class2 = "ts"))
    stop("y object must be numeric or a time series")

  if( anyNA(y) )
    stop("The time series contains missing values")

  if(!is.null(seed))
    set.seed(seed)

  n = length(y)
  if(is.null(pmax)){
    pmax = floor(log(n)^2)
    mod = forecast::auto.arima(y,d = 0,D = 0,max.p = pmax,
                               allowmean = TRUE,max.q = 0,
                               max.P = 0,max.Q = 0)
    p = length(mod$model$phi)
    phi = mod$model$phi
  }
  else{
    p = ifelse(floor(log(n)^2) > pmax, pmax, floor(log(n)^2))

    mod = forecast::auto.arima(y,d = 0,D = 0,max.p = pmax,
                                 allowmean = TRUE,max.q = 0,
                                 max.P = 0,max.Q = 0)

    p = length(mod$model$phi)
    phi = mod$model$phi
  }

  if(p <= 0){
    p = 1
    mod = stats::arima(y,order = c(p,0,0))

    p = length(mod$model$phi)
    phi = mod$model$phi
  }

  N = n + h
  sig = sd(mod$residuals)/sqrt(n-2*p-1)

  x <- parallel::mclapply(1:reps, mc.set.seed = is.null(seed), FUN = function(i){
    N = n + h
    ar_boot <- arima.sim(n = N,model = list(ar = phi,sd = sig))
    c(i,as.vector(ar_boot[(h+1):N]))
  })

  for_err <-matrix(unlist(x),ncol = n+1,byrow = TRUE)

  return(for_err[,-1])
}
#' Normality statistics
#'
#' @importFrom tseries jarque.bera.test
#' @importFrom stats shapiro.test
#' @importFrom nortest ad.test
#' @importFrom nortest cvm.test
#' @noRd
#'
norm.stat = function(y, normality = c("ad","lobato","jb","cvm","shapiro","epps"),
                     c = 1, lambda = c(1, 2)){

  normality = match.arg(normality)

  if(normality == "ad"){
    cc = suppressWarnings(nortest::ad.test(y)$statistic)
  }
  else if(normality == "lobato"){
    cc = suppressWarnings(lobato.statistic(y,c = c))
  }
  else if(normality == "jb"){
    cc = suppressWarnings(tseries::jarque.bera.test(y)$statistic)
  }
  else if(normality == "cvm"){
    cc = suppressWarnings(nortest::cvm.test(y)$statistic)
  }
  else if(normality == "shapiro"){
    cc = suppressWarnings(stats::shapiro.test(y)$statistic)
  }
  else{
    cc = suppressWarnings(epps.statistic(y,lambda = lambda))
  }

  return(cc)
}
