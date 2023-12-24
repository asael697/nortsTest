#' The Sieve Bootstrap Jarque-Bera test for normality.
#'
#' Performs the Jarque Bera test approximating the p-value
#' using Psaradakis and Vavra's procedure. The null
#' hypothesis (H0), is that the given data follows a Gaussian process.
#'
#' @usage jb_bootstrap.test(y, reps = 1000, h = 100, seed = NULL)
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a
#' stationary time series.
#' @param reps an integer with the total bootstrap repetitions.
#' @param h an integer with the first \code{burn-in} sieve bootstrap replicates.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#'
#' @return A list with class \code{"h.test"} containing the following components:
#' \itemize{
#'  \item{statistic }{the sieve bootstrap Jarque Beras' statistic.}
#'  \item{p.value }{the p value for the test.}
#'  \item{alternative }{a character string describing the alternative hypothesis.}
#'  \item{method }{a character string \dQuote{Sieve-Bootstrap Jarque Beras' test}.}
#'  \item{data.name }{a character string giving the name of the data.}
#' }
#'
#' @details
#' Employs Jarque Beras' skewness-kurtosis test approximating the p-value using
#' a sieve-bootstrap procedure, \emph{Psaradakis, Z. and Vávra, M. (2020)}.
#'
#' @export
#'
#' @author Asael Alonzo Matamoros.
#'
#' @seealso \code{\link{vavra.test}},\code{\link{sieve.bootstrap}}
#'
#' @references
#' Psaradakis, Z. and Vávra, M. (2020) Normality tests for dependent
#' data: large-sample and bootstrap approaches. Communications in
#' \emph{Statistics-Simulation and Computation 49 (2)}. ISSN 0361-0918.
#'
#' Bulmann, P. (1997). Sieve Bootstrap for time series. \emph{Bernoulli}.
#' 3(2), 123 -148.
#'
#' J. B. Cromwell, W. C. Labys and M. Terraza (1994): Univariate Tests for Time
#' Series Models, Sage, Thousand Oaks, CA, pages 20–22.
#'
#' @examples
#' # Generating an stationary arma process
#' y = arima.sim(100,model = list(ar = 0.3))
#' jb_bootstrap.test(y)
#'
jb_bootstrap.test = function(y, reps = 1000, h = 100, seed = NULL){

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

  htest = vavra.test(y = y, c = c, normality = "jb",
                     reps = reps, h = h, seed = seed)

  htest$method = "Sieve-Bootstrap Jarque-Bera test"

  return(htest)
}
