#' The Lobato and Velasco's Test for normality.
#'
#' Performs the Lobato and Velasco's test for normality. The null hypothesis (H0),
#' is that the given data follows a Gaussian process.
#'
#' @usage lobato.test(y,c = 1)
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a stationary time series.
#' @param c a positive real value that identifies the total amount of values used in the
#' cumulative sum.
#'
#' @return A list with class \code{"h.test"} containing the following components:
#' \itemize{
#'  \item{statistic }{the Lobato and Velasco's statistic.}
#'  \item{parameter }{the test degrees freedoms.}
#'  \item{p.value }{the p-value for the test.}
#'  \item{alternative }{a character string describing the alternative hypothesis.}
#'  \item{method }{a character string \dQuote{Lobato and Velasco's test}.}
#'  \item{data.name }{a character string giving the name of the data.}
#' }
#'
#' @details
#'
#' This test  proves a normality assumption in correlated data employing the
#' skewness-kurtosis test statistic, but studentized by standard error estimates
#' that are consistent under serial dependence of the observations. The test was
#' proposed by \emph{Lobato, I., & Velasco, C. (2004)} and implemented by
#' \emph{Nieto-Reyes, A., Cuesta-Albertos, J. & Gamboa, F. (2014)}.
#'
#' @export
#'
#' @author Asael Alonzo Matamoros and Alicia Nieto-Reyes.
#'
#' @seealso \code{\link{lobato.statistic}},\code{\link{epps.test}}
#'
#' @references
#' Lobato, I., & Velasco, C. (2004). A simple test of normality in time series.
#' \emph{Journal of econometric theory}. 20(4), 671-689.
#'
#' Nieto-Reyes, A., Cuesta-Albertos, J. & Gamboa, F. (2014). A random-projection
#' based test of Gaussianity for stationary processes. \emph{Computational
#' Statistics & Data Analysis, Elsevier}, vol. 75(C), pages 124-141.
#'
#' @examples
#' # Generating an stationary arma process
#' y = arima.sim(100,model = list(ar = 0.3))
#' lobato.test(y)
#'
lobato.test = function(y, c = 1){

  if( !is.numeric(y) & !is(y,class2 = "ts") )
    stop("y object must be numeric or a time series")

  if( anyNA(y) )
    stop("The time series contains missing values")

  # checking stationarity
  cc = uroot.test(y)
  if(!cc$stationary)
    warning("y has a unit root, lobato.test requires stationary process")

  # checking seasonality
  if(frequency(y) > 1){
    cc = seasonal.test(y)
    if(cc$seasonal)
      warning("y has a seasonal unit root, lobato.test requires stationary process")
  }

  dname = deparse(substitute(y))
  alt = paste(dname,"does not follow a Gaussian Process")
  tstat = lobato.statistic(y, c = c)
  names(tstat) <- "lobato"
  df =  2
  names(df) <- "df"
  pval = pchisq(q = tstat,df = df,lower.tail = FALSE)

  rval <- list(statistic = tstat,
               parameter = df,
               p.value = pval,
               alternative = alt,
               method = "Lobato and Velasco's test",
               data.name = dname)

  class(rval) <- "htest"
  return(rval)
}
#'  Computes the Lobato and Velasco statistic.
#'
#' Computes the Lobato and Velasco's statistic. This test  proves a normality
#' assumption in correlated data employing the skewness-kurtosis test statistic,
#' but studentized by standard error estimates that are consistent under serial
#' dependence of the observations.
#'
#' @usage lobato.statistic(y,c = 1)
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a stationary time series.
#' @param c a positive real value that identifies the total amount of values used in the
#' cumulative sum.
#'
#' @details This function is the equivalent of \code{GestadisticoVn} of \emph{Nieto-Reyes, A.,
#' Cuesta-Albertos, J. & Gamboa, F. (2014)}.
#'
#' @return A real value with the Lobato and Velasco test's statistic.
#'
#' @importFrom stats arima
#' @export
#'
#' @author Alicia Nieto-Reyes and Asael Alonzo Matamoros.
#'
#' @seealso \code{\link{epps.statistic}}
#'
#' @references
#' Lobato, I., & Velasco, C. (2004). A simple test of normality in time series.
#' \emph{Journal of econometric theory}. 20(4), 671-689.
#'
#' Nieto-Reyes, A., Cuesta-Albertos, J. & Gamboa, F. (2014). A random-projection
#' based test of Gaussianity for stationary processes. \emph{Computational
#' Statistics & Data Analysis, Elsevier}, vol. 75(C), pages 124-141.
#'
#' @examples
#' # Generating an stationary arma process
#' y = arima.sim(100,model = list(ar = 0.3))
#' lobato.statistic(y,3)
#'
lobato.statistic = function(y,c = 1){

  if( !is.numeric(y) & !is(y,class2 = "ts") )
    stop("y object must be numeric or a time series")

  if( anyNA(y) )
    stop("The time series contains missing values")

  # Identifiying the process sample statistics
  n = length(y)
  mu1 = mean(y)
  mu2 = var(y)*(n-1)/n
  mu3 = sum((y-mu1)^3)/n
  mu4 = sum((y-mu1)^4)/n

  hn = ceiling(c*sqrt(n)-1)
  a = stats::acf(y,lag.max = hn,plot = FALSE)
  gamma = as.numeric(a$acf)[2:hn]
  gat = rev(gamma)

  F3 = abs(2*sum(gamma*(gamma+gat)^2)+mu2^3)
  F4 = abs(2*sum(gamma*(gamma+gat)^3)+mu2^4)
  G = n*(mu3^2/(6*F3)+(mu4-3*mu2^2)^2/(24*F4))

  return(G)
}
