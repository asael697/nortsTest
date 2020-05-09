#' The Lagrange Multiplier Test for arch effect.
#'
#' Performs the Lagrange Multipliers test for homoscedasticity in  a univariate
#' stationary process. The initial hypothesis (H0), is that Y process is
#' homoscedastic.
#'
#' @usage  Lm.test(y,lag.max = 2,alpha = 0.05)
#'
#' @param y a numeric vector or an object of the ts class containing an univariate
#' time series to be tested.
#' @param lag.max an integer with the number of used lags.
#' @param alpha Level of the test, possible values range from 0.01 to 0.1. By default
#' \code{alpha = 0.05} is used
#'
#' @return a h.test class with the main results of the Lagrage multiplier hypothesis test.
#' The h.test class have the following values
#' \itemize{
#'  \item{"Lm"}{The lagrange multiplier statistic}
#'  \item{"df"}{The test degrees freedoms}
#'  \item{"p.value"}{The p value}
#'  \item{"alternative"}{The alternative hypothesis}
#'  \item{"method"}{The used method}
#'  \item{"data.name"}{The data name}
#' }
#'
#' @details
#' The Lagrange Multiplier test proposed by \emph{Engle (1982)}
#' fits a linear regression model for the squared residuals and
#' examine whether the fitted model is significant. So the null
#' hypothesis is that the squared residuals are a sequence of
#' white noise, namely, the residuals are homoscedastic.
#'
#' @export
#'
#' @author Asael Alonzo Matamoros
#'
#' @seealso \code{\link{arch.test}}
#'
#' @references
#' Engle, R. F. (1982). Autoregressive Conditional Heteroscedasticity
#' with Estimates of the Variance of United Kingdom Inflation.
#' \emph{Econometrica}. 50(4), 987-1007.
#' \code{doi:https://doi.org/0012-9682(198207)50:4<987:ACHWEO>2.0.CO;2-3}
#'
#' McLeod, A. I. and W. K. Li. (1984). Diagnostic Checking ARMA Time Series
#' Models Using Squared-Residual Autocorrelations. \emph{Journal of Time
#' Series Analysis.} 4, 269-273.
#' \code{doi: https://doi.org/10.1111/j.1467-9892.1983.tb00373.x}
#'
#' @examples
#' # generating an stationary arma process
#' y = arima.sim(100,model = list(ar = 0.3))
#' Lm.test(y)
#'
Lm.test = function(y,lag.max = 2,alpha = 0.05){

  if( !is.numeric(y) & !is(y,class2 = "ts") )
    stop("y object must be numeric or a time series")

  if(NA %in% y )
    stop("The time series contains missing values")

  k = lag.max
  if(lag.max < 2) k = 2

  res2 =  y^2
  n = length(res2)
  if (n < 2)
    stop("not enough length of residuals")

  if(length(n) < lag.max )
    k = 2

  dname = deparse(substitute(y))
  rt = embed(res2, k)
  lm.fit = lm(rt[, 1] ~ rt[, -1])
  SSE = sum((rt[, 1] - residuals(lm.fit))^2)
  SST = sum((rt[, 1] - mean(rt[, 1]))^2)
  LM = ((SST - SSE)/k)/(SSE/(n - 2 * k - 1))
  names(LM) = "LM"
  df = k-1
  names(df) = "df"
  pval =  1 - pchisq(LM, df)

  rval = list(statistic =LM, parameter = df, p.value = pval,
              alternative = "y is heteroscedastic",
              method = "Lagrange Multiplier test", data.name = dname)

  rval$hetero = rval$p.value < alpha

  if(rval$hetero)
    rval$Conc = "Conclusion: y is heteroscedastic"
  else
    rval$Conc = "Conclusion: y is homoscedastic"

  class(rval) = "htest"
  return(rval)
}
