#' The Unit root tests function
#'
#' Perform a unit root test to check stationarity in a linear stochastic process
#'
#' @usage  uroot.test(y,unit_root="adf",alpha=0.05)
#'
#' @param y A univariate time series
#' @param unit_root A character string naming the desired unit root test for checking stationarity.
#' Valid values are \code{"adf"} for the Augmented Dickey-Fuller, \code{"pp"} for the Phillips-Perron,
#' \code{"kpss"} for Kwiatkowski, Phillips, Schmidt, and Shin, and \code{"box"} for the Ljung-Box. The default
#' value is \code{"adf"} for the Augmented Dickey-Fuller test.
#' @param alpha Level of the test, possible values range from 0.01 to 0.1. By default \code{alpha = 0.05}
#' is used
#'
#' @return a h.test class with the main results of unit root hypothesis test.
#'
#' @details Several different tests are available:
#' In the  \code{kpss} test, the null hypothesis that \code{y} has a stationary root
#' against a unit-root alternative. In the two remaining tests, the null hypothesis
#' is that \code{y} has a unit root against a stationary root alternative. By default,
#'  \code{alpha = 0.05} is used to select the more likely hypothesis.
#'
#' @importFrom tseries adf.test pp.test kpss.test
#' @importFrom stats Box.test
#' @export
#'
#' @author  Asael Alonzo Matamoros
#'
#' @seealso \code{\link{normal.test}},\code{\link{seasonal.test}}
#'
#' @references
#' Dickey, D. & Fuller, W. (1979). Distribution of the Estimators for
#' Autoregressive Time Series with a Unit Root. \emph{Journal of the American
#' Statistical Association}. 74, 427-431.
#'
#' Kwiatkowski, D., Phillips, P., Schmidt, P. & Shin, Y. (1992). Testing the Null
#' Hypothesis of Stationarity against the Alternative of a Unit Root,
#' \emph{Journal of Econometrics}. 54, 159-178.
#'
#' Phillips, P. & Perron, P. (1988). Testing for a unit root in time series regression,
#' \emph{Biometrika}. 72(2), 335-346.
#'
#' Ljung, G. M. & Box, G. E. P. (1978). On a measure of lack of fit in time series models.
#' \emph{Biometrika}. 65, 297-303. \code{doi: 10.2307/2335207}.
#'
#' @keywords unit root test, time series, stationary
#'
#' @examples
#' #  stationary  ar process
#' y = arima.sim(100,model = list(ar = 0.3))
#' uroot.test(y)
#'
#' # a random walk process
#' y = cumsum(y)
#' uroot.test(y,unit_root = "pp")
#'
uroot.test = function(y,unit_root = "adf",alpha = 0.05){

  if( !is.numeric(y) & !is(y,class2 = "ts") )
    stop("y object must be numeric or a time series")

  if(NA %in% y )
    stop("The time series contains missing values")

  if(unit_root == "kpss"){
    cc = suppressWarnings(tseries::kpss.test(y))
    cc$stationary = cc$p.value > alpha
    cc$alternative = "unit root"

    if(cc$stationary)
      cc$Conc = "Conclusion: y is stationary"
    else
      cc$Conc = "Conclusion: y has a unit root"
  }
  else if(unit_root == "pp"){
    cc = suppressWarnings(tseries::pp.test(y))
    cc$stationary = cc$p.value < alpha

    if(cc$stationary)
      cc$Conc = "Conclusion: y is stationary"
    else
      cc$Conc = "Conclusion: y has a unit root"
  }
  else if(unit_root == "box"){
    cc = suppressWarnings(stats::Box.test(y,lag = 1,type = "Ljung"))
    cc$alternative = "unit root"
    cc$stationary = cc$p.value > alpha

    if(cc$stationary)
      cc$Conc = "Conclusion: y is stationary"
    else
      cc$Conc = "Conclusion: y has a unit root"
  }
  else{
    cc = suppressWarnings(tseries::adf.test(y))
    cc$stationary = cc$p.value < alpha

    if(cc$stationary)
      cc$Conc = "Conclusion: y is stationary"
    else
      cc$Conc = "Conclusion: y has a unit root"
  }
  return(cc)
}
#' The goodness of fit tests for normal distribution.
#'
#' Perform goodness of fit test to check Gaussian distribution in stationarity stochastic process.
#'
#' @usage  normal.test(y,normality="epps",alpha=0.05)
#'
#' @param y A univariate time series
#' @param normality A character string naming the desired test for checking gaussian distribution.
#' Valid values are \code{"epps"} for the Epps, \code{"lobato"} for Lobato and Velasco's,\code{"vavra"} for
#' the Psaradakis and  Vavra, \code{"rp"} for the random projections, \code{"jb"} for the Jarque and Beras,
#' \code{"ad"} for Anderson Darling test, and \code{"shapiro"} for the Shapiro-Wilk's test. The default value
#' is \code{"epps"} test.
#' @param alpha Level of the test, possible values range from 0.01 to 0.1. By default \code{alpha = 0.05}
#' is used
#'
#' @return An h.test class with the main results of normal hypothesis test.
#'
#' @details Several different tests are available:
#' \code{"lobato"}, \code{"epps"}, \code{"vavras"} and \code{"rp"} test are for testing normal
#' distributions in stationary time series. \code{"jb"}, \code{"ad"}, and  \code{"shapiro"} tests
#' are for numeric data. In all cases, the alternative hypothesis is that \code{y} has
#' a Gaussian distribution. By default, \code{alpha = 0.05} is used to select the more likely
#' hypothesis.
#'
#' @importFrom tseries jarque.bera.test
#' @importFrom stats shapiro.test
#' @importFrom nortest ad.test
#' @export
#'
#' @author  Asael Alonzo Matamoros
#'
#' @seealso \code{\link{uroot.test}},\code{\link{seasonal.test}}
#'
#' @references
#' Epps, T.W. (1987). Testing that a stationary time series is Gaussian. \emph{The
#' Annals of Statistic}. 15(4), 1683-1698.\url{http://www.jstor.org/stable/2336512}.
#' \code{doi:10.1214/aos/1176350618}
#'
#' Lobato, I., & Velasco, C. (2004). A simple test of normality in time series.
#' \emph{Journal of econometric theory}. 20(4), 671-689.
#' \url{doi:10.1017/S0266466604204030}
#'
#' Psaradakis, Z. & Vavra, M. (2017). A distance test of normality for a wide class
#' of stationary process. \emph{Journal of Econometrics and Statistics}. 2, 50-60.
#' \url{http://www.sciencedirect.com/science/article/pii/S2452306216300296}.
#' \code{doi: https://doi.org/10.1016/j.ecosta.2016.11.005}.
#'
#' Nieto-Reyes, A., Cuesta-Albertos, J. & Gamboa, F. (2014). A random-projection
#' based test of Gaussianity for stationary processes. \emph{Computational
#' Statistics & Data Analysis, Elsevier}, vol. 75(C), pages 124-141.
#' \url{http://www.sciencedirect.com/science/article/pii/S0167947314000243}
#' \code{doi:https://doi.org/10.1016/j.csda.2014.01.013}.
#'
#' Patrick, R. (1982). An extension of Shapiro and Wilk's W test for
#' normality to large samples. \emph{Journal of Applied Statistics}.
#' 31, 115-124. \code{doi: 10.2307/2347973}.
#'
#' Cromwell, J. B., Labys, W. C. & Terraza, M. (1994). Univariate Tests for
#' Time Series Models. \emph{Sage, Thousand Oaks, CA}. 20-22.
#'
#' @keywords time series, stationarity, normal distribution
#'
#' @examples
#' #  stationary  ar process
#' y = arima.sim(100,model = list(ar = 0.3))
#' normal.test(y) # epps test
#'
#' # normal random sample
#' y = rnorm(100)
#' normal.test(y,normality = "shapiro")
#'
#' # exponential random sample
#' y = rexp(100)
#' normal.test(y,normality = "ad")
#'
normal.test = function(y,normality = "epps",alpha = 0.05){

  if( !is.numeric(y) & !is(y,class2 = "ts") )
    stop("y object must be numeric or a time series")

  if(NA %in% y )
    stop("The time series contains missing values")

  if(normality == "lobato"){
    cc = suppressWarnings(lobato.test(y))
    cc$Gaussian = cc$p.value > alpha
    if(cc$Gaussian)
      cc$Conc = "Conclusion: y has a normal distribution"
    else
      cc$Conc = "Conclusion: y does not have a Gaussian distribution"
  }
  else if(normality == "vavra"){
    cc = suppressWarnings(vavra.test(y))
    cc$Gaussian = cc$p.value > alpha
    if(cc$Gaussian)
      cc$Conc = "Conclusion: y has a normal distribution"
    else
      cc$Conc = "Conclusion: y does not have a Gaussian distribution"
  }
  else if(normality == "rp"){
    cc = suppressWarnings(rp.test(y))
    cc$Gaussian = cc$p.value > alpha
    if(cc$Gaussian)
      cc$Conc = "Conclusion: y has a normal distribution"
    else
      cc$Conc = "Conclusion: y does not have a Gaussian distribution"
  }
  else if(normality == "jb"){
    cc = suppressWarnings(tseries::jarque.bera.test(y))
    cc$alternative = "y is not Gaussian"
    cc$Gaussian = cc$p.value > alpha
    if(cc$Gaussian)
      cc$Conc = "Conclusion: y has a normal distribution"
    else
      cc$Conc = "Conclusion: y does not have a Gaussian distribution"
  }
  else if(normality == "ad"){
    cc = suppressWarnings(nortest::ad.test(y))
    cc$alternative = "y is not Gaussian"
    cc$Gaussian = cc$p.value > alpha
    if(cc$Gaussian)
      cc$Conc = "Conclusion: y has a normal distribution"
    else
      cc$Conc = "Conclusion: y does not have a Gaussian distribution"
  }
  else if(normality == "shapiro"){
    cc = suppressWarnings(stats::shapiro.test(y))
    cc$alternative = "y is not Gaussian"
    cc$Gaussian = cc$p.value > alpha
    if(cc$Gaussian)
      cc$Conc = "Conclusion: y has a normal distribution"
    else
      cc$Conc = "Conclusion: y does not have a Gaussian distribution"
  }
  else{
    cc = suppressWarnings(epps.test(y))
    cc$Gaussian = cc$p.value > alpha
    if(cc$Gaussian)
      cc$Conc = "Conclusion: y has a normal distribution"
    else
      cc$Conc = "Conclusion: y does not have a Gaussian distribution"
  }
  return(cc)
}
#' The Seasonal unit root tests function
#'
#' Perform a seasonal unit root test to check seasonality in a linear stochastic process
#'
#' @usage  seasonal.test(y,seasonal="ocsb",alpha=0.05)
#'
#' @param y A univariate time series
#' @param seasonal A character string naming the desired seasonal unit root test for checking seasonality.
#' Valid values are \code{"ocsb"} for the Osborn, Chui, Smith, and Birchenhall, \code{"ch"} for the
#' Canova and Hansen, and \code{"hegy"} for Hylleberg, Engle, Granger, and Yoo. The default value is
#' \code{"ocsb"} for the Osborn, Chui, Smith, and Birchenhall test.
#' @param alpha Level of the test, possible values range from 0.01 to 0.1. By default \code{alpha = 0.05}
#' is used
#'
#' @return a h.test class with the main results of unit root hypothesis test.
#'
#' @details Several different tests are available:
#' In the  \code{kpss} test, the null hypothesis that \code{y} has a stationary root
#' against a unit-root alternative. In the two remaining tests, the null hypothesis
#' is that \code{y} has a unit root against a stationary root alternative. By default,
#'  \code{alpha = 0.05} is used to select the more likely hypothesis.
#'
#' @importFrom uroot ch.test hegy.test
#' @importFrom forecast ocsb.test
#' @export
#'
#' @author  Asael Alonzo Matamoros
#'
#' @seealso \code{\link{normal.test}},\code{\link{uroot.test}}
#'
#' @references
#' Osborn, D., Chui, A., Smith, J., & Birchenhall, C. (1988). Seasonality and the
#' order of integration for consumption. \emph{Oxford Bulletin of Economics
#' and Statistics}. 50(4), 361-377.
#'
#' Canova, F. & Hansen, B. (1995). Are Seasonal Patterns Constant over Time? A Test
#' for Seasonal Stability. \emph{Journal of Business and Economic Statistics}.
#' 13(3), 237-252.
#'
#' Hylleberg, S., Engle, R., Granger, C. & Yoo, B. (1990). Seasonal integration
#' and cointegration. \emph{Journal of Econometrics} 44(1), 215-238.
#'
#' @keywords unit root test, time series, seasonality
#'
#' @examples
#' #  stationary  ar process
#' y = ts(rnorm(100),frequency = 6)
#' seasonal.test(y)
#'
seasonal.test = function(y,seasonal = "ocsb",alpha = 0.05){

  if( !is.numeric(y) & !is(y,class2 = "ts") )
    stop("y object must be numeric or a time series")

  if(NA %in% y )
    stop("The time series contains missing values")

  if(seasonal == "ch"){
    cc = suppressWarnings(uroot::ch.test(y))
    cc$seasonal = any(cc$pvalues < 0.05)
    cc$alternative = "seasonality"

    if(cc$seasonal)
      cc$Conc = "Conclusion: y has a seasonal unit root"
    else
      cc$Conc = "Conclusion: y is stationary"
  }
  else if(seasonal == "hegy"){
    cc = suppressWarnings(uroot::hegy.test(y))
    cc$seasonal = any(cc$pvalues > 0.05)
    cc$alternative = "stationarity"

    if(cc$seasonal)
      cc$Conc = "Conclusion: y has a seasonal unit root"
    else
      cc$Conc = "Conclusion: y is stationary"
  }
  else{
    cc = suppressWarnings(forecast::ocsb.test(y,maxlag = 3,lag.method = "AIC"))
    cc$seasonal= cc$statistics > cc$critical

    if(cc$seasonal)
      cc$Conc = "Conclusion: y has a seasonal unit root"
    else
      cc$Conc = "Conclusion: y is stationary"
  }
  return(cc)
}
#' The ARCH effect test function.
#'
#' Performs the Pormanteau Q and Lagrange Multipliers test for homoscedasticity in  a univariate
#' stationary process. The null hypothesis (H0), is that the process is homoscedastic.
#'
#' @usage  arch.test(y,arch="box",alpha=0.05,lag.max = 2)
#'
#' @param y A univariate time series
#' @param arch A character string naming the desired  test for checking stationarity. Valid values are
#' \code{"box"} for the Ljung-Box, and \code{"Lm"} for the Lagrange Multiplier test. The default
#' value is \code{"box"} for the Augmented Ljung-Box test.
#' @param lag.max an integer with the number of used lags.
#' @param alpha Level of the test, possible values range from 0.01 to 0.1. By default \code{alpha = 0.05}
#' is used
#'
#' @return a h.test class with the main results of unit root hypothesis test.
#'
#' @details Several different tests are available:
#' Performs Portmanteau Q and Lagrange Multiplier tests for the null hypothesis that the residuals of
#' an ARIMA model are homoscedastic. The ARCH test is based on the fact that if the residuals (defined
#' as \code{e(t)}) are heteroscedastic, the squared residuals (e^2[t]) are autocorrelated.
#' The first type of test is to examine whether the squares of residuals are a sequence of white noise,
#' which is called the Portmanteau Q test, and similar to the Ljung-Box test on the squared residuals.
#' By default, \code{alpha = 0.05} is used to select the more likely hypothesis.
#'
#' @importFrom stats Box.test
#' @export
#'
#' @author  Asael Alonzo Matamoros
#'
#' @seealso \code{\link{normal.test}},\code{\link{seasonal.test}},\code{\link{uroot.test}}
#'
#' @references
#' Engle, R. F. (1982). Autoregressive Conditional Heteroscedasticity
#' with Estimates of the Variance of United Kingdom Inflation.
#' \emph{Econometrica}. 50(4), 987-1007.
#' \code{doi:https://doi.org/0012-9682(198207)50:4<987:ACHWEO>2.0.CO;2-3}
#'
#' McLeod, A. I. & W. K. Li. (1984). Diagnostic Checking ARMA Time Series
#' Models Using Squared-Residual Autocorrelations. \emph{Journal of Time
#' Series Analysis.} 4, 269-273.
#' \code{doi: https://doi.org/10.1111/j.1467-9892.1983.tb00373.x}
#'
#'
#' @keywords time series, stationary, arch effect, heteroscedastic
#'
#' @examples
#' #  stationary  ar process
#' y = arima.sim(100,model = list(ar = 0.3))
#' arch.test(y)
#'
arch.test = function(y,arch = "box",alpha = 0.05,lag.max = 2){

  if( !is.numeric(y) & !is(y,class2 = "ts") )
    stop("y object must be numeric or a time series")

  if(NA %in% y )
    stop("The time series contains missing values")

  if(arch == "Lm"){
    cc = suppressWarnings(Lm.test(y,lag.max = lag.max,alpha = alpha))
  }
  else{
    cc = suppressWarnings(stats::Box.test(y^2,lag = lag.max,type = "Ljung"))
    cc$alternative = "y is heteroscedastic"
    cc$hetero = cc$p.value < alpha

    if(cc$hetero)
      cc$Conc = "Conclusion: y is heteroscedastic"
    else
      cc$Conc = "Conclusion: y is homoscedastic"
  }
  return(cc)
}
