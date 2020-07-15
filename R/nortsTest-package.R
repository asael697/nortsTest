#' 'Assessing Normality of Stationary Process.'
#'
#' @docType package
#' @name nortsTest-package
#' @aliases nortsTest
#'
#' @description
#' Despite that several tests for normality in stationary processes have been proposed
#' in the literature, consistent implementations of these tests in programming languages
#' are limited. Four normality test are implemented. The Lobato and Velasco's, Epps,
#' Psaradakis and  Vavra, and the random projections tests for stationary process.
#' Some other diagnostics such as, unit root test for stationarity, seasonal tests for
#' seasonality, and arch effect test for volatility; are also performed. The package also
#' offers residual diagnostic for linear time series models developed in several packages.
#'
#' @details
#' We present four main functions, for testing the hypothesis of
#' normality in stationary process, the \code{epps.test}, \code{lobato.test},
#' \code{rp.test}, and \code{varvra.test}. Additionally, we provide functions
#' for unit root, seasonality and ARCH effects tests for stationary, and other additional
#' methods for visual checks using the \pkg{ggplot2} and \pkg{forecast} packages.
#'
#' @import methods ggplot2 gridExtra forecast nortest stats tseries uroot MASS
#'
#'
#' @references
#' Epps, T.W. (1987). Testing that a stationary time series is Gaussian. \emph{The
#' Annals of Statistic}. 15(4), 1683-1698.\url{https://projecteuclid.org/euclid.aos/1176350618}.
#'
#' Lobato, I., & Velasco, C. (2004). A simple test of normality in time series.
#' \emph{Journal of econometric theory}. 20(4), 671-689.
#' \code{doi:https://doi.org/10.1017/S0266466604204030}.
#'
#' Psaradakis, Z. & Vavra, M. (2017). A distance test of normality for a wide class
#' of stationary process. \emph{Journal of Econometrics and Statistics}. 2, 50-60.
#' \code{doi:https://doi.org/10.1016/j.ecosta.2016.11.005}
#'
#' Nieto-Reyes, A., Cuesta-Albertos, J. & Gamboa, F. (2014). A random-projection
#' based test of Gaussianity for stationary processes. \emph{Computational
#' Statistics & Data Analysis, Elsevier}, vol. 75(C), pages 124-141.
#'
#' Hyndman, R. & Khandakar, Y. (2008). Automatic time series forecasting: the
#' forecast package for \code{R}. \emph{Journal of Statistical Software}. 26(3),
#' 1-22.\code{doi:	10.18637/jss.v027.i03}.
#'
#' Wickham, H. (2008). ggplot2: Elegant Graphics for Data Analysis.
#' \emph{Springer-Verlag New York}.
#'
NULL
