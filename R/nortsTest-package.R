#' The 'nortsTest' package.
#'
#' @docType package
#' @name nortsTest-package
#' @aliases nortsTest
#'
#' @description
#' This package is for testing normal distribution of stationary time series,
#' it performs the Lobato and Velasco's, Epps, Psaradakis and  Vavra, and the
#' random projection tests for stationary univariate time  series. Some other
#' diagnostics such as unit root test for stationarity, seasonal tests for
#' seasonality, and ARCH effect tests for volatility; are also performed. The
#' package also offers residual diagnostic for linear time series model developed
#' in several time series analysis packages.
#'
#' @details
#' The \pkg{nortsTest} package has four main functions, for testing the hypothesis of
#' normal distribution in stationary time series, the \code{epps.test}, \code{lobato.test},
#' \code{rp.test}, and \code{varvra.test()}. Additionally, we provide functions
#' for unit root, seasonality and ARCH effects tests for stationary, and other additional
#' methods for visual checks using the \pkg{ggplot2} and \pkg{forecast} packages.
#'
#' @import methods ggplot2 gridExtra forecast nortest stats tseries uroot MASS
#'
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
#' \url{http://www.sciencedirect.com/science/article/pii/S2452306216300296"}.
#' \code{doi: https://doi.org/10.1016/j.ecosta.2016.11.005}.
#'
#' Nieto-Reyes, A., Cuesta-Albertos, J. & Gamboa, F. (2014). A random-projection
#' based test of Gaussianity for stationary processes. \emph{Computational
#' Statistics & Data Analysis, Elsevier}, vol. 75(C), pages 124-141.
#' \url{http://www.sciencedirect.com/science/article/pii/S0167947314000243}
#' \code{doi:https://doi.org/10.1016/j.csda.2014.01.013}.
#'
#' Hyndman, R. & Khandakar, Y. (2008). Automatic time series forecasting: the
#' forecast package for \code{R}. \emph{Journal of Statistical Software}. 26(3),
#' 1-22.\url{http://www.jstatsoft.org/article/view/v027i03}.
#'
#' Wickham, H. (2008). ggplot2: Elegant Graphics for Data Analysis.
#' \emph{Springer-Verlag New York} \url{https://ggplot2.tidyverse.org}
#'
NULL
