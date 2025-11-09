#' Assessing Normality of a Stationary Process
#'
#' Despite several tests for normality in stationary processes being proposed
#' in the literature, consistent implementations in programming languages are limited.
#' This package implements seven normality tests: the asymptotic Lobato and Velasco,
#' asymptotic Epps, Psaradakis & Vávra, Lobato & Velasco's sieve bootstrap approximation,
#' Elbouch et al., Epps sieve bootstrap approximation, and the random projections test
#' for univariate stationary processes. Additional diagnostics such as unit root tests,
#' seasonality tests, and ARCH effect tests are also included. The Elbouch test additionally
#' performs bivariate normality testing. Residual diagnostics for linear time series models
#' from several packages are also available.
#'
#' @details
#' Functions provided for univariate normality tests include:
#' \code{epps.test}, \code{lobato.test}, \code{rp.test},
#' \code{lobato-bootstrap.test}, \code{epps-bootstrap.test},
#' \code{elbouch.test}, and \code{varvra.test}. The \code{elbouch.test} function can
#' perform a bivariate test if a second time series is supplied. Model diagnostics
#' functions include unit root, seasonality, and ARCH effect tests. Visual checks can
#' be done with \pkg{ggplot2} and \pkg{forecast}.
#'
#' @import methods ggplot2 gridExtra forecast nortest stats tseries uroot MASS
#'
#' @references
#' Epps, T.W. (1987). Testing that a stationary time series is Gaussian. \emph{The
#' Annals of Statistic}, 15(4), 1683-1698. \url{https://projecteuclid.org/euclid.aos/1176350618}
#'
#' Lobato, I., & Velasco, C. (2004). A simple test of normality in time series.
#' \emph{Journal of Econometric Theory}, 20(4), 671-689.
#' \code{doi:https://doi.org/10.1017/S0266466604204030}
#'
#' Psaradakis, Z. & Vávra, M. (2017). A distance test of normality for a wide class
#' of stationary processes. \emph{Journal of Econometrics and Statistics}, 2, 50-60.
#' \code{doi:https://doi.org/10.1016/j.ecosta.2016.11.005}
#'
#' Nieto-Reyes, A., Cuesta-Albertos, J., & Gamboa, F. (2014). A random-projection
#' based test of Gaussianity for stationary processes. \emph{Computational
#' Statistics & Data Analysis}, 75(C), 124-141.
#'
#' Hyndman, R. & Khandakar, Y. (2008). Automatic time series forecasting: the
#' forecast package for \code{R}. \emph{Journal of Statistical Software}, 26(3), 1-22.
#' \code{doi:10.18637/jss.v027.i03}
#'
#' Wickham, H. (2008). ggplot2: Elegant Graphics for Data Analysis.
#' \emph{Springer-Verlag New York}.
"_PACKAGE"