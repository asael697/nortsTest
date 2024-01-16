#' Computes El Bouch, et al.'s test for normality of multivariate dependent samples.
#'
#' Computes the El Bouch, Michel, & Comon's  test for normality of a bivariate dependent samples.
#'
#' @usage elbouch.test(y, x = NULL)
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a stationary time series.
#' @param x a numeric vector or an object of the \code{ts} class containing a stationary time series.
#'
#' @return A list with class \code{"h.test"} containing the following components:
#'  \item{statistic:}{the El Bouch Z statistic.}
#'  \item{p.value:}{the p value for the test.}
#'  \item{alternative:}{a character string describing the alternative hypothesis.}
#'  \item{method:}{a character string \dQuote{El Bouch, Michel & Comon's test}.}
#'  \item{data.name:}{a character string giving the name of the data.}
#'
#' @details This function computes El Bouch, et al. (2022) test for normality of
#' bivariate dependent samples. If `x` is set to `NULL`, the test computes the univariate
#' counterpart. This test is a correction of Mardia's, (1970) multivariate skewness
#' and kurtosis test for multivariate samples.
#'
#' @export
#'
#' @author Asael Alonzo Matamoros.
#'
#' @seealso \code{\link{lobato.test}}
#'
#' @references
#' El Bouch, S., Michel, O. & Comon, P.  (2022). A normality test for Multivariate dependent
#' samples. \emph{Journal of Signal Processing}. Volume 201.
#'
#' Mardia, K. (1970). Measures of multivariate skewness and kurtosis with applications.
#' \emph{Biometrika}, 57 519-530
#'
#' Lobato, I., & Velasco, C. (2004). A simple test of normality in time series.
#' \emph{Journal of econometric theory}. 20(4), 671-689.
#'
#' @examples
#' # Generate an univariate stationary arma process
#' y = arima.sim(100,model = list(ar = 0.3))
#' elbouch.test(y)
#'
#' # Generate a bivariate Gaussian random vector
#' x = rnorm(200)
#' y = rnorm(200)
#' elbouch.test(y = y, x = x)
#'
elbouch.test <- function(y, x = NULL){
  if( !is.numeric(y) & !is(y, class2 = "ts") )
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

  if(!is.null(x)){
    if( !is.numeric(x) & !is(x, class2 = "ts") )
      stop("x object must be numeric or a time series")

    if( anyNA(x) )
      stop("The time series contains missing values")

    # checking seasonality
    if(frequency(x) > 1){
      cc = seasonal.test(x)
      if(cc$seasonal)
        warning("x has a seasonal unit root, lobato.test requires stationary process")
    }

    # checking stationarity
    cc = uroot.test(x)
    if(!cc$stationary)
      warning("x has a unit root, lobato.test requires stationary process")
  }

  if(is.null(x)){
    dname = deparse(substitute(y))
    tstat = elbouch.statistic(y - mean(y))[4]
  }else{
    dname = "w = (y, x)"
    tstat = elbouch.statistic(y - mean(y), x - mean(x))[4]
  }

  alt = paste(dname,"does not follow a Gaussian Process")
  names(tstat) <- "Z"

  pval = pnorm(q = tstat, lower.tail = FALSE)

  rval <- list(statistic = tstat,
               p.value = pval,
               alternative = alt,
               method = "El Bouch, Michel & Comon's test",
               data.name = dname)

  class(rval) <- "htest"
  return(rval)
}

#' Computes El Bouch, et al.'s z statistic.
#'
#' Computes the El Bouch, Michel, & Comon's z test statistic for normality of
#' a univariate or bivariate time series.
#'
#' @usage elbouch.statistic(y, x = NULL)
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a stationary time series.
#' @param x a numeric vector or an object of the \code{ts} class containing a stationary time series.
#'
#' @details This function computes Mardia's standardized `z = (B - E_B)/ sd_B` statistic
#' corrected by El Bouch, et al. (2022) for stationary bivariate time series. Where:
#' `B` is the square of a quadratic form of the process `c(y, x)`; `E_B` and `sd_B` are
#' the estimator's expected value and standard error respectively.  If `x` is set to
#' `NULL`, the test computes the univariate counterpart.
#'
#' @return a real value with El Bouch test's statistic.
#'
#' @export
#'
#' @author Asael Alonzo Matamoros.
#'
#' @seealso \code{\link{lobato.statistic}}
#'
#' @references
#' El Bouch, S., Michel, O. & Comon, P.  (2022). A normality test for Multivariate dependent
#' samples. \emph{Journal of Signal Processing}. Volume 201.
#'
#' Mardia, K. (1970). Measures of multivariate skewness and kurtosis with applications.
#' \emph{Biometrika}, 57 519-530
#'
#' Lobato, I., & Velasco, C. (2004). A simple test of normality in time series.
#' \emph{Journal of econometric theory}. 20(4), 671-689.
#'
#' @examples
#' # Genere an univariate stationary ARMA process
#' y = arima.sim(100,model = list(ar = 0.3))
#' elbouch.statistic(y)
#'
#' # Generate a bivariate Gaussian random vector
#' x = rnorm(200)
#' y = rnorm(200)
#' elbouch.statistic(y = y, x = x)
#'
elbouch.statistic <- function(y, x = NULL){

  if(!is.null(x)){

    n = length(y) - 1
    S = cov(cbind(y,x))
    S1 = solve(S)

    B = sum( apply(cbind(y,x), 1, function(x) t(x) %*% S1 %*% x )^2 ) / (n + 1)

    s12 = ccf(x = y, y = x, lag.max = n, type = "covariance", plot = FALSE)$acf[1:n, , ]
    s21 = ccf(x = x, y = y, lag.max = n, type = "covariance", plot = FALSE)$acf[1:n, , ]
    s11 = acf(y, lag.max = n, type = "covariance", plot = FALSE)$acf[1:n, , ]
    s22 = acf(x, lag.max = n, type = "covariance", plot = FALSE)$acf[1:n, , ]

    s = matrix(c(s11, s22, s12, s21), ncol = 4)
    Q1 = apply(s, 1, compute_q1_temp, S = S) * seq(1, n) / ( det(S)^2 )
    Q2 = apply(s, 1, compute_q2_temp, S = S) * seq(1, n) / ( det(S)^4 )

    mu = 8 - 16/(n+1) - 4  * sum(Q1) / (n+1)^2
    sigma2 = 64/(n+1) + 16 * sum(Q2) / (n+1)^2

    z = ( B - mu ) / sqrt(sigma2)
  }
  else{
    n  = length(y)
    S  = var(y)
    S1 = diag(2) / S
    B  = sum( apply(cbind(y,rep(0, n)), 1, function(x) t(x) %*% S1 %*% x )^2 ) / (n)

    s1 = acf(y, lag.max = n, type = "covariance", plot = FALSE)$acf[2:n, , ]

    Q1 = 12 * sum(s1^2  * seq(n-1, 1) / S  ) / (n ^ 2 )
    Q2 =  2 * sum(s1^4  * seq(n-1, 1) / S ^ 2) / (n ^ 2 )

    mu     = 3 - 6 / n - Q1
    sigma2 = 24/ n  + 24* Q2

    z = ( B - mu ) / sqrt(sigma2)
  }

  return(c(mu, sigma2, pnorm(z,lower.tail = FALSE), z))
}

#'
#' @noRd
#'
compute_q1_temp <- function(s, S){

  p1 = S[1,1] * S[2,2] * ( (s[3] + s[4] ) ^ 2 - 4 * s[1] * s[2] )
  p2 = S[1,2] * S[1,2] * ( 2 * (s[3] + s[4] ) ^ 2 + 4 * s[1] * s[2] )
  p3 = S[2,2] * S[1,2] * s[1] * ( s[3] + s[4] )
  p4 = S[1,1] * S[1,2] * s[2] * ( s[3] + s[4] )
  p5 = ( S[1,1] * s[2] ) ^ 2 + (S[2,2] * s[1] ) ^ 2

  p1 + p2 + 6 * (p5 - p4 - p3)
}

#'
#' @noRd
#'
compute_q2_temp <- function(s, S){
  s1s2  = s[1] * s[2]
  s12   = s[3] * s[4]
  s1_s2 = s[3] + s[4]

  p1 = (S[1,1]*S[2,2])^2 * ( 2*(s1s2^2) - 16*(s1s2*s12) + 3*(s[3]^2 + s[4]^2)^2
                             + 12*(s1s2*s1_s2^2) -(2*s12)^2 )

  p2 = (S[1,1]*S[1,2])^2 * ( 8*s1s2 + 3*(5*s1s2 + s12) * (s1_s2^2)
                             -4*s12*(s[2]^2 - s12) )

  p3 = (S[2,2]*S[1,2])^2 * ( 8*s1s2 + 3*(5*s1s2 + s12) * (s1_s2^2)
                             -4*s12*(s[1]^2 - s12) )

  p4 = (S[1,1]*s[2])^4 + (S[2,2]*s[1])^4

  p5 = (S[1,2]^4) * ( (s1s2)^2 + 4*s1s2*s12 + (s12^2) )

  p6 = s1_s2 * ( S[1,1]*S[2,2]*(2*s1s2 + s[3]^2 + s[4]^2) +2*(S[1,2]^2)*(s1s2 + s12) )

  p7 = S[1,1]*S[1,2]*s[2] + S[2,2]*S[1,2]*s[1]

  p1 + (2 * (p2 + p3)) + (3 * p4) + (8 * p5) - (12 * p7 * p6)
}
