#' The k random projections test for normality.
#'
#' Performs the random projection test for normality. The null hypothesis (H0)
#' is that the given data follows a stationary Gaussian process.
#'
#' @usage rp.test(y, k = 1, FDR = TRUE, pars1 = c(100,1), pars2  = c(2,7),
#'                seed = NULL)
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a
#' stationary time series.
#' @param k an integer that determines the `2k` random projections are used for
#' every type of test. The `pars1` argument generates the first `k` projections,
#' and `pars2` generates the later `k` projections. By default, \code{k = 1}.
#' @param FDR a logical value for mixing the p.values using a False discovery
#' rate method. If \code{FDR = TRUE}, then the p.values are mixed using Benjamin
#' and Yekutieli (2001) False discovery Rate method for dependent procedures, on
#' the contrary, it applies Benjamini and Hochberg (1995) procedure.
#' By default, \code{FDR = TRUE}.
#' @param pars1 an optional real vector with the shape parameters of the beta
#' distribution used for the first `k` random projections By default,
#' \code{pars1 = c(100,1)} where, \code{shape1 = 100} and \code{shape2 = 1}.
#' @param pars2 an optional real vector with the shape parameters of the beta
#' distribution used to compute the last `k` random projections. By default,
#' \code{pars2 = c(2,7)} where, \code{shape1 = 2} and \code{shape2 = 7}.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#'
#' @return A list with class \code{"h.test"} containing the following components:
#'  \item{statistic:}{an integer value with the amount of projections per test.}
#'  \item{parameter:}{a text that specifies the p.value mixing FDR method.}
#'  \item{p.value:}{the FDR mixed p-value for the test.}
#'  \item{alternative:}{a character string describing the alternative hypothesis.}
#'  \item{method:}{a character string \dQuote{k random projections test}.}
#'  \item{data.name:}{a character string giving the name of the data.}
#'
#' @details
#' The random projection test generates `2k` random projections of `y`. Applies
#' Epps statistics to the odd projections and Lobato and Velasco’s statistics to
#' the even ones. Computes the `2k` p.values using an asymptotic chi-square distribution
#' with two degrees of freedom. Finally, mixes the p.values using a false discover
#' rate procedure. By default, mixes the p.values using Benjamin and Yekutieli’s
#' (2001) method.
#'
#' The function uses beta distributions for generating the `2k` random projections.
#' By default, uses a \code{beta(shape1 = 100,shape = 1)} distribution contained
#' in \code{pars1} argument to generate the first `k` projections. For the later
#' `k` projections the functions uses a \code{beta(shape1 = 2,shape = 7)} distribution
#' contained in \code{pars2} argument.
#'
#' The test was proposed by \emph{Nieto-Reyes, A.,Cuesta-Albertos, J. & Gamboa,
#' F. (2014)}.
#'
#' @export
#'
#' @author Asael Alonzo Matamoros and Alicia Nieto-Reyes.
#'
#' @seealso \code{\link{lobato.test}}, \code{\link{epps.test}}
#'
#' @references
#' Nieto-Reyes, A., Cuesta-Albertos, J. & Gamboa, F. (2014). A random-projection
#' based test of Gaussianity for stationary processes. \emph{Computational
#' Statistics & Data Analysis, Elsevier}, vol. 75(C), pages 124-141.
#'
#' Lobato, I., & Velasco, C. (2004). A simple test of normality in time series.
#' \emph{Journal of econometric theory}. 20(4), 671-689.
#'
#' Benjamini, Y., and Yekutieli, D. (2001). The control of the false discovery rate in
#' multiple testing under dependency. \emph{Annals of Statistics}. 29, 1165–1188.
#' Doi:10.1214/aos/1013699998.
#'
#' Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of significance.
#' \emph{Biometrika}. 75, 800–803. Doi:10.2307/2336325.
#'
#' Epps, T.W. (1987). Testing that a stationary time series is Gaussian. \emph{The
#' Annals of Statistic}. 15(4), 1683-1698.
#'
#' @examples
#' # Generating an stationary arma process
#' y = arima.sim(100,model = list(ar = 0.3))
#' rp.test(y)
#'
rp.test = function(y, k = 1, FDR = TRUE, pars1 = c(100,1), pars2 = c(2,7), seed = NULL){

  if( !is.numeric(y) & !is(y,class2 = "ts") )
    stop("y object must be numeric or a time series")

  if( anyNA(y) )
    stop("The time series contains missing values")

  # checking stationarity
  cc = uroot.test(y)
  if(!cc$stationary)
    warning("y has a unit root, epps.test requires stationary process")

  # checking seasonality
  if(frequency(y) > 1){
    cc = seasonal.test(y)
    if(cc$seasonal)
      warning("y has a seasonal unit root, epps.test requires stationary process")
  }

  if (!is.null(seed))
    set.seed(seed)

  rps = rp.sample(as.numeric(y), k = k, seed = seed, pars1 = pars1, pars2 = pars2)
  F1 = pchisq(q = c(rps$lobato, rps$epps), df = 2, lower.tail = FALSE)

  if(FDR)
    F1 = min(p.adjust(F1, method = "BY"))
  else
    F1 = min(p.adjust(F1, method = "BH"))

  dname = deparse(substitute(y))
  alt = paste(dname,"does not follow a Gaussian Process")
  stat = k
  names(stat) = "k"

  # tests parameters
  parameters  = ifelse(FDR, "Benjamini & Yekutieli", "Benjamini & Hocheberg")
  names(parameters) = "p.value adjust"
  names(F1) = "fdr.value"

  #htest class
  rval <- list(statistic = stat,
               parameter = parameters,
               p.value = F1,
               alternative = alt,
               method = "k random projections test.",
               data.name = dname)

  class(rval) <- "htest"

  return(rval)
}
#' Generates a test statistics sample of random projections.
#'
#' Generates a 2k sample of test statistics  projecting the stationary process
#' using the random projections procedure.
#'
#' @usage rp.sample(y, k = 1, pars1 = c(100,1), pars2 = c(2,7), seed = NULL)
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a
#' stationary time series.
#' @param k an integer that determines the `2k` random projections are used for
#' every type of test. The `pars1` argument generates the first `k` projections,
#' and `pars2` generates the later `k` projections. By default, \code{k = 1}.
#' @param pars1 an optional real vector with the shape parameters of the beta
#' distribution used for the first `k` random projections By default,
#' \code{pars1 = c(100,1)} where, \code{shape1 = 100} and \code{shape2 = 1}.
#' @param pars2 an optional real vector with the shape parameters of the beta
#' distribution used to compute the last `k` random projections. By default,
#' \code{pars2 = c(2,7)} where, \code{shape1 = 2} and \code{shape2 = 7}.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#'
#' @return A list with 2 real value vectors:
#'  \item{lobato:}{A vector with the Lobato and Velasco's statistics sample.}
#'  \item{epps:}{A vector with the Epps statistics sample.}
#'
#' @details
#' The \code{rp.sample} function generates `2k` tests statistics by projecting
#' the time series using `2k` stick breaking processes. First, the function
#' samples `k` stick breaking processes using \code{pars1} argument. Then, projects
#' the time series using the sampled stick processes. Later, applies the Epps
#' statistics to the odd projections and the Lobato and Velasco’s statistics to
#' the even ones. Analogously, the function performs the three steps using also
#' \code{pars2} argument
#'
#' The function uses beta distributions for generating the `2k` random projections.
#' By default, uses a \code{beta(shape1 = 100,shape = 1)} distribution contained
#' in \code{pars1} argument to generate the first `k` projections. For the later
#' `k` projections the functions uses a \code{beta(shape1 = 2,shape = 7)} distribution
#' contained in \code{pars2} argument.
#'
#' The test was proposed by \emph{Nieto-Reyes, A.,Cuesta-Albertos, J. & Gamboa, F. (2014)}.
#'
#' @export
#'
#' @author Alicia Nieto-Reyes and Asael Alonzo Matamoros
#'
#' @seealso \code{\link{lobato.test}}, \code{\link{epps.test}}
#'
#' @references
#' Nieto-Reyes, A., Cuesta-Albertos, J. & Gamboa, F. (2014). A random-projection
#' based test of Gaussianity for stationary processes. \emph{Computational
#' Statistics & Data Analysis, Elsevier}, vol. 75(C), pages 124-141.
#'
#' Epps, T.W. (1987). Testing that a stationary time series is Gaussian. \emph{The
#' Annals of Statistic}. 15(4), 1683-1698.
#'
#' Lobato, I., & Velasco, C. (2004). A simple test of normality in time series.
#' \emph{Journal of econometric theory}. 20(4), 671-689.
#'
#' @examples
#' # Generating an stationary ARMA process
#' y = arima.sim(100,model = list(ar = 0.3))
#' rp.sample(y)
#'
rp.sample = function(y, k = 1, pars1 = c(100,1), pars2 = c(2,7), seed = NULL){

  if( !is.numeric(y) & !is(y,class2 = "ts") )
    stop("y object must be numeric or a time series")

  if( anyNA(y) )
    stop("The time series contains missing values")

  if (!is.null(seed))
    set.seed(seed)

  if(k <= 1)
    k = 1

  k2 = ifelse(k %% 2 == 0, k/2,(k+1)/2)

  x1 = parallel::mclapply(1:k2, FUN = function(i){
    yh1 = random.projection(as.numeric(y),shape1 = pars1[1],shape2 = pars1[2])
    yh2 = random.projection(as.numeric(y),shape1 = pars1[1],shape2 = pars1[2])
    c(i,epps.statistic(yh1),lobato.statistic(yh2))
  })

  x2 = parallel::mclapply(1:k2, FUN = function(i){
    yh1 = random.projection(as.numeric(y),shape1 = pars2[1],shape2 = pars2[2])
    yh2 = random.projection(as.numeric(y),shape1 = pars2[1],shape2 = pars2[2])
    c(i,epps.statistic(yh1),lobato.statistic(yh2))
  })

  x = rbind(matrix(unlist(x1),ncol = 3,byrow = TRUE),
            matrix(unlist(x2),ncol = 3,byrow = TRUE))

  rp.sample = list(lobato = x[1:k, 2],epps = x[1:k ,3])

  return(rp.sample)
}
#' Generate a random projection.
#'
#' Generates a random projection of a univariate stationary stochastic process. Using
#' a beta(shape1,shape2) distribution.
#'
#' @usage random.projection(y,shape1,shape2,seed = NULL)
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a stationary time series.
#' @param shape1 an optional real value with the first shape parameters of the beta
#' distribution.
#' @param shape2 an optional real value with the second shape parameters of the beta
#' distribution.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#'
#' @return a real vector with the projected stochastic process.
#'
#' @details
#' Generates one random projection of a stochastic process using a beta distribution.
#' For more details, see: \emph{Nieto-Reyes, A.,Cuesta-Albertos, J. & Gamboa, F. (2014)}.
#'
#' @export
#'
#' @author Alicia Nieto-Reyes and Asael Alonzo Matamoros.
#'
#' @seealso \code{\link{lobato.test}} \code{\link{epps.test}}
#'
#' @references
#' Nieto-Reyes, A., Cuesta-Albertos, J. & Gamboa, F. (2014). A random-projection
#' based test of Gaussianity for stationary processes. \emph{Computational
#' Statistics & Data Analysis, Elsevier}, vol. 75(C), pages 124-141.Result
#'
#' Epps, T.W. (1987). Testing that a stationary time series is Gaussian. \emph{The
#' Annals of Statistic}. 15(4), 1683-1698.
#'
#' Lobato, I., & Velasco, C. (2004). A simple test of normality in time series.
#' \emph{Journal of econometric theory}. 20(4), 671-689.
#'
#' @examples
#' # Generating an stationary arma process
#' y = arima.sim(100,model = list(ar = 0.3))
#' rp.test(y)
#'
random.projection = function(y,shape1,shape2,seed = NULL){

  if( !is.numeric(y) & !is(y,class2 = "ts") )
    stop("y object must be numeric or a time series")

  if( anyNA(y) )
    stop("The time series contains missing values")

  if (!is.null(seed))
    set.seed(seed)

  n = length(y);yp = rep(0,n);ch = 1; C = n
  H = rbeta(n,shape1,shape2)

  while(ch > 10^(-15) && C > 1){

    a = ch*H[C]; ch=ch-a;
    if (C == n)
      H[C] = sqrt(a)
    else
      H[C] = sqrt(a)/(n-C)

    C = C-1
  }
  H[C] = sqrt(ch)/(n-C); h = H[C:n]; k = length(h);

  yp[1:(k-1)] = sapply(1:(k-1), function(j) sum(y[1:j]*h[(k-j+1):k]))

  yp[k:n] = sapply(k:n, function(j) sum(y[(j-k+1):j]*h))

  return(yp)
}
