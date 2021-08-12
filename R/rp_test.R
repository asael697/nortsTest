#' The k random projections test for normality
#'
#' Performs the random projection test for normality. The null hypothesis (H0) is that the given data
#' follows a stationary Gaussian process, and k is the number of used random projections.
#'
#' @usage  rp.test(y,k = 16,FDR = TRUE,pars1 = c(100,1),pars2  = c(2,7),seed = NULL)
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a stationary time series.
#' @param k an integer with the number of random projections to be used, by default
#' \code{k = 16}.
#' @param FDR a logical value for mixing the p-values using a dependent False discovery
#' rate method. If \code{FDR =TRUE}, then the p-values are mixed using a False discovery Rate method,
#' on the contrary it applies the Benjamin and Yekuteli (2001) procedure.  By default \code{FDR = TRUE}.
#' @param pars1 an optional real vector with the shape parameters of the beta distribution
#' used for the odd number random projection. By default, \code{pars1 = c(100,1)} where,
#' \code{shape1 = 100} and \code{shape2 = 1}.
#' @param pars2 an optional real vector with the shape parameters of the beta distribution
#' used for the even number random projection. By default, \code{pars2 = c(2,7)} where,
#' \code{shape1 = 2} and \code{shape2 = 7}.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#'
#' @return a h.test class with the main results of the Epps hypothesis test. The
#' h.test class have the following values:
#' \itemize{
#'  \item{"k"}{The number of used projections}
#'  \item{"lobato"}{The average Lobato and Velasco's test statistics of the k projected samples}
#'  \item{"epps"}{The average Epps test statistics of the k projected samples}
#'  \item{"p.value"}{The mixed p value}
#'  \item{"alternative"}{The alternative hypothesis}
#'  \item{"method"}{The used method: rp.test}
#'  \item{"data.name"}{The data name.}
#' }
#'
#' @details
#' The random projection test generates k independent random projections of the process.
#' A Lobato and Velasco's test are applied to the first half of the projections, and an
#' Epps test for the other half. Then, a False discovery rate is used for mixing the
#' obtained p.values
#'
#' For generating the k random projections a beta distribution is used. By default a
#' \code{beta(shape1 = 100,shape = 1)} and a \code{beta(shape1 = 2,shape = 7)} are used
#' to generate the odd and even projections respectively. For using a different parameter
#' set, change \code{pars1} or \code{pars2}.
#'
#' The test was proposed by \emph{Nieto-Reyes, A.,Cuesta-Albertos, J. & Gamboa, F. (2014)}.
#'
#' @export
#'
#' @author Asael Alonzo Matamoros and Alicia Nieto-Reyes.
#'
#' @seealso \code{\link{lobato.test}} \code{\link{epps.test}}
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
#' # Generating an stationary arma process
#' y = arima.sim(100,model = list(ar = 0.3))
#' rp.test(y,k = 4)
#'
rp.test = function(y,k = 16,FDR = TRUE,pars1 = c(100,1),pars2 = c(2,7),seed = NULL){

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

  rps = rp.sample(as.numeric(y),k = k,seed = seed,pars1 = pars1,pars2 = pars2)
  F1 = pchisq(q = c(rps$lobato,rps$epps),df = 2,lower.tail = FALSE)

  if(FDR)
    F1 = min(p.adjust(F1,method = "fdr"))
  else
    F1 = min(p.adjust(F1,method = "BY"))

  dname = deparse(substitute(y))
  alt = paste(dname,"does not follow a Gaussian Process")
  stat = k
  names(stat) = "k"

  # tests parameters
  lobato.stat = ifelse(is.null(rps$lobato),0,mean(rps$lobato))
  epps.stat   = ifelse(is.null(rps$epps),0,mean(rps$epps))
  parameters  = c(lobato.stat,epps.stat)
  names(parameters) = c("lobato","epps")

  #htest class
  rval <- list(statistic = stat,parameters = parameters, p.value = F1,
               alternative = alt,
               method = "k random projections test", data.name = dname)
  class(rval) <- "htest"

  return(rval)
}
#' Generates a test statistics sample of random projections
#'
#' Generates a sample of test statistics using k independent random projections
#' of a stationary process. The first half values of the sample, are estimated
#' using a Lobato and Velasco's statistic test. The last half values with an Epps
#' statistic test.
#'
#' @usage  rp.sample(y,k = 16,pars1 = c(100,1),pars2 = c(2,7),seed = NULL)
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a stationary time series.
#' @param k an integer with the number of random projections to be used, by default \code{k = 16}.
#' @param pars1 an optional real vector with the shape parameters of the beta distribution
#' used for the odd number random projection. By default, \code{pars1 = c(100,1)} where,
#' \code{shape1 = 100} and \code{shape2 = 1}.
#' @param pars2 an optional real vector with the shape parameters of the beta distribution
#' used for the even number random projection. By default, \code{pars2 = c(2,7)} where,
#'\code{shape1 = 2} and \code{shape2 = 7}.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#'
#' @return A list with 2 real value vectors:
#' \itemize{
#'  \item{"lobato"}{A vector with the Lobato and Velasco's statistics sample}
#'  \item{"epps"}{A vector with the Epps statistics sample.}
#' }
#'
#' @details
#' The \code{rp.sample} function generates k independent random projections of the process.
#' A Lobatos and Velasco's test is applied to the first half of the projections. And an
#' Epps test for the other half.
#'
#' For generating the k random projections a beta distribution is used. By default a
#' \code{beta(shape1 = 100,shape = 1)} and a \code{beta(shape1 = 2,shape = 7)} are used
#' to generate the odd and even projections respectively. For using a different parameter
#' set, change \code{pars1} or \code{pars2} values.
#'
#' The test was proposed by \emph{Nieto-Reyes, A.,Cuesta-Albertos, J. & Gamboa, F. (2014)}.
#'
#' @export
#'
#' @author Alicia Nieto-Reyes and Asael Alonzo Matamoros
#'
#' @seealso \code{\link{lobato.test}} \code{\link{epps.test}}
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
#' # Generating an stationary arma process
#' y = arima.sim(100,model = list(ar = 0.3))
#' rp.test(y,k = 4)
#'
rp.sample = function(y,k = 16,pars1 = c(100,1),pars2 = c(2,7),seed = NULL){

  if( !is.numeric(y) & !is(y,class2 = "ts") )
    stop("y object must be numeric or a time series")

  if( anyNA(y) )
    stop("The time series contains missing values")

  if (!is.null(seed))
    set.seed(seed)

  if(k <= 1)
    k = 1

  n = length(y);T1 = NULL;T2 = NULL

  for(j in 1:k){
    yh1 = random.projection(as.numeric(y),shape1 = pars1[1],shape2 = pars1[2])
    yh2 = random.projection(as.numeric(y),shape1 = pars2[1],shape2 = pars2[2])

    T1 = c(T1,lobato.statistic(yh1) )
    T1 = c(T1,lobato.statistic(yh2) )

    T2 = c(T2,epps.statistic(yh1) )
    T2 = c(T2,epps.statistic(yh2) )
  }
  rp.sample = list(lobato = T1,epps = T2)

  return(rp.sample)
}
#' Generate a random projection
#'
#' Generates a random projection of a univariate stationary stochastic process. Using
#' a beta(shape1,shape2) distribution.
#'
#' @usage  random.projection(y,shape1,shape2,seed = NULL)
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
#' @author Alicia Nieto-Reyes and Asael Alonzo Matamoros
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

  for (j in 1:(k-1))
    yp[j] = sum(y[1:j]*h[(k-j+1):k])

  for (j in k:n)
    yp[j] = sum(y[(j-k+1):j]*h)

  return(yp)
}
