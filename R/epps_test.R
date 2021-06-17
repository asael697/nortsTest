#' The Epps and Pulley Test for normality.
#'
#' Performs the Epps test for normality. The null hypothesis (H0) is that the given data
#' follows a stationary Gaussian process.
#'
#' @usage  epps.test(y)
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a stationary time series.
#'
#' @return a h.test class with the main results of the Epps hypothesis test. The
#' h.test class have the following values:
#' \itemize{
#'  \item{"epps"}{The Epps statistic}
#'  \item{"df"}{The test degrees freedoms}
#'  \item{"p.value"}{The p value}
#'  \item{"alternative"}{The alternative hypothesis}
#'  \item{"method"}{The used method}
#'  \item{"data.name"}{The data name.}
#' }
#'
#' @details
#' The Epps test minimize the process' empirical characteristic function using a
#' quadratic loss in terms of the process two first moments. The test was proposed
#' by \emph{Epps, T.W. (1987)} and implemented by \emph{Nieto-Reyes, A.,
#' Cuesta-Albertos, J. & Gamboa, F. (2014)} using the \code{amoebam()} function of
#' \emph{Press, W.H., Teukolsky, S.A., Vetterling, W.T. and  Flannery, B.P. (2007)}.
#'
#' @export
#'
#' @author Asael Alonzo Matamoros and Alicia Nieto-Reyes.
#'
#' @seealso \code{\link{lobato.test}}
#'
#' @references
#' Epps, T.W. (1987). Testing that a stationary time series is Gaussian. \emph{The
#' Annals of Statistic}. 15(4), 1683-1698.
#'
#' Nieto-Reyes, A., Cuesta-Albertos, J. & Gamboa, F. (2014). A random-projection
#' based test of Gaussianity for stationary processes. \emph{Computational
#' Statistics & Data Analysis, Elsevier}, vol. 75(C), pages 124-141.
#'
#' Press, W.H., Teukolsky, S.A., Vetterling, W.T. and  Flannery, B.P. (2007).
#' Numerical Recipes. The Art of Scientific Computing. \emph{Cambridge
#' University Press}.
#'
#' @examples
#' # Generating an stationary arma process
#' y = arima.sim(100,model = list(ar = 0.3))
#' epps.test(y)
#'
epps.test = function(y){

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

  dname = deparse(substitute(y))
  alt = paste(dname,"does not follow a Gaussian Process")
  tstat = epps.statistic(y)
  names(tstat) <- "epps"
  df = 2
  names(df) <- "df"
  pval = pchisq(q = tstat,df = 2,lower.tail = FALSE)

  rval <- list(statistic =tstat, parameter = df, p.value = pval,
                 alternative = alt,
                 method = "Epps test", data.name = dname)
  class(rval) <- "htest"
  return(rval)
}
#' Estimates the Epps statistic
#'
#' Estimates the Epps statistic minimizing the quadratic loss of the process'
#' characteristic function in terms of the first two moments.
#'
#' @usage  epps.statistic(y)
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a stationary time series.
#'
#' @details This function is the equivalent of \code{Sub} in \emph{Nieto-Reyes, A.,
#' Cuesta-Albertos, J. & Gamboa, F. (2014)}. This function uses a quadratic
#' optimization solver implemented by \emph{Press, W.H., Teukolsky, S.A.,
#' Vetterling, W.T. and  Flannery, B.P. (2007)}.
#'
#' @return a real value with the Epps test's statistic.
#'
#' @importFrom MASS ginv
#' @export
#'
#' @author Alicia Nieto-Reyes and Asael Alonzo Matamoros.
#'
#' @seealso \code{\link{lobato.statistic}}
#'
#' @references
#' Epps, T.W. (1987). Testing that a stationary time series is Gaussian. \emph{The
#' Annals of Statistic}. 15(4), 1683-1698.
#'
#' Nieto-Reyes, A., Cuesta-Albertos, J. & Gamboa, F. (2014). A random-projection
#' based test of Gaussianity for stationary processes. \emph{Computational
#' Statistics & Data Analysis, Elsevier}, vol. 75(C), pages 124-141.
#'
#' Press, W.H., Teukolsky, S.A., Vetterling, W.T. and  Flannery, B.P. (2007).
#' Numerical Recipes. The Art of Scientific Computing. \emph{Cambridge
#' University Press}.
#'
#' @examples
#' # Generating an stationary arma process
#' y = arima.sim(100,model = list(ar = 0.3))
#' epps.statistic(y)
#'
epps.statistic =  function(y){
  n = length(y);dN=4;N = 2
  deviSt = sd(y)*(n-1)/n
  lambda = c(1,2)/deviSt
  rn = floor(n^.4)
  gmatrix = matrix(0,n,dN)
  for (i in 1:n) {
    ly = lambda*y[i];
    co = cos(ly);
    si = sin(ly);

    for (j in 1:N) gmatrix[i,((j*2)-1):(j*2)] = c(co[j], si[j]);
  }
  gn = apply(gmatrix, 2, mean)
  dpifcero = matrix(0,dN,dN)
  de2 = matrix(0,dN,dN)

  for (j in 1:n) {
    zj = gmatrix[j,]-gn;
    dpifcero = dpifcero+(cbind(zj) %*% zj)
  }
  for (r in 1:rn) {
    de1 = matrix(0,dN,dN)
    for (j in 1:(n-r)) de1 = de1+(cbind(gmatrix[j,]-gn) %*% (gmatrix[j+r,]-gn))
    de2 = de2+de1*(1-r/rn)
  }

  dpifcero = (dpifcero+2*de2)/n;
  Gm = MASS::ginv(dpifcero);
  me = mean(y); sts = deviSt/sqrt(n); ts2 = sqrt(2/n); Va = deviSt^2;
  P  = matrix(c(me-sts,me+sts,me,Va*(1-ts2),Va*(1-ts2),Va*(1+ts2)),nrow=3,ncol=2)
  Y  = cbind(rep(0,3))

  for (i in 1:3) Y[i] = Quadratic(P[i,],gn,lambda,Gm,N,dN)

  return(amoebam(P,Y,n,gn,lambda,Gm,N,dN))
}
#'
#' @noRd
#'
Quadratic = function(m, gn, lambda, Gm, N, dN){
  mu = m[1]; sigma = m[2];
  ml = mu*lambda
  e = exp(-sigma*(lambda^2)/2)
  re = e*cos(ml); im = e*sin(ml); gms = rep(0,dN);
  for (j in 1:N) { j2 = j*2; gms[(j2-1):j2] = c(re[j], im[j]); }
  g = gn-gms
  q = (g %*% Gm %*% cbind(g))
  return(q[1,1])
}
#'
#' @noRd
#'
amoebam = function(P,Y,n,gn,lambda,Gm,N,dN){
  NDIM = 2; FTOL = .0001; NMAX = 20; ALPHA = 1;
  BETA = 0.5; GAMMA = 2; ITMAX = 500; PR = rep(0,NMAX);
  PRR = rep(0,NMAX); MPTS = NDIM+1; ITER = 0; IXXX = 1;
  while (IXXX==1) {
    ILO = 1; if(Y[1]>Y[2]){ IHI = 1; INHI = 2;} else { IHI = 2; INHI = 1; }
    for (I in 1:MPTS) {
      if(Y[I]<Y[ILO]){ ILO = I }
      if(Y[I]>Y[IHI]){ INHI = IHI; IHI = I; }
      else if(Y[I]>Y[INHI] && I != IHI){ INHI = I }
    }
    RTOL = 2*abs(Y[IHI]-Y[ILO])/(abs(Y[IHI])+abs(Y[ILO]))
    if(!is.na(RTOL) && RTOL<FTOL) { return(n*min(Y)) }
    if(ITER==ITMAX) { return(n*min(Y)) }
    ITER = ITER+1; PBAR = rep(0,NDIM);
    for (I in 1:MPTS) {
      if(I != IHI){ for (J in 1:NDIM) { PBAR[J] = PBAR[J]+P[I,J] } }
    }
    for (J in 1:NDIM) {
      PBAR[J] = PBAR[J]/NDIM
      PR[J] = (1.+ALPHA)*PBAR[J]-ALPHA*P[IHI,J]
    }
    YPR = Quadratic(PR,gn,lambda,Gm,N,dN)
    if (YPR<=Y[ILO]){
      for (J in 1:NDIM) { PRR[J] = GAMMA*PR[J]+(1.-GAMMA)*PBAR[J] }
      YPRR = Quadratic(PRR,gn,lambda,Gm,N,dN)
      if (YPRR < Y[ILO]){
        for (J in 1:NDIM) { P[IHI,J] = PRR[J] }
        Y[IHI] = YPRR
      }
      else {
        for (J in 1:NDIM) { P[IHI,J] = PR[J] }; Y[IHI] = YPR;
      }
    }
    else if(YPR >= Y[INHI]){
      if(YPR < Y[IHI]){
        for (J in 1:NDIM) { P[IHI,J] = PR[J] }
        Y[IHI] = YPR
      }
      for (J in 1:NDIM) { PRR[J] = BETA*P[IHI,J]+(1.-BETA)*PBAR[J] }
      YPRR = Quadratic(PRR,gn,lambda,Gm,N,dN)
      if(YPRR<Y[IHI]){
        for (J in 1:NDIM) { P[IHI,J] = PRR[J] }
        Y[IHI] = YPRR
      }
      else{
        for (I in 1:MPTS) {
          if(I != ILO){
            for (J in 1:NDIM) {
              PR[J] = 0.5*(P[I,J]+P[ILO,J]); P[I,J] = PR[J];
            }
            Y[I] = Quadratic(PR,gn,lambda,Gm,N,dN)
          }
        }
      }
    } else{
      for (J in 1:NDIM) { P[IHI,J] = PR[J] }; Y[IHI] = YPR;
    }
  }
}
