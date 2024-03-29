% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epps_test.R
\name{epps_bootstrap.test}
\alias{epps_bootstrap.test}
\title{The Sieve Bootstrap Epps and Pulley test for normality.}
\usage{
epps_bootstrap.test(y, lambda = c(1,2), reps = 500, h = 100, seed = NULL)
}
\arguments{
\item{y}{a numeric vector or an object of the \code{ts} class containing a stationary
time series.}

\item{lambda}{a numeric vector for evaluating the characteristic function.}

\item{reps}{an integer with the total bootstrap repetitions.}

\item{h}{an integer with the first \code{burn-in} sieve bootstrap replicates.}

\item{seed}{An optional \code{\link[=set.seed]{seed}} to use.}
}
\value{
A list with class \code{"h.test"} containing the following components:
 \item{statistic:}{the sieve bootstrap Epps and Pulley's statistic.}
 \item{p.value:}{the p value for the test.}
 \item{alternative:}{a character string describing the alternative hypothesis.}
 \item{method:}{a character string \dQuote{Sieve-Bootstrap Epps' test}.}
 \item{data.name:}{a character string giving the name of the data.}
}
\description{
Performs the approximated Epps and Pulley's test of normality for univariate time series.
Computes the p-value using Psaradakis and Vavra's (2020) sieve bootstrap procedure.
}
\details{
The Epps test minimize the process' empirical characteristic function using a
quadratic loss in terms of the process two first moments, \emph{Epps, T.W. (1987)}.
Approximates the p-value using a sieve-bootstrap procedure \emph{Psaradakis, Z.
and Vávra, M. (2020)}.
}
\examples{
# Generating an stationary arma process
y = arima.sim(300, model = list(ar = 0.3))
epps_bootstrap.test(y, reps = 1000)

}
\references{
Psaradakis, Z. and Vávra, M. (2020) Normality tests for dependent
data: large-sample and bootstrap approaches. Communications in
\emph{Statistics-Simulation and Computation 49 (2)}. ISSN 0361-0918.

Nieto-Reyes, A., Cuesta-Albertos, J. & Gamboa, F. (2014). A random-projection
based test of Gaussianity for stationary processes. \emph{Computational
Statistics & Data Analysis, Elsevier}, vol. 75(C), pages 124-141.

Epps, T.W. (1987). Testing that a stationary time series is Gaussian. \emph{The
Annals of Statistic}. 15(4), 1683-1698.
}
\seealso{
\code{\link{lobato.statistic}}, \code{\link{epps.test}}
}
\author{
Asael Alonzo Matamoros and Alicia Nieto-Reyes.
}
