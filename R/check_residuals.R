#'
#' @export
#'
check_residuals<- function(y,...) {
  UseMethod("check_residuals")
}
#' Generic functions for checking residuals in time series models
#'
#' Generic function for residuals check analysis, these methods are inspired in the \code{check.residuals} function
#' provided by the \code{forecast} package.
#'
#' @aliases check_residuals check_residuals.ts check_residuals.arima0 check_residuals.Arima
#' check_residuals.fGarch check_residuals.numeric check_residuals.lm  check_residuals.HoltWinters
#' check_residuals.ets check_residuals.forecast
#'
#' @rdname check_residuals
#'
#' @param y 	Either a time series model,the supported classes are \code{arima0}, \code{Arima}, \code{sarima},
#' \code{fGarch}, or a time series (assumed to be residuals).
#' @param normality A character string naming the desired test for checking gaussian distribution.
#' Valid values are \code{"epps"} for the Epps, \code{"lobato"} for Lobato and Velasco's,\code{"vavras"} for
#' the Psaradakis and  Vavra, \code{"rp"} for the random projections, \code{"jb"} for the Jarque and Beras,
#' \code{"ad"} for Anderson Darling test, and \code{"shapiro"} for the Shapiro-Wilk's test. The default value
#' is \code{"epps"} test.
#' @param unit_root A character string naming the desired unit root test for checking stationarity.
#' Valid values are \code{"adf"} for the Augmented Dickey-Fuller, \code{"pp"} for the Phillips-Perron,
#' and \code{"kpss"} for Kwiatkowski, Phillips, Schmidt, and Shin. The default value is \code{"adf"} for the
#' Augmented Dickey-Fuller test.
#' @param seasonal A character string naming the desired unit root test for checking seasonality.
#' Valid values are \code{"ocsb"} for the Osborn, Chui, Smith, and Birchenhall, \code{"ch"} for the
#' Canova and Hansen, and \code{"hegy"} for Hylleberg, Engle, Granger, and Yoo. The default value is
#' \code{"ocsb"} for the Osborn, Chui, Smith, and Birchenhall test.
#' @param arch A character string naming the desired  test for checking stationarity. Valid values are
#' \code{"box"} for the Ljung-Box, and \code{"Lm"} for the Lagrange Multiplier test. The default
#' value is \code{"box"} for the Augmented Ljung-Box test.
#' @param alpha Level of the test, possible values range from 0.01 to 0.1. By default \code{alpha = 0.05}
#' is used
#' @param plot A boolean value. If \code{TRUE}, will produce produces a time plot of the residuals,
#' the corresponding ACF, and a histogram.
#' @param ... Other testing parameters
#'
#' @details The function performs a residuals analysis, it prints a unit root and seasonal test to check
#' stationarity, and a normality test for checking  Gaussian distribution. In addition, if the plot option is
#' \code{TRUE} a time plot, ACF, and histogram of the series are presented.
#'
#' @return The function does not return any value
#'
#' @author Asael Alonzo Matamoros
#'
#' @importFrom forecast  ggtsdisplay
#' @method check_residuals ts
#'
#' @export
#'
#' @references
#' Dickey, D. & Fuller, W. (1979). Distribution of the Estimators for
#' Autoregressive Time Series with a Unit Root. \emph{Journal of the American
#' Statistical Association}. 74, 427-431.
#'
#' Epps, T.W. (1987). Testing that a stationary time series is Gaussian. \emph{The
#' Annals of Statistic}. 15(4), 1683-1698.\url{http://www.jstor.org/stable/2336512}.
#' \code{doi:10.1214/aos/1176350618}
#'
#' Osborn, D., Chui, A., Smith, J., & Birchenhall, C. (1988). Seasonality and the
#' order of integration for consumption. \emph{Oxford Bulletin of Economics
#' and Statistics}. 50(4), 361-377.
#'
#' @examples
#' # Generating an stationary arma process
#' y = arima.sim(100,model = list(ar = 0.3))
#' check_residuals(y,unit_root = "adf")
#'
check_residuals.ts = function(y,normality = "epps",
                              unit_root = NULL,
                              seasonal = NULL,
                              arch = NULL,
                              alpha = 0.05,plot = FALSE,...){

  if( !is(y,class2 = "ts") )
    stop("The object y is not a time series")

  if( anyNA(y) )
    stop("The time series contains missing values")

  # Unit root test
  if(!is.null(unit_root)){
    cc = uroot.test(y = y,unit_root = unit_root,alpha = alpha)
    cat("\n *************************************************** \n")
    cat("\n Unit root test for stationarity: \n")
    print(cc)
    cat("\n",cc$Conc)
    cat("\n *************************************************** \n")
  }

  # seasonal test
  if(frequency(y) > 1){
    if(!is.null(seasonal)){
      cc = seasonal.test(y = y,seasonal = seasonal,alpha = alpha)

      cat("\n Unit root test for seasonality: \n")
      print(cc)
      cat("\n",cc$Conc)
      cat("\n \n *************************************************** \n")
    }
  }

  # arch test
  if(!is.null(arch)){
    cc = arch.test(y = y,arch = arch,alpha = alpha,...)
    cat("\n arch test for heteroscedasticity: \n")
    print(cc)
    cat("\n",cc$Conc)
    cat("\n \n *************************************************** \n")
  }

  # Normality test
  cc = normal.test(y = y,normality =normality,alpha = alpha)
  cat("\n Goodness of fit test for Gaussian Distribution: \n")
  print(cc)
  cat("\n",cc$Conc)
  cat("\n \n *************************************************** \n")

  if (plot) {
   suppressMessages(suppressWarnings(check_plot.ts(y)))
  }
}
#'
#' @method check_residuals numeric
#' @export
#'
check_residuals.numeric = function(y,normality = "epps", unit_root = NULL,seasonal = NULL,
                                   arch = NULL,alpha = 0.05,plot = FALSE,...){
  if( !is.numeric(y) )
    stop("The object y is not a numeric array")

  check_residuals.ts(ts(y,frequency = 1),unit_root = unit_root,
                     normality = normality,
                     seasonal = seasonal,arch = arch,
                     alpha = 0.05,
                     plot = plot,...)
}
#' @method check_residuals arima0
#' @export
#'
check_residuals.arima0 = function(y,normality = "epps",
                                  unit_root = NULL,
                                  seasonal = NULL,
                                  arch = NULL,
                                  alpha = 0.05,plot = FALSE,...){
  if( !is(y,class2 = "arima0") )
    stop("The object y is not an arima0 class")

  check_residuals.ts(y$residuals,unit_root = unit_root,
                     normality = normality,
                     seasonal = seasonal,arch = arch,
                     alpha = 0.05,
                     plot = plot,...)
}
#'
#' @method check_residuals Arima
#' @export
#'
check_residuals.Arima = function(y,normality = "epps",
                                 unit_root = NULL,
                                 seasonal = NULL,
                                 arch = NULL,
                                 alpha = 0.05,plot = FALSE,...){
  if( !is(y,class2 = "Arima") )
    stop("The object y is not a Arima class")

  check_residuals.ts(y$residuals,unit_root = unit_root,
                     normality = normality,
                     seasonal = seasonal,arch = arch,
                     alpha = 0.05,
                     plot = plot,...)
}
#'
#' @method check_residuals fGARCH
#' @export
#'
check_residuals.fGARCH = function(y,normality = "epps",
                                  unit_root = NULL,
                                  seasonal = NULL,
                                  arch = NULL,
                                  alpha = 0.05,plot = FALSE,...){
  if( !is(y,class2 = "fGARCH") )
    stop("The object y is not a fGARCH class")

  check_residuals.numeric(y@residuals,unit_root = unit_root,
                          normality = normality,
                          seasonal = seasonal,arch = arch,
                          alpha = 0.05,
                          plot = plot,...)
}
#'
#' @method check_residuals lm
#' @export
#'
check_residuals.lm = function(y,normality = "epps",
                              unit_root = NULL,
                              seasonal = NULL,
                              arch = NULL,
                              alpha = 0.05,plot = FALSE,...){
  if( !is(y,class2 = "lm") )
    stop("The object y is not a lm class")

  check_residuals.numeric(y$residuals,unit_root = unit_root,
                          normality = normality,
                          seasonal = seasonal,arch = arch,
                          alpha = 0.05,
                          plot = plot,...)
}
#' @method check_residuals HoltWinters
#' @export
#'
check_residuals.HoltWinters = function(y,normality = "epps",
                                       unit_root = NULL,
                                       seasonal = NULL,
                                       arch = NULL,
                                       alpha = 0.05,plot = FALSE,...){
  if( !is(y,class2 = "HoltWinters") )
    stop("The object y is not an HoltWinters class")

  y1 = residuals(y)

  check_residuals.ts(y1,unit_root = unit_root,
                     normality = normality,
                     seasonal = seasonal,arch = arch,
                     alpha = 0.05,
                     plot = plot,...)
}
#'
#' @method check_residuals ets
#' @export
#'
check_residuals.ets = function(y,normality = "epps",
                               unit_root = NULL,
                               seasonal = NULL,
                               arch = NULL,
                               alpha = 0.05,plot = FALSE,...){
  if( !is(y,class2 = "ets") )
    stop("The object y is not a ets class")

  check_residuals.ts(y$residuals,
                     unit_root = unit_root,
                     normality = normality,
                     seasonal = seasonal,arch = arch,
                     alpha = 0.05,
                     plot = plot,...)
}
#'
#' @method check_residuals forecast
#' @export
#'
check_residuals.forecast = function(y,normality = "epps",
                                    unit_root = NULL,
                                    seasonal = NULL,
                                    arch = NULL,
                                    alpha = 0.05,plot = FALSE,...){
  if( !is(y,class2 = "forecast") )
    stop("The object y is not a forecast class")

  check_residuals.ts(y$residuals,
                     unit_root = unit_root,
                     normality = normality,
                     seasonal = seasonal,arch = arch,
                     alpha = 0.05,
                     plot = plot,...)
}
