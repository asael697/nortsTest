**Packages NEWS and Updates**
============

**nortsTest 1.0.3 Date: 12/06/2021**
----------------------------------

### Features:

### Changes:

### Fixes:
- Fix documentation.
- Update p-values mixing method using a False discovery rates, and Benjamin and 
  Yekuteli (2001) procedures.
- update the random projections procedure for the rp.test.

----------------------------------

**nortsTest 1.0.1 Date: 17/06/2021**
----------------------------------

### Features:

- Add GPL-2 license.

### Changes:

- False discovery rate set as default for the random projections test

- Change "NA %in$ y" to "anyNA(y)"

- use match.args() function for more flexibility

### Fixes:

- Bug fix for the False discovery rate.


**nortsTest 1.0.0 Date: 08/07/2020**
----------------------------------

### Fixes:

- Bug fix for CRAN submission.


**nortsTest 1.0.0 Date: 12/06/2020**
----------------------------------

### Fixes:

- Bug fix for CRAN submission.


**nortsTest 1.0.0 Date: 09/05/2020**
----------------------------------

### Features:

- vavra.test() function for the Psaradaki and Vavra test.

- sieve.bootstrap() function for bootstrap sub-sample in stationary time series.

- vavra.sample() function for the Anderson Darling sample statistics for the Psaradaki and Vavra test.

- rp.test() function for the random projections test.

- rp.sample() function for the random projections statistics samples.

- random.projection() function to generate a random projection of a stationary process.

- gghist(), ggnorm(), ggacf() and ggpacf() function for visualization.

### Improvements:

- epps.statistic() bug fixed.

- package documentation, errors, warnings and notes corrected.

### Changes:

- plot.compare() function deleted

- epps.statistic() using the **poweR** package deleted.

### Fixes:

-   The amoebam algorithm for epps.statistc() is corrected.


**normality 0.0.1.000 Date: 12/03/2020**
----------------------------------

### Features:

-  normal.test() function  with all the normality test available

-  uroot.test() function with all the unit root test available

-  seasonal.test() function with all the seasonal unit root test available

-  arch.test() function with the arch effect test for stationarity

- check_plot() methods implemented for residual diagnostics visualization for lm, glm, ets, forecast, Arima, arima0, fgarch, HoltWinters, ts, and numeric classes

- check_residuals() methods implemented for residual diagnostics for lm, glm, ets, forecast, Arima, arima0, fgarch, HoltWinters, ts, and numeric classes

- the autoplot() methods are overloaded for plotting time series (*ts*) and multivariate (*mts*) classes

### Improvements:

-   lobato.test() is implemented as a htest method

-   epps.test() is implemented as a htest method

### Changes:

-   LV.statistic() function change to lobato.statistic()

### Fixes:

-   The amoebam algorithm for Epps.statistc() is changed for the one implemented by the **PoweR** package


**normality 0.0.1.000 Date: 23/10/2020**
----------------------------------

### Features:

-  Updated the epps.statistic function

-  Updated the LV.statistic function

-  Perform documentation and package description


**normality 0.0.1.000 Date: 22/10/2020**
----------------------------------

### Features:

-  Package Repository created

-  Discussion of the previous work and code homogeneity check README.md

-  Incorporation of references and latex thesis algorithm

-  Incorporation of the article structure
