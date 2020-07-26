## Test environments
* local R installation, R 4.0.0
* ubuntu 16.04 (on travis-ci), R 4.0.0
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a new release.

## Comments previous CRAN 

- cat() function is used only in check_residuals() the objctive of 
this function is to **print** a summary of the test applyed to the data. 
Therefore, it can not be replaced for warning or message.


- The theoric reference for most of the tests was included already except for the random projections test, and that is because the DOI is not available and the link gives an error message.

