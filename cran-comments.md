## This is the third resubmission

* Actions taken since previous submission:
  * Added `krieger()` function to compute the Index of Concentration at the Extremes (ICE) based on [Feldman et al. (2015)](https://www.doi.org/10.1136/jech-2015-205728) and [Krieger et al. (2016)](https://www.doi.org/10.2105/AJPH.2015.302955) for specified counties/tracts 2009-2020. 
  * Added `df` argument for the `messer()` and `powell_wiley()` functions to specify a pre-formatted dataset input for the NDI computation
  * Added `DCtracts2020` a testing dataset for the `ndi` package and its documentation
  * Fixed bug in `powell_wiley()` function where the internal PCA will now run properly if only one factor has an eigenvalue above 1 
  * Optimized the code to calculate missingness in all functions thanks to a contribution by 
  * Cleaned-up output formatting in functions
  * `usethis` is now Suggests and `LazyData` is set to 'true'
  * Updated tests for the `df` argument in `messer()` and `powell_wiley()` functions
  * Updated vignette and README for new features
  * Fixed typos throughout documentation
  * Updated Description in DESCRIPTION and fixed typos
  * Updated 'package.R' with new details
  * Updated CITATION with new citations for the additional metric

* Documentation for DESCRIPTION, README, NEWS, and vignette references the following DOIs, which throws a NOTE but are a valid URL:
  * <https://doi.org/10.1111/j.1749-6632.2009.05333.x>
  * <https://doi.org/10.2307/2223319>
  
* Some tests and examples for `anthopolos()`, `bravo()`, `gini()`, `messer()` and `powell_wiley()` functions require require a Census API key so they are skipped if NULL or not run

## Test environments
* local OS X install, R 4.2.1
* win-builder, (devel, release, oldrelease)
* Rhub
  * Fedora Linux, R-devel, clang, gfortran
  * Ubuntu Linux 20.04.1 LTS, R-release, GCC
  * Windows Server 2022, R-devel, 64 bit
  * Windows Server 2008 R2 SP1, R-release, 32⁄64 bit
  * Oracle Solaris 10, x86, 32 bit, R-release
  * macOS 10.13.6 High Sierra, R-release, CRAN's setup

## R CMD check results
0 errors | 0 warnings | 0 notes

## Submitted by Maintainer
