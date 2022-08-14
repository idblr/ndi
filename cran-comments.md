## This is the second resubmission

* Actions taken since previous submission:
* Added `anthopolos()` function to compute the Racial Isolation Index (RI) based on based on [Anthopolos et al. (2011)](https://www.doi.org/10.1016/j.sste.2011.06.002) for specified counties/tracts 2009-2020
* Added `bravo()` function to compute the Educational Isolation Index (EI) based on based on [Bravo et al. (2021)](https://www.doi.org/10.3390/ijerph18179384) for specified counties/tracts 2009-2020
* Added `gini()` function to retrieve the Gini Index based on [Gini (1921)](https://www.doi.org/10.2307/2223319) for specified counties/tracts 2009-2020
* `Matrix` and `sf` are now Depends
* Updated vignette and README for new features
* Fixed typos throughout documentation
* Updated Description in DESCRIPTION
* Updated 'package.R' with new details and section
* Updated CITATION with new citations for the additional metrics

* Documentation for DESCRIPTION, README, NEWS, and vignette references the following DOIs, which throws a NOTE but are a valid URL:
  * https://doi.org/10.1111/j.1749-6632.2009.05333.x
  * https://doi.org/10.2307/2223319
  
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
