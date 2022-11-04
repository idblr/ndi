## This is the fourth resubmission

* Actions taken since previous submission:
  * Added `duncan()` function to compute the Dissimilarity Index (DI) based on [Duncan & Duncan (1955)](https://doi.org/10.2307/2088328) for specified counties/tracts 2009 onward.
  * Fixed bug in `bravo()` function where ACS-5 data (2005-2009) are from the "B15002" question and "B06009" after
  * Fixed bug in missingness warning for all metrics
  * `utils` is now Imports
  * Updated vignette and README with new features
  * Updated Description in DESCRIPTION
  * Updated tests
  * Updated CITATION with new citation for the additional metric

* Documentation for DESCRIPTION, README, NEWS, and vignette references the following DOIs, which throws a NOTE but are a valid URL:
  * <https://doi.org/10.1111/j.1749-6632.2009.05333.x>
  * <https://doi.org/10.2307/2223319>
  * <https://doi.org/10.2307/2088328>
  
* Some tests and examples for `anthopolos()`, `bravo()`, `duncan()`, `gini()`, `krieger()`, `messer()`, and `powell_wiley()` functions require a Census API key so they are skipped if NULL or not run

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
