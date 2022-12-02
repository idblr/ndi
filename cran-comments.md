## This is the fifth resubmission

* Actions taken since previous submission:
  * Fixed bug in reverse dependency check failure for `anthopolos()` and `bravo()` functions removing `returnValue()` when data are not missing
  * Thank you, [Roger Bivand](https://github.com/rsbivand), for the catch. Relates to [ndi Issue #5](https://github.com/idblr/ndi/issues/5)
  * Updated `duncan()`, `gini()`, `krieger()`, `messer()`, and `powell_wiley()` for consistency in messaging when data are not missing
  * Updated tests for `anthopolos()` and `bravo()` if `Sys.getenv("CENSUS_API_KEY") != ""`

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
