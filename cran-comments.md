## This is the fifth resubmission

* Actions taken since previous submission:
  * Added `atkinson()` function to compute the Atkinson Index (AI) based on [Atkinson (1970)](https://doi.org/10.1016/0022-0531(70)90039-6) for specified counties/tracts 2009 onward
  * Fixed bug in reverse dependency check failure for `anthopolos()` and `bravo()` functions removing `returnValue()` when data are not missing
  * Thank you, [Roger Bivand](https://github.com/rsbivand), for the catch. Relates to [ndi Issue #5](https://github.com/idblr/ndi/issues/5)
  * Updated `duncan()`, `gini()`, `krieger()`, `messer()`, and `powell_wiley()` for consistency in messaging when data are not missing
  * Updated tests for `anthopolos()` and `bravo()` if `Sys.getenv("CENSUS_API_KEY") != ""`
  * Added `omit_NAs` argument in `duncan()` function to choose if NA values will be included in its computation
  * In `duncan()` function, if any smaller geographic unit has zero counts the output for its larger geographic unit will be NA
  * Fixed bug in `duncan()` function for multiple `subgroup` and `subgroup_ref` selections
  * Updated documentation throughout

* Documentation for DESCRIPTION, README, NEWS, and vignette references the following DOIs, which throws a NOTE but are a valid URL:
  * <https://doi.org/10.1111/j.1749-6632.2009.05333.x>
  * <https://doi.org/10.2307/2223319>
  * <https://doi.org/10.2307/2088328>
  * <https://doi.org/10.2307/270845>
  * <https://doi.org/10.1080/17445647.2020.1750066>
  
* Some tests and examples for `anthopolos()`, `atkinson()`, `bravo()`, `duncan()`, `gini()`, `krieger()`, `messer()`, and `powell_wiley()` functions require a Census API key so they are skipped if NULL or not run

## Test environments
* local Windows install, R 4.2.1
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
