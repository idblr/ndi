## This is the fifth resubmission

* Actions taken since previous submission:
  * Added `atkinson()` function to compute the aspatial income or racial/ethnic Atkinson Index (AI) based on [Atkinson (1970)](https://doi.org/10.1016/0022-0531(70)90039-6) for specified counties/tracts 2009 onward
  * Added `bell()` function to compute the aspatial racial/ethnic Isolation Index (II) based on Shevky & Williams (1949; ISBN-13:978-0837156378) and [Bell (1954)](https://doi.org/10.2307/2574118)
  * Added `white()` function to compute the aspatial racial/ethnic Correlation Ratio based on [Bell (1954)](https://doi.org/10.2307/2574118) and [White (1986)](https://doi.org/10.2307/3644339)
  * Added `sudano()` function to compute the aspatial racial/ethnic Location Quotient based on [Merton (1939)](https://doi.org/10.2307/2084686) and [Sudano et al. (2013)](https://doi.org/10.1016/j.healthplace.2012.09.015)
  * Added `bemanian_beyer()` function to compute the aspatial racial/ethnic Local Exposure and Isolation metric based on [Bemanian & Beyer (2017)](https://doi.org/10.1158/1055-9965.EPI-16-0926)
  * `car` is now Imports
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
  * <https://doi.org/10.2307/3644339>
  * <https://doi.org/10.2307/2084686>
  
* Some tests and examples for `anthopolos()`, `atkinson()`, `bell()`, `bemanian_beyer()`, `bravo()`, `duncan()`, `gini()`, `krieger()`, `messer()`, `powell_wiley()`, `sudano()`, and `white()` functions require a Census API key so they are skipped if NULL or not run

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
