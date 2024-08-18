## This is the seventh resubmission

* Actions taken since previous submission:
  * Added `hoover()` function to compute the aspatial racial/ethnic Delta (*DEL*) based on [Hoover (1941)](https://doi.org/10.1017/S0022050700052980) and Duncan et al. (1961; LC:60007089)
  * Added  `geo_large = 'cbsa'` option for computing Core Based Statistical Areas as the larger geographical unit in `atkinson()`, `bell()`, `bemanian_beyer()`, `duncan()`, `hoover()`, `sudano()`, and `white()` functions.
  * Thank you for the feature suggestions, [Symielle Gaston](https://orcid.org/0000-0001-9495-1592)
  * Fixed bug in `bell()`, `bemanian_beyer()`, `duncan()`, `sudano()`, and `white()` when a smaller geography contains n=0 total population, will assign a value of zero (0) in the internal calculation instead of NA
  * `tigris` is now Imports
  * 'package.R' deprecated. Replaced with 'ndi-package.R'
  * Re-formatted code and documentation throughout for consistent readability
  * Updated documentation about value range of *V* (White) from `{0 to 1}` to `{-Inf to Inf}`
  * Updated examples in vignette (& README) an example for `hoover()` and a larger variety of U.S. states
  * Updated documentation formatting of metric names in most functions

* Documentation for DESCRIPTION, README, NEWS, and vignette references the following DOIs, which throws a NOTE but are a valid URL:
  * <https://doi.org/10.1111/j.1749-6632.2009.05333.x>
  * <https://doi.org/10.2307/2223319>
  * <https://doi.org/10.2307/2088328>
  * <https://doi.org/10.2307/270845>
  * <https://doi.org/10.1080/17445647.2020.1750066>
  * <https://doi.org/10.2307/3644339>
  * <https://doi.org/10.2307/2084686>
  
* Some tests and examples for `anthopolos()`, `atkinson()`, `bell()`, `bemanian_beyer()`, `bravo()`, `duncan()`, `gini()`, `hoover()`, `krieger()`, `messer()`, `powell_wiley()`, `sudano()`, and `white()` functions require a Census API key so they are skipped if NULL or not run

## Test environments
* local Windows install, R 4.4.1
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
