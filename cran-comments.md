## This is the seventh resubmission

* Actions taken since previous submission:
  * Added `hoover()` function to compute the aspatial racial/ethnic Delta (*DEL*) based on [Hoover (1941)](https://doi.org/10.1017/S0022050700052980) and Duncan et al. (1961; LC:60007089)
  * Added `white_blau()` function to compute an index of spatial proximity (*SP*) based on [White (1986)](https://doi.org/10.2307/3644339) and Blau (1977; ISBN-13:978-0-029-03660-0)
  * Added `lieberson()` function to compute he aspatial racial/ethnic Isolation Index (_xPx\*_) based on [White (1986)](https://doi.org/10.2307/3644339) and Blau (1977; ISBN-13:978-0-029-03660-0)
  * Added `geo_large = 'cbsa'` for Core Based Statistical Areas, `geo_large = 'csa'` for Combined Statistical Areas, and `geo_large = 'metro'` for Metropolitan Divisions as the larger geographical unit in `atkinson()`, `bell()`, `bemanian_beyer()`, `duncan()`, `hoover()`, `sudano()`, and `white()`, `white_blau()` functions.
  * Added `holder` argument to `atkinson()` function to toggle the computation with or without the Hölder mean. The function can now compute *A* without the Hölder mean. The default is `holder = TRUE`.
  * `bell()` function computes the Interaction Index (Bell) not the Isolation Index as previously documented. Updated documentation throughout
  * Fixed bug in `bell()`, `bemanian_beyer()`, `duncan()`, `sudano()`, and `white()` functions when a smaller geography contains n=0 total population, will assign a value of zero (0) in the internal calculation instead of NA
  * Renamed *AI* as *A*, *DI* as *D*, *Gini* as *G*, and *II* as _xPy\*_ to align with the definitions from [Massey & Denton (1988)](https://doi.org/10.1093/sf/67.2.281). The output for `atkinson()` now produces `a` instead of `ai`. The output for `duncan()` now produces `d` instead of `ai`. The output for `gini()` now produces `g` instead of `gini`. The output for `bell()` now produces `xPy_star` instead of `II`. The internal functions `ai_fun()`, `di_fun()` and `ii_fun()` were renamed `a_fun()`, `d_fun()` and `xpy_star_fun()`, respectively.
  * `tigris` and `units` are now Imports
  * 'package.R' deprecated. Replaced with 'ndi-package.R'
  * Re-formatted code and documentation throughout for consistent readability
  * Updated documentation about value range of *V* (White) from `{0 to 1}` to `{-Inf to Inf}`
  * Add examples for `hoover()` and `white_blau()` functions in vignette and README
  * Add example for `holder` argument in `atkinson()` function in README
  * Reformatted functions for consistent internal structure
  * Updated examples in vignette to showcase a larger variety of U.S. states
  * Updated examples in functions to better describe the metrics
  * Updated documentation formatting of metric names in all functions

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
