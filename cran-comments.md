## This is the seventh resubmission

* Actions taken since previous submission:
  * Added `duncan_cuzzort()` function to compute the aspatial racial or ethnic Absolute Centralization (*ACE*) based on Duncan, Cuzzort, & Duncan (1961; LC:60007089) and [Massey & Denton (1988)](https://doi.org/10.1093/sf/67.2.281)
  * Added `hoover()` function to compute the aspatial racial or ethnic Delta (*DEL*) based on [Hoover (1941)](https://doi.org/10.1017/S0022050700052980) and Duncan, Cuzzort, & Duncan (1961; LC:60007089)
  * Added `james_taeuber()` function to compute the aspatial racial or ethnic Dissimilarity Index (*D*) based on [James & Taeuber (1985)](https://doi.org/10.2307/270845)
  * Added `lieberson()` function to compute the aspatial racial or ethnic Isolation Index (_xPx\*_) based on Lieberson (1981; ISBN-13:978-1-032-53884-6) and and [Bell (1954)](https://doi.org/10.2307/2574118)
  * Added `theil()` function the aspatial racial or ethnic Entropy (*H*) based on Theil (1972; ISBN:978-0-444-10378-9) and [Theil & Finizza (1971)](https://doi.org/110.1080/0022250X.1971.9989795)
  * Added `white_blau()` function to compute an index of spatial proximity (*SP*) based on [White (1986)](https://doi.org/10.2307/3644339) and Blau (1977; ISBN-13:978-0-029-03660-0)
  * Added `geo_large = 'place'` for census-designated places, `geo_large = 'cbsa'` for core-based statistical areas, `geo_large = 'csa'` for combined statistical areas, and `geo_large = 'metro'` for metropolitan divisions as the larger geographical unit in `atkinson()`, `bell()`, `bemanian_beyer()`, `duncan()`, `duncan_cuzzort()`, `hoover()`, `james_taeuber()`, `lieberson()`, `sudano()`, `theil()`, and `white()`, `white_blau()` functions.
  * Added census block group computation for `anthopolos()` by specifying `geo == 'cbg'` or `geo == 'block group'`
  * Added `holder` argument to `atkinson()` function to toggle the computation with or without the Hölder mean. The function can now compute *A* without the Hölder mean. The default is `holder = FALSE`.
  * Added `crs` argument to `anthopolos()`, `bravo()`, and `white_blau()` functions to provide spatial projection of the distance-based metrics
  * The `gini()` function now computes the aspatial racial or ethnic Gini Index (*G*) based on [Gini (1921)](https://doi.org/10.2307/2223319) as the main outcome. Arguments `geo_large`, `geo_small`, `subgroup`, and `omit_NAs` were added and argument `geo` was deprecated. The `gini()` function still retrieves the original output of the aspatial income Gini Index (*G*) at each smaller geography and is moved from the `g` output to `g_data` output.
  * Specifying census block groups in `geo` or `geo_small` arguments is now `'block group'` or `'cbg'` to match internal `get_acs()` function from the [tidycensus](https://CRAN.R-project.org/package=tidycensus) package
  * `bell()` function computes the Interaction Index (Bell) not the Isolation Index as previously documented. Updated documentation throughout.
  * Fixed bug in `bell()`, `bemanian_beyer()`, `duncan()`, `sudano()`, and `white()` functions when a smaller geography contains n=0 total population, will assign a value of zero (0) in the internal calculation instead of NA
  * Fixed bug in `atkinson()` function to properly compute the income Atkinson Index
  * Renamed *AI* as *A*, *DI* as *D*, *Gini* as *G*, and *II* as _xPy\*_ to align with the definitions from [Massey & Denton (1988)](https://doi.org/10.1093/sf/67.2.281). The output for `atkinson()` now produces `a` instead of `ai`. The output for `duncan()` now produces `d` instead of `ai`. The output for `gini()` now produces `g` instead of `gini`. The output for `bell()` now produces `xPy_star` instead of `II`. The internal functions `ai_fun()`, `di_fun()` and `ii_fun()` were renamed `a_fun()`, `ddd_fun()` and `xpy_star_fun()`, respectively.
  * `tigris` and `units` are now Imports
  * Reformatted functions for consistent internal structure
  * 'package.R' deprecated. Replaced with 'ndi-package.R' and reordered the contents
  * Consolidated DESCRIPTION
  * Re-formatted code and documentation throughout for consistent readability
  * Renamed 'race/ethnicity' or 'racial/ethnic' to 'race or ethnicity' or 'racial or ethnic' throughout documentation to use more modern, inclusive, and appropriate language
  * Updated documentation about value range of *V* (White) from `{0 to 1}` to `{-Inf to Inf}`
  * Split up vignette into three separate vignettes: 'ndi1', 'ndi2', and 'ndi3' for the *NDI*, racial or ethnic residential segregation, and additional socioeconomic disparity indices, respectively
  * Added examples for `atkinson()`, `duncan_cuzzort()`, `gini()`, `hoover()`, `james_taeuber()`, `lieberson()`, `theil()`, and `white_blau()` functions in vignettes and README
  * Added example for `holder` argument in `atkinson()` function in README
  * Reordered the README examples alphabetically
  * Reordered the vignette examples to group the racial or ethnic residential segregation indices
  * Updated examples in vignettes to showcase a larger variety of U.S. states
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
  
* Some tests and examples for `anthopolos()`, `atkinson()`, `bell()`, `bemanian_beyer()`, `bravo()`, `duncan()`, `gini()`, `hoover()`, `james_taeuber()`, `krieger()`, `lieberson()`, `messer()`, `powell_wiley()`, `sudano()`, `theil()`, `white()`, and `white_blau()` functions require a Census API key so they are skipped if NULL or not run

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
