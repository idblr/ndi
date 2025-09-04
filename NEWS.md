# ndi (development version)

## ndi v0.2.1

### New Features
* None

### Updates
* Fixed broken URLs in 'theil.Rd', 'ndi-package.Rd', 'ndi2.html', README, and NEWS

## ndi v0.2.0

### New Features

#### New Functions
* Added `duncan_cuzzort()` function to compute the aspatial racial or ethnic Absolute Centralization (*ACE*) based on Duncan, Cuzzort, & Duncan (1961; LC:60007089) and [Massey & Denton (1988)](https://doi.org/10.1093/sf/67.2.281)
* Added `hoover()` function to compute the aspatial racial or ethnic Delta (*DEL*) based on [Hoover (1941)](https://doi.org/10.1017/S0022050700052980) and Duncan, Cuzzort, & Duncan (1961; LC:60007089)
* Added `james_taeuber()` function to compute the aspatial racial or ethnic Dissimilarity Index (*D*) based on [James & Taeuber (1985)](https://doi.org/10.2307/270845)
* Added `lieberson()` function to compute the aspatial racial or ethnic Isolation Index (_xPx\*_) based on Lieberson (1981; ISBN-13:978-1-032-53884-6) and and [Bell (1954)](https://doi.org/10.2307/2574118)
* Added `theil()` function the aspatial racial or ethnic Entropy (*H*) based on Theil (1972; ISBN:978-0-444-10378-9) and [Theil & Finizza (1971)](https://doi.org/10.1080/0022250X.1971.9989795)
* Added `white_blau()` function to compute an index of spatial proximity (*SP*) based on [White (1986)](https://doi.org/10.2307/3644339) and Blau (1977; ISBN-13:978-0-029-03660-0)
* Thank you for the feature suggestions above, [Symielle Gaston](https://orcid.org/0000-0001-9495-1592)
* Added `denton()` function to compute the aspatial racial or ethnic Relative Clustering (*RCL*) based on [Massey & Denton (1988)](https://doi.org/10.1093/sf/67.2.281)
* Added `denton_cuzzort()` function to compute the aspatial racial or ethnic Relative Concentration (*RCO*) based on [Massey & Denton (1988)](https://doi.org/10.1093/sf/67.2.281) and Duncan, Cuzzort, & Duncan (1961; LC:60007089)
* Added `duncan_duncan()` function to compute the aspatial racial or ethnic Relative Centralization (*RCE*) based on [Duncan & Duncan (1955b)](https://doi.org/10.1086/221609) and [Massey & Denton (1988)](https://doi.org/10.1093/sf/67.2.281)
* Added `massey()` function to compute the aspatial racial or ethnic Absolute Clustering (*ACL*) based on [Massey & Denton (1988)](https://doi.org/10.1093/sf/67.2.281)
* Added `massey_duncan()` function to compute the aspatial racial or ethnic Absolute Concentration (*ACO*) based on [Massey & Denton (1988)](https://doi.org/10.1093/sf/67.2.281) and Duncan, Cuzzort, & Duncan (1961; LC:60007089)
* Added `morgan_denton()` function to compute the aspatial racial or ethnic Distance-Decay Interaction Index (_DPxy\*_) based on [Morgan (1983)](https://www.jstor.org/stable/20001935) and [Massey & Denton (1988)](https://doi.org/10.1093/sf/67.2.281)
* Added `morgan_massey()` function to compute the aspatial racial or ethnic Distance-Decay Isolation Index (_DPxx\*_) based on [Morgan (1983)](https://www.jstor.org/stable/20001935) and [Massey & Denton (1988)](https://doi.org/10.1093/sf/67.2.281)

#### New Function Capabilities
* Added `geo_large = 'place'` for census-designated places, `geo_large = 'cbsa'` for core-based statistical areas, `geo_large = 'csa'` for combined statistical areas, and `geo_large = 'metro'` for metropolitan divisions as the larger geographical unit in `atkinson()`, `bell()`, `bemanian_beyer()`, `denton()`, `denton_cuzzort()`, `duncan()`, `duncan_cuzzort()`, `duncan_duncan()`, `hoover()`, `james_taeuber()`, `lieberson()`, `massey()`, `massey_duncan()`, `morgan_denton()`, `morgan_denton()`, `morgan_massey()`, `sudano()`, `theil()`, and `white()`, `white_blau()` functions.
* Added census block group computation for `anthopolos()` by specifying `geo == 'cbg'` or `geo == 'block group'`
* Added `holder` argument to `atkinson()` function to toggle the computation with or without the Hölder mean. The function can now compute *A* without the Hölder mean. The default is `holder = FALSE`.
* Added `crs` argument to `anthopolos()`, `bravo()`, and `white_blau()` functions to provide spatial projection of the distance-based metrics
* The `gini()` function now computes the aspatial racial or ethnic Gini Index (*G*) based on [Gini (1921)](https://doi.org/10.2307/2223319) as the main outcome. Arguments `geo_large`, `geo_small`, `subgroup`, and `omit_NAs` were added and argument `geo` was deprecated. The `gini()` function still retrieves the original output of the aspatial income Gini Index (*G*) at each smaller geography and is moved from the `g` output to `g_data` output.
* Specifying census block groups in `geo` or `geo_small` arguments is now `'block group'` or `'cbg'` to match internal `get_acs()` function from the [tidycensus](https://CRAN.R-project.org/package=tidycensus) package

### Updates

#### Bug Fixes
* Fixed NOTE in CRAN checks to provide package anchors for Rd `\link{}` targets not in the package itself and the base packages  within 'ndi-package.Rd'
* Updated population-weighted quantile method from `stats::quantile` to `Hmisc::wtd.quantile` in `powell_wiley()` thanks to a contribution (#32) by [Hunter Miller](https://github.com/huntermills707)
* `bell()` function computes the Interaction Index (Bell) not the Isolation Index as previously documented. Updated documentation throughout.
* Fixed bug in `bell()`, `bemanian_beyer()`, `duncan()`, `sudano()`, and `white()` functions when a smaller geography contains n=0 total population, will assign a value of zero (0) in the internal calculation instead of NA
* Fixed bug in `atkinson()` function to properly compute the income Atkinson Index
* Renamed *AI* as *A*, *DI* as *D*, *Gini* as *G*, and *II* as _xPy\*_ to align with the definitions from [Massey & Denton (1988)](https://doi.org/10.1093/sf/67.2.281). The output for `atkinson()` now produces `a` instead of `ai`. The output for `duncan()` now produces `d` instead of `ai`. The output for `gini()` now produces `g` instead of `gini`. The output for `bell()` now produces `xPy_star` instead of `II`. The internal functions `ai_fun()`, `di_fun()` and `ii_fun()` were renamed `a_fun()`, `ddd_fun()` and `xpy_star_fun()`, respectively.
* 'package.R' deprecated. Replaced with 'ndi-package.R'
* Output of racial or ethnic residential segregation indices is now rounded to four significant digits

#### New Dependencies
* `Hmisc`, `tigris`, and `units` are now Imports

#### Updated Documentation
* Split up vignette into three separate vignettes: 'ndi1', 'ndi2', and 'ndi3' for the *NDI*, racial or ethnic residential segregation, and additional socioeconomic disparity indices, respectively
* Consolidated DESCRIPTION
* Reformatted functions for consistent internal structure
* Re-formatted code and documentation throughout for consistent readability
* Renamed 'race/ethnicity' or 'racial/ethnic' to 'race or ethnicity' or 'racial or ethnic' throughout documentation to use more modern, inclusive, and appropriate language
* Updated documentation about value range of *V* (White) from `{0 to 1}` to `{-Inf to Inf}`
* Added examples for `atkinson()`, `duncan_cuzzort()`, `duncan_duncan()`, `gini()`, `hoover()`, `james_taeuber()`, `lieberson()`, `massey()`, `massey_duncan()`, `morgan_massey()`, `theil()`, and `white_blau()` functions in vignettes and README
* Added example for `holder` argument in `atkinson()` function in README
* Added internal and external references within each functions documentation
* Reordered and reformatted CITATION alphabetically by index abbreviation
* Reordered the contents of 'ndi-package.R' thematically
* Reordered the README examples alphabetically
* Reordered the vignette examples to group the racial or ethnic residential segregation indices
* Updated examples in vignettes to showcase a larger variety of U.S. states
* Updated examples in functions to better describe the metrics
* Updated documentation formatting of metric names in all functions

## ndi v0.1.5

### New Features
* None

### Updates
* 'DescTools' is now Suggests to fix Rd cross-references NOTE
* Fixed 'lost braces in \itemize' NOTE for `anthopolos()`, `atkinson()`, `bell()`, `bemanian_beyer()`, `bravo()`, `duncan()`, `krieger()`, `messer()`, `powell_wiley()`, `sudano()`, and `white()` functions
* Fixed 'Moved Permanently' content by replacing the old URL with the new URL
* Fixed citation for [Slotman et al. (2022)](https://doi.org/10.1016/j.dib.2022.108002) in CITATION

## ndi v0.1.4

### New Features
* Added `atkinson()` function to compute the aspatial income or racial or ethnic Atkinson Index (*A*) based on [Atkinson (1970)](https://doi.org/10.1016/0022-0531(70)90039-6) for specified counties/tracts 2009 onward
* Added `bell()` function to compute the aspatial racial or ethnic Interaction Index (_xPy\*_) based on Shevky & Williams (1949; ISBN-13:978-0837156378) and [Bell (1954)](https://doi.org/10.2307/2574118)
* Added `white()` function to compute the aspatial racial or ethnic Correlation Ratio (*V*) based on [Bell (1954)](https://doi.org/10.2307/2574118) and [White (1986)](https://doi.org/10.2307/3644339)
* Added `sudano()` function to compute the aspatial racial or ethnic Location Quotient (*LQ*) based on [Merton (1939)](https://doi.org/10.2307/2084686) and [Sudano et al. (2013)](https://doi.org/10.1016/j.healthplace.2012.09.015)
* Added `bemanian_beyer()` function to compute the aspatial racial or ethnic Local Exposure and Isolation (*LEx/Is*) metric based on [Bemanian & Beyer (2017)](https://doi.org/10.1158/1055-9965.EPI-16-0926)

### Updates
* `car` is now Imports
* Fixed bug in reverse dependency check failure for `anthopolos()` and `bravo()` functions removing `returnValue()` when data are not missing
* Thank you, [Roger Bivand](https://github.com/rsbivand), for the catch. Relates to [*ndi* Issue #5](https://github.com/idblr/ndi/issues/5)
* Updated `duncan()`, `gini()`, `krieger()`, `messer()`, and `powell_wiley()` for consistency in messaging when data are not missing
* Updated tests for `anthopolos()` and `bravo()` if `Sys.getenv('CENSUS_API_KEY') != ''`
* Added `omit_NAs` argument in `duncan()` function to choose if NA values will be included in its computation
* In `duncan()` function, if any smaller geographic unit has zero counts the output for its larger geographic unit will be NA
* Fixed bug in `duncan()` function for multiple `subgroup` and `subgroup_ref` selections
* Updated documentation throughout
* Added GitHub R-CMD-check
* Updated citation style for CITATION file

## ndi v0.1.3

### New Features
* Added `duncan()` function to compute the Dissimilarity Index (*D*) based on [Duncan & Duncan (1955a)](https://doi.org/10.2307/2088328) for specified counties/tracts 2009 onward
* Thank you for the feature suggestion, [Jessica Madrigal](https://orcid.org/0000-0001-5303-5109)
* Added 'utils.R' file with internal `di_fun()` function for `duncan()` function

### Updates
* Fixed bug in `bravo()` function where ACS-5 data (2005-2009) are from the 'B15002' question and 'B06009' after
* Fixed bug in missingness warning for all metrics
* `utils` is now Imports
* Updated vignette and README with new features
* Updated Description in DESCRIPTION
* Updated tests
* Updated CITATION with new citation for the additional metric
* Updated maintainer contact information

## ndi v0.1.2

### New Features
* Added `krieger()` function to compute the Index of Concentration at the Extremes (*ICE*) based on [Feldman et al. (2015)](https://doi.org/10.1136/jech-2015-205728) and [Krieger et al. (2016)](https://doi.org/10.2105/AJPH.2015.302955) for specified counties/tracts 2009 onward 
* Thank you for the feature suggestion, [David Berrigan](https://orcid.org/0000-0002-5333-179X)
* Added `df` argument for the `messer()` and `powell_wiley()` functions to specify a pre-formatted data set input for the NDI computation
* Added `round_output` argument for the `messer()` and `powell_wiley()` functions to provide raw output as the default and rounded output as optional.
* Thank you for the suggested enhancements, [Chris Prener](https://github.com/chris-prener)
* Added `DCtracts2020` a testing data set for the *ndi* package and its documentation

### Updates
* Fixed bug in `powell_wiley()` function where the internal PCA will now run properly if only one factor has an eigenvalue above 1 
* Optimized the code to calculate missingness in all functions
* Thank you for the suggested bug fixes, [Jacob Englert](https://github.com/jacobenglert)
* Fixed bug in `powell_wiley()` function where 'PctNoPhone' before 2015 is 'DP04_0074PE' and 'DP04_0075PE' after
* Thank you for alerting this issue, [Jessica Gleason](https://orcid.org/0000-0001-9877-7931)
* Relaxed `year` argument in functions to include any year after 2009 or 2010 for the indices
* Cleaned-up output formatting in functions
* `usethis` is now Suggests and `LazyData` is set to 'true'
* Updated tests for the `df` argument in `messer()` and `powell_wiley()` functions
* Updated vignette and README with new features
* Fixed typos throughout documentation
* Updated Description in DESCRIPTION and fixed typos
* Updated 'package.R' with new details
* Updated CITATION with new citations for the additional metric

## ndi v0.1.1

### New Features
* Added `anthopolos()` function to compute the Racial Isolation Index (*RI*) based on based on [Anthopolos et al. (2011)](https://doi.org/10.1016/j.sste.2011.06.002) for specified counties/tracts 2009 onward
* Added `bravo()` function to compute the Educational Isolation Index (*EI*) based on based on [Bravo et al. (2021)](https://doi.org/10.3390/ijerph18179384) for specified counties/tracts 2009 onward
* Added `gini()` function to retrieve the Gini Index (*G*) based on [Gini (1921)](https://doi.org/10.2307/2223319) for specified counties/tracts 2009 onward
* Thank you for the feature suggestions, [Jessica Madrigal](https://orcid.org/0000-0001-5303-5109)

### Updates
* `Matrix` and `sf` are now Imports
* Updated vignette and README for new features
* Fixed typos throughout documentation
* Updated Description in DESCRIPTION
* Updated 'package.R' with new details and section
* Updated CITATION with new citations for the additional metrics

## ndi v0.1.0
* Fixed invalid URL and typos in package README.md

## ndi v0.0.1
* Initial CRAN submission
