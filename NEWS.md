# ndi (development version)

## ndi v0.1.6.9000

### New Features
* Added `hoover()` function to compute the aspatial racial/ethnic Delta (DEL) based on [Hoover (1941)](https://doi.org/10.1017/S0022050700052980) and Duncan et al. (1961; LC:60007089)
* Thank you for the feature suggestion, [Symielle Gaston](https://orcid.org/0000-0001-9495-1592)

### Updates
* Fixed bug in `bell()`, `bemanian_beyer()`, `duncan()`, `sudano()`, and `white()` when a smaller geography contains n=0 total population, will assign a value of zero (0) in the internal calculation instead of NA
* 'package.R' deprecated. Replaced with 'ndi-package.R'.

## ndi v0.1.5

### New Features
* None

### Updates
* 'DescTools' is now Suggests to fix Rd cross-references NOTE
* Fixed 'lost braces in \itemize' NOTE for `anthopolos()`, `atkinson()`, `bell()`, `bemanian_beyer()`, `bravo()`, `duncan()`, `krieger()`, `messer()`, `powell_wiley()`, `sudano()`, and `white()` functions
* Fixed 'Moved Permanently' content by replacing the old URL with the new URL
* Fixed citation for Slotman et al. (2022) in CITATION

## ndi v0.1.4

### New Features
* Added `atkinson()` function to compute the aspatial income or racial/ethnic Atkinson Index (AI) based on [Atkinson (1970)](https://doi.org/10.1016/0022-0531(70)90039-6) for specified counties/tracts 2009 onward
* Added `bell()` function to compute the aspatial racial/ethnic Isolation Index (II) based on Shevky & Williams (1949; ISBN-13:978-0837156378) and [Bell (1954)](https://doi.org/10.2307/2574118)
* Added `white()` function to compute the aspatial racial/ethnic Correlation Ratio (V) based on [Bell (1954)](https://doi.org/10.2307/2574118) and [White (1986)](https://doi.org/10.2307/3644339)
* Added `sudano()` function to compute the aspatial racial/ethnic Location Quotient (LQ) based on [Merton (1939)](https://doi.org/10.2307/2084686) and [Sudano et al. (2013)](https://doi.org/10.1016/j.healthplace.2012.09.015)
* Added `bemanian_beyer()` function to compute the aspatial racial/ethnic Local Exposure and Isolation (LEx/Is) metric based on [Bemanian & Beyer (2017)](https://doi.org/10.1158/1055-9965.EPI-16-0926)

### Updates
* `car` is now Imports
* Fixed bug in reverse dependency check failure for `anthopolos()` and `bravo()` functions removing `returnValue()` when data are not missing
* Thank you, [Roger Bivand](https://github.com/rsbivand), for the catch. Relates to [ndi Issue #5](https://github.com/idblr/ndi/issues/5)
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
* Added `duncan()` function to compute the Dissimilarity Index (DI) based on [Duncan & Duncan (1955)](https://doi.org/10.2307/2088328) for specified counties/tracts 2009 onward
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
* Added `krieger()` function to compute the Index of Concentration at the Extremes (ICE) based on [Feldman et al. (2015)](https://doi.org/10.1136/jech-2015-205728) and [Krieger et al. (2016)](https://doi.org/10.2105/AJPH.2015.302955) for specified counties/tracts 2009 onward 
* Thank you for the feature suggestion, [David Berrigan](https://orcid.org/0000-0002-5333-179X)
* Added `df` argument for the `messer()` and `powell_wiley()` functions to specify a pre-formatted data set input for the NDI computation
* Added `round_output` argument for the `messer()` and `powell_wiley()` functions to provide raw output as the default and rounded output as optional.
* Thank you for the suggested enhancements, [Chris Prener](https://github.com/chris-prener)
* Added `DCtracts2020` a testing data set for the `ndi` package and its documentation

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
* Added `anthopolos()` function to compute the Racial Isolation Index (RI) based on based on [Anthopolos et al. (2011)](https://doi.org/10.1016/j.sste.2011.06.002) for specified counties/tracts 2009 onward
* Added `bravo()` function to compute the Educational Isolation Index (EI) based on based on [Bravo et al. (2021)](https://doi.org/10.3390/ijerph18179384) for specified counties/tracts 2009 onward
* Added `gini()` function to retrieve the Gini Index based on [Gini (1921)](https://doi.org/10.2307/2223319) for specified counties/tracts 2009 onward
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
