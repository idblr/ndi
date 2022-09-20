# ndi (development version)

# ndi v0.1.2.9002

### New Features
* Added `krieger()` function to compute the Index of Concentration at the Extremes (ICE) based on [Feldman et al. (2015)](https://www.doi.org/10.1136/jech-2015-205728) and [Krieger et al. (2016)](https://www.doi.org/10.2105/AJPH.2015.302955) for specified counties/tracts 2009 onward. 
* Thank you for the feature suggestion, [David Berrigan](https://orcid.org/0000-0002-5333-179X).
* Added `df` argument for the `messer()` and `powell_wiley()` functions to specify a pre-formatted dataset input for the NDI computation
* Added `round_output` argument for the `messer()` and `powell_wiley()` functions to provide raw output as the default and rounded output as optional.
* Thank you for the suggested enhancements, [Chris Prener](https://github.com/chris-prener)
* Added `DCtracts2020` a testing dataset for the `ndi` package and its documentation

### Updates
* Fixed bug in `powell_wiley()` function where the internal PCA will now run properly if only one factor has an eigenvalue above 1 
* Optimized the code to calculate missingness in all functions
* Thanks you for the suggested bug fixes, [Jacob Englert](https://github.com/jacobenglert)
* Relaxed `year` argument in functions to include any year after 2009 or 2010 for the indices
* Cleaned-up output formatting in functions
* `usethis` is now Suggests and `LazyData` is set to 'true'
* Updated tests for the `df` argument in `messer()` and `powell_wiley()` functions
* Updated vignette and README with new features
* Fixed typos throughout documentation
* Updated Description in DESCRIPTION and fixed typos
* Updated 'package.R' with new details
* Updated CITATION with new citations for the additional metric

# ndi v0.1.1

### New Features
* Added `anthopolos()` function to compute the Racial Isolation Index (RI) based on based on [Anthopolos et al. (2011)](https://www.doi.org/10.1016/j.sste.2011.06.002) for specified counties/tracts 2009 onward
* Added `bravo()` function to compute the Educational Isolation Index (EI) based on based on [Bravo et al. (2021)](https://www.doi.org/10.3390/ijerph18179384) for specified counties/tracts 2009 onward
* Added `gini()` function to retrieve the Gini Index based on [Gini (1921)](https://www.doi.org/10.2307/2223319) for specified counties/tracts 2009 onward
* Thank you for the suggestions, [Jessica Madrigal](https://orcid.org/0000-0001-5303-5109).

### Updates
* `Matrix` and `sf` are now Depends
* Updated vignette and README for new features
* Fixed typos throughout documentation
* Updated Description in DESCRIPTION
* Updated 'package.R' with new details and section
* Updated CITATION with new citations for the additional metrics

# ndi v0.1.0
* Fixed invalid URL and typos in package README.md

# ndi v0.0.1
* Initial CRAN submission
