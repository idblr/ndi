ndi: Neighborhood Deprivation Indices <img src="man/figures/ndi.png" width="120" align="right" />
===================================================

<!-- badges: start -->
<!-- [![CRAN
version](https://www.r-pkg.org/badges/version-ago/ndi)](https://cran.r-project.org/package=ndi)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/ndi?color=blue)](https://r-pkg.org/pkg/ndi) -->
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
<!-- badges: end -->

**Date repository last updated**: August 04, 2022

### Overview

The `ndi` package is a suite of `R` functions to compute various geospatial neighborhood deprivation indices (NDI) in the United States. Two types of NDI are available in the initial repository: (1) based on [Messer et al. (2006)](https://doi.org/10.1007/s11524-006-9094-x) and (2) based on [Andrews et al. (2020)](https://doi.org/10.1080/17445647.2020.1750066) and [Slotman et al. (2022)](https://doi.org/10.1016/j.dib.2022.108002) who uses variables chosen by [Roux and Mair (2010)](https://doi.org/10.1111/j.1749-6632.2009.05333.x). Both are a decomposition of various demographic characteristics from the U.S. Census Bureau American Community Survey 5-year estimates pulled by the [tidycensus](https://CRAN.R-project.org/package=tidycensus) package.

### Installation

To install the release version from CRAN:

    install.packages("ndi")

To install the development version from GitHub:

    devtools::install_github("idblr/ndi")

### Available functions

<table>
<colgroup>
<col width="30%" />
<col width="70%" />
</colgroup>
<thead>
<tr class="header">
<th>Function</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<td><code>messer</code></td>
<td>Compute NDI based on [Messer et al. (2006)](https://doi.org/10.1007/s11524-006-9094-x).</td>
</tr>
<td><code>powell_wiley</code></td>
<td>Compute NDI based on [Andrews et al. (2020)](https://doi.org/10.1080/17445647.2020.1750066) and [Slotman et al. (2022)](https://doi.org/10.1016/j.dib.2022.108002) with variables chosen by [Roux and Mair (2010)](https://doi.org/10.1111/j.1749-6632.2009.05333.x).</td>
</tr>
</tbody>
<table>

The repository also includes the code to create the project hexsticker.

### Author

* **Ian D. Buller** - *Occupational and Environmental Epidemiology Branch, Division of Cancer Epidemiology and Genetics, National Cancer Institute, National Institutes of Health, Rockville, Maryland* - [GitHub](https://github.com/idblr) - [ORCID](https://orcid.org/0000-0001-9477-8582)

See also the list of [contributors](https://github.com/idblr/ndi/graphs/contributors) who participated in this package.

### Getting Started

* Step 1: You must obtain a unique access key from the U.S. Census Bureau. Follow [this link](http://api.census.gov/data/key_signup.html) to obtain one.
* Step 2: Specify your access key in the `messer()` or `powell_wiley()` functions using the `key` argument or by using the `census_api_key()` function from the `tidycensus` package before running the `messer()` or `powell_wiley()` functions (see an example below).

### Usage

``` r
# ------------------ #
# Necessary packages #
# ------------------ #

library(ndi)
library(ggplot2)
library(sf)
library(tidycensus) # a dependency for the "ndi"" package
library(tigris) # a dependency for the "ndi"" package

# -------- #
# Settings #
# -------- #

## Access Key for census data download
### Obtain one at http://api.census.gov/data/key_signup.html
tidycensus::census_api_key("...") # INSERT YOUR OWN KEY FROM U.S. CENSUS API

# ---------------------- #
# Calculate NDI (Messer) #
# ---------------------- #

# Compute the NDI (Messer) values (2016-2020 5-year ACS) for Washington, D.C. census tracts
messer2020DC <- ndi::messer(state = "DC", year = 2020)

# ------------------------------ #
# Outputs from messer() function #
# ------------------------------ #

# A tibble containing the identification, geographic name, NDI (Messer) values, NDI (Messer) quartiles, and raw census characteristics for each tract
messer2020DC$ndi

# The results from the principal component analysis used to compute the NDI (Messer) values
messer2020DC$pca

# A tibble containing a breakdown of the missingingness of the census characteristics used to compute the NDI (Messer) values
messer2020DC$missing

# -------------------------------------- #
# Visualize the messer() function output #
# -------------------------------------- #

# Obtain the 2020 census tracts from the "tigris" package
tract2020DC <- tigris::tracts(state = "DC", year = 2020, cb = TRUE)

# Join the NDI (Messer) values to the census tract geometry
DC2020messer <- merge(tract2020DC, messer2020DC$ndi, by = "GEOID")

# Visualize the NDI (Messer) values (2016-2020 5-year ACS) for Washington, D.C. census tracts

## Continuous Index
ggplot2::ggplot() + 
  ggplot2::geom_sf(data = DC2020messer, 
                   ggplot2::aes(fill = NDI),
                   color = "white") +
  ggplot2::theme_minimal() + 
  ggplot2::scale_fill_viridis_c() +
  ggplot2::labs(fill = "Index (Continuous)",
                caption = "Source: U.S. Census ACS 2016-2020 estimates")+
  ggplot2::ggtitle("Neighborhood Deprivation Index\nContinuous (Messer, non-imputed)", subtitle = "Washington, D.C. tracts as the referent")

## Categorical Index (Quartiles)
### Rename "9-NDI not avail" level as NA for plotting
DC2020messer$NDIQuartNA <- factor(replace(as.character(DC2020messer$NDIQuart),
                                          DC2020messer$NDIQuart == "9-NDI not avail",
                                          NA),
                                  c(levels(DC2020messer$NDIQuart)[-5], NA))

ggplot2::ggplot() + 
  ggplot2::geom_sf(data = DC2020messer, 
                   ggplot2::aes(fill = NDIQuartNA),
                   color = "white") +
  ggplot2::theme_minimal() + 
  ggplot2::scale_fill_viridis_d(guide = ggplot2::guide_legend(reverse = TRUE),
                                na.value = "grey50") +
  ggplot2::labs(fill = "Index (Categorical)",
                caption = "Source: U.S. Census ACS 2016-2020 estimates") +
  ggplot2::ggtitle("Neighborhood Deprivation Index\nQuartiles (Messer, non-imputed)", subtitle = "Washington, D.C. tracts as the referent")
```
![](man/figures/messer1.png)
![](man/figures/messer2.png)

``` r
# ---------------------------- #
# Calculate NDI (Powell-Wiley) #
# ---------------------------- #

# Compute the NDI (Powell-Wiley) values (2016-2020 5-year ACS) for Washington, D.C. census tracts
powell_wiley2020DC <- powell_wiley(state = "DC", year = 2020)
powell_wiley2020DCi <- powell_wiley(state = "DC", year = 2020, imp = TRUE) # impute missing values

# ------------------------------------ #
# Outputs from powell_wiley() function #
# ------------------------------------ #

# A tibble containing the identification, geographic name, NDI (Powell-Wiley) value, and raw census characteristics for each tract
powell_wiley2020DC$ndi

# The results from the principal component analysis used to compute the NDI (Powell-Wiley) values
powell_wiley2020DC$fa

# A tibble containing a breakdown of the missingingness of the census characteristics used to compute the NDI (Powell-Wiley) values
powell_wiley2020DC$missing

# -------------------------------------------- #
# Visualize the powell_wiley() function output #
# -------------------------------------------- #

# Obtain the 2020 census tracts from the "tigris" package
tract2020DC <- tigris::tracts(state = "DC", year = 2020, cb = TRUE)

# Join the NDI (powell_wiley) values to the census tract geometry
DC2020powell_wiley <- merge(tract2020DC, powell_wiley2020DC$ndi, by = "GEOID")
DC2020powell_wiley <- merge(DC2020powell_wiley, powell_wiley2020DCi$ndi, by = "GEOID")

# Visualize the NDI (Powell-Wiley) values (2016-2020 5-year ACS) for Washington, D.C. census tracts

## Non-imputed missing tracts (Continuous)
ggplot2::ggplot() + 
  ggplot2::geom_sf(data = DC2020powell_wiley, 
                   ggplot2::aes(fill = NDI.x),
                   color = "white") +
  ggplot2::theme_minimal() + 
  ggplot2::scale_fill_viridis_c() +
  ggplot2::labs(fill = "Index (Continuous)",
                caption = "Source: U.S. Census ACS 2016-2020 estimates")+
  ggplot2::ggtitle("Neighborhood Deprivation Index\nContinuous (Powell-Wiley, non-imputed)", subtitle = "Washington, D.C. tracts as the referent")

## Non-imputed missing tracts (Categorical quintiles)
### Rename "9-NDI not avail" level as NA for plotting
DC2020powell_wiley$NDIQuintNA.x <- factor(replace(as.character(DC2020powell_wiley$NDIQuint.x),
                                                  DC2020powell_wiley$NDIQuint.x == "9-NDI not avail",
                                                  NA),
                                          c(levels(DC2020powell_wiley$NDIQuint.x)[-6], NA))

ggplot2::ggplot() + 
  ggplot2::geom_sf(data = DC2020powell_wiley, 
                   ggplot2::aes(fill = NDIQuintNA.x),
                   color = "white") +
  ggplot2::theme_minimal() + 
  ggplot2::scale_fill_viridis_d(guide = ggplot2::guide_legend(reverse = TRUE),
                                na.value = "grey50") +
  ggplot2::labs(fill = "Index (Categorical)",
                caption = "Source: U.S. Census ACS 2016-2020 estimates")+
  ggplot2::ggtitle("Neighborhood Deprivation Index\nPopulation-weighted Quintiles (Powell-Wiley, non-imputed)",
                   subtitle = "Washington, D.C. tracts as the referent")
```

![](man/figures/powell_wiley1.png)
![](man/figures/powell_wiley2.png)

``` r
## Imputed missing tracts (Continuous)
ggplot2::ggplot() + 
  ggplot2::geom_sf(data = DC2020powell_wiley, 
                   ggplot2::aes(fill = NDI.y),
                   color = "white") +
  ggplot2::theme_minimal() + 
  ggplot2::scale_fill_viridis_c() +
  ggplot2::labs(fill = "Index (Continuous)",
                caption = "Source: U.S. Census ACS 2016-2020 estimates")+
  ggplot2::ggtitle("Neighborhood Deprivation Index\nContinuous (Powell-Wiley, imputed)", subtitle = "Washington, D.C. tracts as the referent")

## Imputed missing tracts (Categorical quintiles)
### Rename "9-NDI not avail" level as NA for plotting
DC2020powell_wiley$NDIQuintNA.y <- factor(replace(as.character(DC2020powell_wiley$NDIQuint.y), 
                                                  DC2020powell_wiley$NDIQuint.y == "9-NDI not avail",
                                                  NA), 
                                          c(levels(DC2020powell_wiley$NDIQuint.y)[-6], NA))

ggplot2::ggplot() + 
  ggplot2::geom_sf(data = DC2020powell_wiley, 
                   ggplot2::aes(fill = NDIQuintNA.y),
                   color = "white") +
  ggplot2::theme_minimal() + 
  ggplot2::scale_fill_viridis_d(guide = ggplot2::guide_legend(reverse = TRUE),
                                na.value = "grey50") +
  ggplot2::labs(fill = "Index (Categorical)",
                caption = "Source: U.S. Census ACS 2016-2020 estimates")+
  ggplot2::ggtitle("Neighborhood Deprivation Index\nPopulation-weighted Quintiles (Powell-Wiley, imputed)",
                   subtitle = "Washington, D.C. tracts as the referent")
```

![](man/figures/powell_wiley3.png)
![](man/figures/powell_wiley4.png)

``` r
# --------------------------- #
# Compare the two NDI metrics #
# --------------------------- #

# Merge the two NDI metrics (Messer and Powell-Wiley, imputed)
ndi2020DC <- merge(messer2020DC$ndi, powell_wiley2020DCi$ndi, by = "GEOID", suffixes = c(".messer", ".powell_wiley"))

# Check the correlation the two NDI metrics (Messer and Powell-Wiley, imputed) as continuous values
cor(ndi2020DC$NDI.messer, ndi2020DC$NDI.powell_wiley, use = "complete.obs") # Pearsons r = 0.975

# Check the similarity of the two NDI metrics (Messer and Powell-Wiley, imputed) as quartiles
table(ndi2020DC$NDIQuart, ndi2020DC$NDIQuint)
```

### Funding

Package was developed while the author was a postdoctoral fellow supported by the [Cancer Prevention Fellowship Program](https://cpfp.cancer.gov/) at the [National Cancer Institute](https://www.cancer.gov/).

### Acknowledgments

The `messer()` function functionalizes the code found in [Hruska et al. (2022)](https://doi.org/10.1016/j.janxdis.2022.102529) available on an [OSF repository](https://doi.org/10.17605/OSF.IO/M2SAV), but with percent with income less than $30K added to the computation based on [Messer et al. (2006)](https://doi.org/10.1007/s11524-006-9094-x). The `messer()` function also allows for the computation of NDI (Messer) of each year between 2010-2020 (when the U.S. census characteristics are available to-date). There was no code companion to compute NDI (Powell-Wiley) included in [Andrews et al. (2020)](https://doi.org/10.1080/17445647.2020.1750066) or [Slotman et al. (2022)](https://doi.org/10.1016/j.dib.2022.108002), but the package maintainer worked directly with the latter manuscript authors to replicate their `SAS` code in `R` for the `powell-wiley()` function. When citing this package for publication, please follow:

    citation("ndi")

### Questions? Feedback?

For questions about the package please contact the maintainer [Dr. Ian D. Buller](mailto:ian.buller@nih.gov) or [submit a new issue](https://github.com/idblr/ndi/issues). Confirmation of the computation, feedback, and feature collaboration is welcomed, especially from the authors of references cited above.
