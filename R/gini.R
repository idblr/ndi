#' Gini Index based on Gini (1921)
#'
#' Retrieve the aspatial Gini Index of income inequality.
#'
#' @param geo Character string specifying the geography of the data either census tracts \code{geo = 'tract'} (the default) or counties \code{geo = 'county'}.
#' @param year Numeric. The year to compute the estimate. The default is 2020, and the years 2009 onward are currently available.
#' @param quiet Logical. If TRUE, will display messages about potential missing census information
#' @param ... Arguments passed to \code{\link[tidycensus]{get_acs}} to select state, county, and other arguments for census characteristics
#'
#' @details This function will retrieve the aspatial Gini Index of U.S. census tracts or counties for a specified geographical extent (e.g., the entire U.S. or a single state) based on Gini (1921) \doi{10.2307/2223319}.
#'
#' The function uses the \code{\link[tidycensus]{get_acs}} function to obtain U.S. Census Bureau 5-year American Community Survey estimates of the Gini Index for income inequality (ACS: B19083). The estimates are available for 2009 onward when ACS-5 data are available but are available from other U.S. Census Bureau surveys.
#'
#' Use the internal \code{state} and \code{county} arguments within the \code{\link[tidycensus]{get_acs}} function to specify geographic extent of the data output.
#'
#' According to the U.S. Census Bureau \url{https://www.census.gov/topics/income-poverty/income-inequality/about/metrics/gini-index.html}: 'The Gini Index is a summary measure of income inequality. The Gini coefficient incorporates the detailed shares data into a single statistic, which summarizes the dispersion of income across the entire income distribution. The Gini coefficient ranges from 0, indicating perfect equality (where everyone receives an equal share), to 1, perfect inequality (where only one recipient or group of recipients receives all the income). The Gini is based on the difference between the Lorenz curve (the observed cumulative income distribution) and the notion of a perfectly equal income distribution.'
#'
#' @return An object of class 'list'. This is a named list with the following components:
#'
#' \describe{
#' \item{\code{gini}}{An object of class 'tbl' for the GEOID, name, and Gini index of specified census geographies.}
#' \item{\code{missing}}{An object of class 'tbl' of the count and proportion of missingness for the Gini index.}
#' }
#'
#' @import dplyr
#' @importFrom stringr str_trim
#' @importFrom tidycensus get_acs
#' @importFrom tidyr pivot_longer separate
#' @export
#'
#' @seealso \code{\link[tidycensus]{get_acs}} for additional arguments for geographic extent selection (i.e., \code{state} and \code{county}).
#'
#' @examples
#' \dontrun{
#' # Wrapped in \dontrun{} because these examples require a Census API key.
#'
#'   # Tract-level metric (2020)
#'   gini(geo = 'tract', state = 'GA', year = 2020)
#'
#'   # County-level metric (2020)
#'   gini(geo = 'county', state = 'GA', year = 2020)
#'
#' }
#'
gini <- function(geo = 'tract',
                 year = 2020,
                 quiet = FALSE,
                 ...) {
  
  # Check arguments
  match.arg(geo, choices = c('county', 'tract'))
  stopifnot(is.numeric(year), year >= 2009) # the gini variable is available before and after 2009 but constrained for consistency with out indices (for now)
  
  # Select census variable
  vars <- c(gini = 'B19083_001')
  
  # Acquire Gini Index
  gini_data <- suppressMessages(suppressWarnings(
    tidycensus::get_acs(
      geography = geo,
      year = year,
      output = 'wide',
      variables = vars,
      ...
    )
  ))
    
  if (geo == 'tract') {
    gini_data <- gini_data %>%
      tidyr::separate(NAME, into = c('tract', 'county', 'state'), sep = ',') %>%
      dplyr::mutate(tract = gsub('[^0-9\\.]', '', tract))
  } else {
    gini_data <- gini_data %>% 
      tidyr::separate(NAME, into = c('county', 'state'), sep = ',')
  }
  
  gini_data <- gini_data %>%
    dplyr::mutate(gini = giniE)
  
  # Warning for missingness of census characteristics
  missingYN <- gini_data %>%
    dplyr::select(gini)  %>%
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = 'variable',
      values_to = 'val'
    ) %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(
      total = dplyr::n(),
      n_missing = sum(is.na(val)),
      percent_missing = paste0(round(mean(is.na(val)) * 100, 2), ' %')
    )
  
  if (quiet == FALSE) {
    # Warning for missing census data
    if (sum(missingYN$n_missing) > 0) {
      message('Warning: Missing census data')
    }
  }
  
  if (geo == 'tract') {
    gini <- gini_data %>%
      dplyr::select(GEOID,  state, county, tract, gini)
  } else {
    gini <- gini_data %>%
      dplyr::select(GEOID,  state, county, gini)
  }
  
  gini <- gini %>%
    dplyr::mutate(
      state = stringr::str_trim(state),
      county = stringr::str_trim(county)
    ) %>%
    dplyr::arrange(GEOID) %>%
    dplyr::as_tibble()
  
  out <- list(gini = gini, missing = missingYN)
  
  return(out)
}
