#' Gini Index based on Gini (1921)
#'
#' Retrieve the aspatial Gini Index of income inequality.
#'
#' @param geo Character string specifying the geography of the data either census tracts \code{geo = 'tract'} (the default) or counties \code{geo = 'county'}.
#' @param year Numeric. The year to compute the estimate. The default is 2020, and the years 2009 onward are currently available.
#' @param quiet Logical. If TRUE, will display messages about potential missing census information
#' @param ... Arguments passed to \code{\link[tidycensus]{get_acs}} to select state, county, and other arguments for census characteristics
#'
#' @details This function will retrieve the aspatial Gini Index (\emph{G}) of U.S. census tracts or counties for a specified geographical extent (e.g., the entire U.S. or a single state) based on Gini (1921) \doi{10.2307/2223319}.
#'
#' The function uses the \code{\link[tidycensus]{get_acs}} function to obtain U.S. Census Bureau 5-year American Community Survey estimates of \emph{G} for income inequality (ACS: B19083). The estimates are available for 2009 onward when ACS-5 data are available but are available from other U.S. Census Bureau surveys.
#'
#' Use the internal \code{state} and \code{county} arguments within the \code{\link[tidycensus]{get_acs}} function to specify geographic extent of the data output.
#'
#' According to the U.S. Census Bureau \url{https://www.census.gov/topics/income-poverty/income-inequality/about/metrics/gini-index.html}: 'The Gini Index is a summary measure of income inequality. The Gini coefficient incorporates the detailed shares data into a single statistic, which summarizes the dispersion of income across the entire income distribution. The Gini coefficient ranges from 0, indicating perfect equality (where everyone receives an equal share), to 1, perfect inequality (where only one recipient or group of recipients receives all the income). The Gini is based on the difference between the Lorenz curve (the observed cumulative income distribution) and the notion of a perfectly equal income distribution.'
#'
#' @return An object of class 'list'. This is a named list with the following components:
#'
#' \describe{
#' \item{\code{gini}}{An object of class 'tbl' for the GEOID, name, and \emph{G} of specified census geographies.}
#' \item{\code{missing}}{An object of class 'tbl' of the count and proportion of missingness for \emph{G}.}
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
#'   # Gini Index of income inequality
#'   ## of census tracts within Georgia, U.S.A., counties (2020)
#'   gini(geo = 'tract', state = 'GA', year = 2020)
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
  tmp_dat <- suppressMessages(suppressWarnings(
    tidycensus::get_acs(
      geography = geo,
      year = year,
      output = 'wide',
      variables = vars,
      ...
    )
  ))
    
  if (geo == 'tract') {
    tmp_dat <- tmp_dat %>%
      tidyr::separate(NAME, into = c('tract', 'county', 'state'), sep = ',') %>%
      dplyr::mutate(tract = gsub('[^0-9\\.]', '', tract))
  } else {
    tmp_dat <- tmp_dat %>% 
      tidyr::separate(NAME, into = c('county', 'state'), sep = ',')
  }
  
  tmp_dat <- tmp_dat %>%
    dplyr::mutate(gini = giniE)
  
  # Warning for missingness of census characteristics
  missingYN <- tmp_dat %>%
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
    out <- tmp_dat %>%
      dplyr::select(GEOID,  state, county, tract, gini)
  } else {
    out <- tmp_dat %>%
      dplyr::select(GEOID,  state, county, gini)
  }
  
  out <- out %>%
    dplyr::mutate(
      state = stringr::str_trim(state),
      county = stringr::str_trim(county)
    ) %>%
    dplyr::arrange(GEOID) %>%
    dplyr::as_tibble()
  
  out <- list(gini = out, missing = missingYN)
  
  return(out)
}
