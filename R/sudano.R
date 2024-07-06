#' Location Quotient (LQ) based on Merton (1938) and Sudano et al. (2013)
#'
#' Compute the aspatial Location Quotient (Sudano) of a selected racial/ethnic subgroup(s) and U.S. geographies.
#'
#' @param geo_large Character string specifying the larger geographical unit of the data. The default is counties \code{geo_large = 'county'}.
#' @param geo_small Character string specifying the smaller geographical unit of the data. The default is census tracts \code{geo_large = 'tract'}.
#' @param year Numeric. The year to compute the estimate. The default is 2020, and the years 2009 onward are currently available.
#' @param subgroup Character string specifying the racial/ethnic subgroup(s). See Details for available choices.
#' @param omit_NAs Logical. If FALSE, will compute index for a larger geographical unit only if all of its smaller geographical units have values. The default is TRUE.
#' @param quiet Logical. If TRUE, will display messages about potential missing census information. The default is FALSE.
#' @param ... Arguments passed to \code{\link[tidycensus]{get_acs}} to select state, county, and other arguments for census characteristics
#'
#' @details This function will compute the aspatial Location Quotient (LQ) of selected racial/ethnic subgroups and U.S. geographies for a specified geographical extent (e.g., the entire U.S. or a single state) based on Merton (1939) \doi{10.2307/2084686} and Sudano et al. (2013) \doi{10.1016/j.healthplace.2012.09.015}. This function provides the computation of LQ for any of the U.S. Census Bureau race/ethnicity subgroups (including Hispanic and non-Hispanic individuals).
#'
#' The function uses the \code{\link[tidycensus]{get_acs}} function to obtain U.S. Census Bureau 5-year American Community Survey characteristics used for the aspatial computation. The yearly estimates are available for 2009 onward when ACS-5 data are available but are available from other U.S. Census Bureau surveys. The twenty racial/ethnic subgroups (U.S. Census Bureau definitions) are:
#' \itemize{
#'  \item **B03002_002**: not Hispanic or Latino \code{'NHoL'}
#'  \item **B03002_003**: not Hispanic or Latino, white alone \code{'NHoLW'}
#'  \item **B03002_004**: not Hispanic or Latino, Black or African American alone \code{'NHoLB'}
#'  \item **B03002_005**: not Hispanic or Latino, American Indian and Alaska Native alone \code{'NHoLAIAN'}
#'  \item **B03002_006**: not Hispanic or Latino, Asian alone \code{'NHoLA'}
#'  \item **B03002_007**: not Hispanic or Latino, Native Hawaiian and Other Pacific Islander alone \code{'NHoLNHOPI'}
#'  \item **B03002_008**: not Hispanic or Latino, Some other race alone \code{'NHoLSOR'}
#'  \item **B03002_009**: not Hispanic or Latino, Two or more races \code{'NHoLTOMR'}
#'  \item **B03002_010**: not Hispanic or Latino, Two races including Some other race \code{'NHoLTRiSOR'}
#'  \item **B03002_011**: not Hispanic or Latino, Two races excluding Some other race, and three or more races \code{'NHoLTReSOR'}
#'  \item **B03002_012**: Hispanic or Latino \code{'HoL'}
#'  \item **B03002_013**: Hispanic or Latino, white alone \code{'HoLW'}
#'  \item **B03002_014**: Hispanic or Latino, Black or African American alone \code{'HoLB'}
#'  \item **B03002_015**: Hispanic or Latino, American Indian and Alaska Native alone \code{'HoLAIAN'}
#'  \item **B03002_016**: Hispanic or Latino, Asian alone \code{'HoLA'}
#'  \item **B03002_017**: Hispanic or Latino, Native Hawaiian and Other Pacific Islander alone \code{'HoLNHOPI'}
#'  \item **B03002_018**: Hispanic or Latino, Some other race alone \code{'HoLSOR'}
#'  \item **B03002_019**: Hispanic or Latino, Two or more races \code{'HoLTOMR'}
#'  \item **B03002_020**: Hispanic or Latino, Two races including Some other race \code{'HoLTRiSOR'}
#'  \item **B03002_021**: Hispanic or Latino, Two races excluding Some other race, and three or more races \code{'HoLTReSOR'}
#' }
#'
#' Use the internal \code{state} and \code{county} arguments within the \code{\link[tidycensus]{get_acs}} function to specify geographic extent of the data output.
#'
#' LQ is some measure of relative racial homogeneity of each smaller geography within a larger geography. LQ can range in value from 0 to infinity because it is ratio of two proportions in which the numerator is the proportion of subgroup population in a smaller geography and the denominator is the proportion of subgroup population in its larger geography. For example, a smaller geography with an LQ of 5 means that the proportion of the subgroup population living in the smaller geography is five times the proportion of the subgroup population in its larger geography.
#'
#' Larger geographies available include state \code{geo_large = 'state'}, county \code{geo_large = 'county'}, and census tract \code{geo_large = 'tract'} levels. Smaller geographies available include, county \code{geo_small = 'county'}, census tract \code{geo_small = 'tract'}, and census block group \code{geo_small = 'block group'} levels. If a larger geographical area is comprised of only one smaller geographical area (e.g., a U.S county contains only one census tract), then the LQ value returned is NA.
#'
#' @return An object of class 'list'. This is a named list with the following components:
#'
#' \describe{
#' \item{\code{lq}}{An object of class 'tbl' for the GEOID, name, and LQ at specified smaller census geographies.}
#' \item{\code{lq_data}}{An object of class 'tbl' for the raw census values at specified smaller census geographies.}
#' \item{\code{missing}}{An object of class 'tbl' of the count and proportion of missingness for each census variable used to compute LQ.}
#' }
#'
#' @import dplyr
#' @importFrom sf st_drop_geometry
#' @importFrom stats complete.cases
#' @importFrom tidycensus get_acs
#' @importFrom tidyr pivot_longer separate
#' @importFrom utils stack
#' @export
#'
#' @seealso \code{\link[tidycensus]{get_acs}} for additional arguments for geographic extent selection (i.e., \code{state} and \code{county}).
#'
#' @examples
#' \dontrun{
#' # Wrapped in \dontrun{} because these examples require a Census API key.
#'
#'   # Isolation of non-Hispanic Black populations
#'   ## of census tracts within Georgia, U.S.A., counties (2020)
#'   sudano(
#'     geo_large = 'state',
#'     geo_small = 'county',
#'     state = 'GA',
#'     year = 2020,
#'     subgroup = 'NHoLB'
#'    )
#'
#' }
#'
sudano <- function(geo_large = 'county',
                   geo_small = 'tract',
                   year = 2020,
                   subgroup,
                   omit_NAs = TRUE,
                   quiet = FALSE,
                   ...) {
  
  # Check arguments
  match.arg(geo_large, choices = c('state', 'county', 'tract'))
  match.arg(geo_small, choices = c('county', 'tract', 'block group'))
  stopifnot(is.numeric(year), year >= 2009) # all variables available 2009 onward
  match.arg(
    subgroup,
    several.ok = TRUE,
    choices = c(
      'NHoL',
      'NHoLW',
      'NHoLB',
      'NHoLAIAN',
      'NHoLA',
      'NHoLNHOPI',
      'NHoLSOR',
      'NHoLTOMR',
      'NHoLTRiSOR',
      'NHoLTReSOR',
      'HoL',
      'HoLW',
      'HoLB',
      'HoLAIAN',
      'HoLA',
      'HoLNHOPI',
      'HoLSOR',
      'HoLTOMR',
      'HoLTRiSOR',
      'HoLTReSOR'
    )
  )
  
  # Select census variables
  vars <- c(
    TotalPop = 'B03002_001',
    NHoL = 'B03002_002',
    NHoLW = 'B03002_003',
    NHoLB = 'B03002_004',
    NHoLAIAN = 'B03002_005',
    NHoLA = 'B03002_006',
    NHoLNHOPI = 'B03002_007',
    NHoLSOR = 'B03002_008',
    NHoLTOMR = 'B03002_009',
    NHoLTRiSOR = 'B03002_010',
    NHoLTReSOR = 'B03002_011',
    HoL = 'B03002_012',
    HoLW = 'B03002_013',
    HoLB = 'B03002_014',
    HoLAIAN = 'B03002_015',
    HoLA = 'B03002_016',
    HoLNHOPI = 'B03002_017',
    HoLSOR = 'B03002_018',
    HoLTOMR = 'B03002_019',
    HoLTRiSOR = 'B03002_020',
    HoLTReSOR = 'B03002_021'
  )
  
  selected_vars <- vars[c('TotalPop', subgroup)]
  out_names <- names(selected_vars) # save for output
  in_subgroup <- paste(subgroup, 'E', sep = '')
  
  # Acquire LQ variables and sf geometries
  lq_data <- suppressMessages(suppressWarnings(
    tidycensus::get_acs(
      geography = geo_small,
      year = year,
      output = 'wide',
      variables = selected_vars,
      geometry = TRUE,
      keep_geo_vars = TRUE,
      ...
    )
  ))
    
  
  # Format output
  if (geo_small == 'county') {
    lq_data <- lq_data %>%
      sf::st_drop_geometry() %>%
      tidyr::separate(NAME.y, into = c('county', 'state'), sep = ',')
  }
  if (geo_small == 'tract') {
    lq_data <- lq_data %>%
      sf::st_drop_geometry() %>%
      tidyr::separate(NAME.y, into = c('tract', 'county', 'state'), sep = ',') %>%
      dplyr::mutate(tract = gsub('[^0-9\\.]', '', tract))
  }
  if (geo_small == 'block group') {
    lq_data <- lq_data %>%
      sf::st_drop_geometry() %>%
      tidyr::separate(NAME.y, into = c('block.group', 'tract', 'county', 'state'), sep = ',') %>%
      dplyr::mutate(
        tract = gsub('[^0-9\\.]', '', tract),
        block.group = gsub('[^0-9\\.]', '', block.group)
      )
  }
  
  # Grouping IDs for R computation
  if (geo_large == 'tract') {
    lq_data <- lq_data %>%
      dplyr::mutate(
        oid = paste(.$STATEFP, .$COUNTYFP, .$TRACTCE, sep = ''),
        state = stringr::str_trim(state),
        county = stringr::str_trim(county)
      )
  }
  if (geo_large == 'county') {
    lq_data <- lq_data %>%
      dplyr::mutate(
        oid = paste(.$STATEFP, .$COUNTYFP, sep = ''),
        state = stringr::str_trim(state),
        county = stringr::str_trim(county)
      )
  }
  if (geo_large == 'state') {
    lq_data <- lq_data %>%
      dplyr::mutate(
        oid = .$STATEFP,
        state = stringr::str_trim(state)
      )
  }
  
  # Count of racial/ethnic subgroup populations
  ## Count of racial/ethnic comparison subgroup population
  if (length(in_subgroup) == 1) {
    lq_data <- lq_data %>%
      dplyr::mutate(subgroup = .[, in_subgroup])
  } else {
    lq_data <- lq_data %>%
      dplyr::mutate(subgroup = rowSums(.[, in_subgroup]))
  }
  
  # Compute LQ
  ## From Sudano (2013) https://doi.org/10.1016/j.healthplace.2012.09.015
  ## LQ_{im} = (x_{im}/X_{i})/(X_{m}/X)
  ## for:
  ## i smaller geography and subgroup m
  
  ## Compute
  LQtmp <- lq_data %>%
    split(., f = list(lq_data$oid)) %>%
    lapply(., FUN = lq_fun, omit_NAs = omit_NAs) %>%
    do.call('rbind', .)
  
  # Warning for missingness of census characteristics
  missingYN <- lq_data[, c('TotalPopE', in_subgroup)]
  names(missingYN) <- out_names
  missingYN <- missingYN %>%
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
  
  # Format output
  lq <- lq_data %>%
    dplyr::left_join(LQtmp, by = dplyr::join_by(GEOID))
  
  if (geo_small == 'state') {
    lq <- lq %>%
      dplyr::select(GEOID, state, LQ)
  }
  if (geo_small == 'county') {
    lq <- lq %>%
      dplyr::select(GEOID, state, county, LQ)
  }
  if (geo_small == 'tract') {
    lq <- lq %>%
      dplyr::select(GEOID, state, county, tract, LQ)
  }
  if (geo_small == 'block group') {
    lq <- lq %>%
      dplyr::select(GEOID, state, county, tract, block.group, LQ)
  }
  
  lq <- lq %>%
    unique(.) %>%
    .[.$GEOID != 'NANA',] %>%
    dplyr::arrange(GEOID) %>%
    dplyr::as_tibble()
  
  lq_data <- lq_data %>%
    dplyr::arrange(GEOID) %>%
    dplyr::as_tibble()
  
  out <- list(lq = lq, lq_data = lq_data, missing = missingYN)
  
  return(out)
}
