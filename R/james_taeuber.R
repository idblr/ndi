#' Dissimilarity Index based on James & Taeuber (1985)
#'
#' Compute the aspatial Dissimilarity Index (James & Taeuber) of selected racial or ethnic subgroup(s) and U.S. geographies
#'
#' @param geo_large Character string specifying the larger geographical unit of the data. The default is counties \code{geo_large = 'county'}.
#' @param geo_small Character string specifying the smaller geographical unit of the data. The default is census tracts \code{geo_large = 'tract'}.
#' @param year Numeric. The year to compute the estimate. The default is 2020, and the years 2009 onward are currently available.
#' @param subgroup Character string specifying the racial or ethnic subgroup(s) as the comparison population. See Details for available choices.
#' @param omit_NAs Logical. If FALSE, will compute index for a larger geographical unit only if all of its smaller geographical units have values. The default is TRUE.
#' @param quiet Logical. If TRUE, will display messages about potential missing census information. The default is FALSE.
#' @param ... Arguments passed to \code{\link[tidycensus]{get_acs}} to select state, county, and other arguments for census characteristics
#'
#' @details This function will compute the aspatial Dissimilarity Index (\emph{D}) of selected racial or ethnic subgroups and U.S. geographies for a specified geographical extent (e.g., the entire U.S. or a single state) based on James & Taeuber (1985) \doi{10.2307/270845}. This function provides the computation of \emph{D} for any of the U.S. Census Bureau race or ethnicity subgroups (including Hispanic and non-Hispanic individuals).
#'
#' The function uses the \code{\link[tidycensus]{get_acs}} function to obtain U.S. Census Bureau 5-year American Community Survey characteristics used for the aspatial computation. The yearly estimates are available for 2009 onward when ACS-5 data are available (2010 onward for \code{geo_large = 'cbsa'} and 2011 onward for \code{geo_large = 'csa'} or \code{geo_large = 'metro'}) but may be available from other U.S. Census Bureau surveys. The twenty racial or ethnic subgroups (U.S. Census Bureau definitions) are:
#' \itemize{
#'  \item \strong{B03002_002}: not Hispanic or Latino \code{'NHoL'}
#'  \item \strong{B03002_003}: not Hispanic or Latino, white alone \code{'NHoLW'}
#'  \item \strong{B03002_004}: not Hispanic or Latino, Black or African American alone \code{'NHoLB'}
#'  \item \strong{B03002_005}: not Hispanic or Latino, American Indian and Alaska Native alone \code{'NHoLAIAN'}
#'  \item \strong{B03002_006}: not Hispanic or Latino, Asian alone \code{'NHoLA'}
#'  \item \strong{B03002_007}: not Hispanic or Latino, Native Hawaiian and Other Pacific Islander alone \code{'NHoLNHOPI'}
#'  \item \strong{B03002_008}: not Hispanic or Latino, Some other race alone \code{'NHoLSOR'}
#'  \item \strong{B03002_009}: not Hispanic or Latino, Two or more races \code{'NHoLTOMR'}
#'  \item \strong{B03002_010}: not Hispanic or Latino, Two races including Some other race \code{'NHoLTRiSOR'}
#'  \item \strong{B03002_011}: not Hispanic or Latino, Two races excluding Some other race, and three or more races \code{'NHoLTReSOR'}
#'  \item \strong{B03002_012}: Hispanic or Latino \code{'HoL'}
#'  \item \strong{B03002_013}: Hispanic or Latino, white alone \code{'HoLW'}
#'  \item \strong{B03002_014}: Hispanic or Latino, Black or African American alone \code{'HoLB'}
#'  \item \strong{B03002_015}: Hispanic or Latino, American Indian and Alaska Native alone \code{'HoLAIAN'}
#'  \item \strong{B03002_016}: Hispanic or Latino, Asian alone \code{'HoLA'}
#'  \item \strong{B03002_017}: Hispanic or Latino, Native Hawaiian and Other Pacific Islander alone \code{'HoLNHOPI'}
#'  \item \strong{B03002_018}: Hispanic or Latino, Some other race alone \code{'HoLSOR'}
#'  \item \strong{B03002_019}: Hispanic or Latino, Two or more races \code{'HoLTOMR'}
#'  \item \strong{B03002_020}: Hispanic or Latino, Two races including Some other race \code{'HoLTRiSOR'}
#'  \item \strong{B03002_021}: Hispanic or Latino, Two races excluding Some other race, and three or more races \code{'HoLTReSOR'}
#' }
#'
#' Use the internal \code{state} and \code{county} arguments within the \code{\link[tidycensus]{get_acs}} function to specify geographic extent of the data output.
#'
#' \emph{D} is a measure of the evenness of racial or ethnic residential segregation when comparing smaller geographical areas to larger ones within which the smaller geographical areas are located. \emph{D} can range in value from 0 to 1 and represents the proportion of racial or ethnic subgroup members that would have to change their area of residence to achieve an even distribution within the larger geographical area under conditions of maximum segregation.
#'
#' Larger geographies available include state \code{geo_large = 'state'}, county \code{geo_large = 'county'}, census tract \code{geo_large = 'tract'}, Core Based Statistical Area \code{geo_large = 'cbsa'}, Combined Statistical Area \code{geo_large = 'csa'}, and Metropolitan Division \code{geo_large = 'metro'} levels. Smaller geographies available include, county \code{geo_small = 'county'}, census tract \code{geo_small = 'tract'}, and census block group \code{geo_small = 'block group'} levels. If a larger geographical area is comprised of only one smaller geographical area (e.g., a U.S county contains only one census tract), then the \emph{D} value returned is NA. If the larger geographical unit is Combined Based Statistical Areas \code{geo_large = 'csa'} or Core Based Statistical Areas \code{geo_large = 'cbsa'}, only the smaller geographical units completely within a larger geographical unit are considered in the \emph{D} computation (see internal \code{\link[sf]{st_within}} function for more information) and recommend specifying all states within which the interested larger geographical unit are located using the internal \code{state} argument to ensure all appropriate smaller geographical units are included in the \emph{D} computation.
#' 
#' @return An object of class 'list'. This is a named list with the following components:
#'
#' \describe{
#' \item{\code{di}}{An object of class 'tbl' for the GEOID, name, and \emph{D} at specified larger census geographies.}
#' \item{\code{d_data}}{An object of class 'tbl' for the raw census values at specified smaller census geographies.}
#' \item{\code{missing}}{An object of class 'tbl' of the count and proportion of missingness for each census variable used to compute \emph{D}.}
#' }
#'
#' @import dplyr
#' @importFrom sf st_drop_geometry st_within
#' @importFrom stats complete.cases
#' @importFrom stringr str_trim
#' @importFrom tidycensus get_acs
#' @importFrom tidyr pivot_longer separate
#' @importFrom tigris combined_statistical_areas core_based_statistical_areas metro_divisions
#' @importFrom utils stack
#' @export
#'
#' @seealso \code{\link[tidycensus]{get_acs}} for additional arguments for geographic extent selection (i.e., \code{state} and \code{county}).
#'
#' @examples
#' \dontrun{
#' # Wrapped in \dontrun{} because these examples require a Census API key.
#'
#'   # Dissimilarity Index (James & Taeuber) 
#'   ## of Black populations
#'   ## of census tracts within counties within Georgia, U.S.A., counties (2020)
#'   james_taeuber(
#'     geo_large = 'county',
#'     geo_small = 'tract',
#'     state = 'GA',
#'     year = 2020,
#'     subgroup = c('NHoLB', 'HoLB')
#'    )
#'
#' }
#'
james_taeuber <- function(geo_large = 'county',
                          geo_small = 'tract',
                          year = 2020,
                          subgroup,
                          omit_NAs = TRUE,
                          quiet = FALSE,
                          ...) {
  
  # Check arguments
  match.arg(geo_large, choices = c('state', 'county', 'tract', 'cbsa', 'csa', 'metro'))
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
  
  # Select census variable
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
  in_subgroup <- paste0(subgroup, 'E')
  
  # Acquire D variables and sf geometries
  out_dat <- suppressMessages(suppressWarnings(
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
    out_dat <- out_dat %>%
      tidyr::separate(NAME.y, into = c('county', 'state'), sep = ',')
  }
  if (geo_small == 'tract') {
    out_dat <- out_dat %>%
      tidyr::separate(NAME.y, into = c('tract', 'county', 'state'), sep = ',') %>%
      dplyr::mutate(tract = gsub('[^0-9\\.]', '', tract))
  }
  if (geo_small == 'block group') {
    out_dat <- out_dat %>%
      tidyr::separate(NAME.y, into = c('block.group', 'tract', 'county', 'state'), sep = ',') %>%
      dplyr::mutate(
        tract = gsub('[^0-9\\.]', '', tract),
        block.group = gsub('[^0-9\\.]', '', block.group)
      )
  }
  
  # Grouping IDs for D computation
  if (geo_large == 'state') {
    out_dat <- out_dat %>%
      dplyr::mutate(
        oid = STATEFP,
        state = stringr::str_trim(state)
      ) %>%
      sf::st_drop_geometry()
  }
  if (geo_large == 'tract') {
    out_dat <- out_dat %>%
      dplyr::mutate(
        oid = paste0(STATEFP, COUNTYFP, TRACTCE),
        state = stringr::str_trim(state),
        county = stringr::str_trim(county)
      ) %>%
      sf::st_drop_geometry()
  }
  if (geo_large == 'county') {
    out_dat <- out_dat %>%
      dplyr::mutate(
        oid = paste0(STATEFP, COUNTYFP),
        state = stringr::str_trim(state),
        county = stringr::str_trim(county)
      ) %>%
      sf::st_drop_geometry()
  }
  if (geo_large == 'cbsa') {
    stopifnot(is.numeric(year), year >= 2010) # CBSAs only available 2010 onward
    lgeom <- suppressMessages(suppressWarnings(tigris::core_based_statistical_areas(year = year)))
    wlgeom <- sf::st_within(out_dat, lgeom)
    out_dat <- out_dat %>%
      dplyr::mutate(
        oid = lapply(wlgeom, function(x) { 
          tmp <- lgeom[x, 3] %>% sf::st_drop_geometry()
          lapply(tmp, function(x) { if (length(x) == 0) NA else x })
        }) %>% 
          unlist(),
        cbsa = lapply(wlgeom, function(x) { 
          tmp <- lgeom[x, 4] %>% sf::st_drop_geometry()
          lapply(tmp, function(x) { if (length(x) == 0) NA else x })
        }) %>% 
          unlist()
      ) %>% 
      sf::st_drop_geometry()
  }
  if (geo_large == 'csa') {
    stopifnot(is.numeric(year), year >= 2011) # CSAs only available 2011 onward
    lgeom <- suppressMessages(suppressWarnings(tigris::combined_statistical_areas(year = year)))
    wlgeom <- sf::st_within(out_dat, lgeom)
    out_dat <- out_dat %>%
      dplyr::mutate(
        oid = lapply(wlgeom, function(x) { 
          tmp <- lgeom[x, 2] %>% sf::st_drop_geometry()
          lapply(tmp, function(x) { if (length(x) == 0) NA else x })
        }) %>% 
          unlist(),
        csa = lapply(wlgeom, function(x) { 
          tmp <- lgeom[x, 3] %>% sf::st_drop_geometry()
          lapply(tmp, function(x) { if (length(x) == 0) NA else x })
        }) %>% 
          unlist()
      ) %>% 
      sf::st_drop_geometry()
  }
  if (geo_large == 'metro') {
    stopifnot(is.numeric(year), year >= 2011) # Metro Divisions only available 2011 onward
    lgeom <- suppressMessages(suppressWarnings(tigris::metro_divisions(year = year)))
    wlgeom <- sf::st_within(out_dat, lgeom)
    out_dat <- out_dat %>%
      dplyr::mutate(
        oid = lapply(wlgeom, function(x) { 
          tmp <- lgeom[x, 4] %>% sf::st_drop_geometry()
          lapply(tmp, function(x) { if (length(x) == 0) NA else x })
        }) %>% 
          unlist(),
        metro = lapply(wlgeom, function(x) { 
          tmp <- lgeom[x, 5] %>% sf::st_drop_geometry()
          lapply(tmp, function(x) { if (length(x) == 0) NA else x })
        }) %>% 
          unlist()
      ) %>% 
      sf::st_drop_geometry()
  }
  
  # Count of racial or ethnic subgroup populations
  ## Count of racial or ethnic comparison subgroup population
  if (length(in_subgroup) == 1) {
    out_dat <- out_dat %>%
      dplyr::mutate(subgroup = .[, in_subgroup])
  } else {
    out_dat <- out_dat %>%
      dplyr::mutate(subgroup = rowSums(.[, in_subgroup]))
  }

  # Compute D
  ## From James & Taeuber (1985) https://doi.org/10.2307/270845
  ## D = \frac{\sum_{n}^{i=1}t_{i}\left|p_{i}-P\right|}{2TP(1-P)}
  ## Where for i smaller geographies:
  ## t_{i} is the total population of area i 
  ## p_{i} is the proportion of the subgroup population of area i
  ## T is the total population of all areas
  ## P is the proportion of the subgroup population of all areas
  
  ## Compute
  out_tmp <- out_dat %>%
    split(., f = list(out_dat$oid)) %>%
    lapply(., FUN = djt_fun, omit_NAs = omit_NAs) %>%
    utils::stack(.) %>%
    dplyr::mutate(
      D = values,
      oid = ind
    ) %>%
    dplyr::select(D, oid)
  
  # Warning for missingness of census characteristics
  missingYN <- out_dat[, c('TotalPopE', in_subgroup)]
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
  if (geo_large == 'state') {
    out <- out_dat %>%
      dplyr::left_join(out_tmp, by = dplyr::join_by(oid)) %>%
      dplyr::select(oid, state, D) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, state, D) %>%
      .[.$GEOID != 'NANA',]
  }
  if (geo_large == 'county') {
    out <- out_dat %>%
      dplyr::left_join(out_tmp, by = dplyr::join_by(oid)) %>%
      dplyr::select(oid, state, county, D) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, state, county, D) %>%
      .[.$GEOID != 'NANA',]
  }
  if (geo_large == 'tract') {
    out <- out_dat %>%
      dplyr::left_join(out_tmp, by = dplyr::join_by(oid)) %>%
      dplyr::select(oid, state, county, tract, D) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, state, county, tract, D) %>%
      .[.$GEOID != 'NANA',]
  }
  if (geo_large == 'cbsa') {
    out <- out_dat %>%
      dplyr::left_join(out_tmp, by = dplyr::join_by(oid)) %>%
      dplyr::select(oid, cbsa, D) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, cbsa, D) %>%
      .[.$GEOID != 'NANA', ] %>%
      dplyr::distinct(GEOID, .keep_all = TRUE) %>%
      dplyr::filter(stats::complete.cases(.))
  }
  if (geo_large == 'csa') {
    out <- out_dat %>%
      dplyr::left_join(out_tmp, by = dplyr::join_by(oid)) %>%
      dplyr::select(oid, csa, D) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, csa, D) %>%
      .[.$GEOID != 'NANA', ] %>%
      dplyr::distinct(GEOID, .keep_all = TRUE) %>%
      dplyr::filter(stats::complete.cases(.))
  }
  if (geo_large == 'metro') {
    out <- out_dat %>%
      dplyr::left_join(out_tmp, by = dplyr::join_by(oid)) %>%
      dplyr::select(oid, metro, D) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, metro, D) %>%
      .[.$GEOID != 'NANA', ] %>%
      dplyr::distinct(GEOID, .keep_all = TRUE) %>%
      dplyr::filter(stats::complete.cases(.))
  }
  
  out <- out %>%
    dplyr::arrange(GEOID) %>%
    dplyr::as_tibble()
  
  out_dat <- out_dat %>%
    dplyr::arrange(GEOID) %>%
    dplyr::as_tibble()
  
  out <- list(d = out, d_data = out_dat, missing = missingYN)
  
  return(out)
}
