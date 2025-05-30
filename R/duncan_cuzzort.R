#' Absolute Centralization based on Duncan, Cuzzort, & Duncan (1961) and Massey & Denton (1988)
#' 
#' Compute the aspatial Absolute Centralization (Duncan & Cuzzort) of a selected racial or ethnic subgroup(s) and U.S. geographies.
#'
#' @param geo_large Character string specifying the larger geographical unit of the data. The default is counties \code{geo_large = 'county'}.
#' @param geo_small Character string specifying the smaller geographical unit of the data. The default is census tracts \code{geo_small = 'tract'}.
#' @param year Numeric. The year to compute the estimate. The default is 2020, and the years 2009 onward are currently available.
#' @param subgroup Character string specifying the racial or ethnic subgroup(s) as the comparison population. See Details for available choices.
#' @param crs Numeric or character string specifying the coordinate reference system to compute the distance-based metric. The default is Albers North America \code{crs = 'ESRI:102008'}.
#' @param omit_NAs Logical. If FALSE, will compute index for a larger geographical unit only if all of its smaller geographical units have values. The default is TRUE.
#' @param quiet Logical. If TRUE, will display messages about potential missing census information. The default is FALSE.
#' @param ... Arguments passed to \code{\link[tidycensus]{get_acs}} to select state, county, and other arguments for census characteristics
#'
#' @details This function will compute the aspatial Absolute Centralization (\emph{ACE}) of selected racial or ethnic subgroups and U.S. geographies for a specified geographical extent (e.g., the entire U.S. or a single state) based on Duncan, Cuzzort, & Duncan (1961; LC:60007089) and Massey & Denton (1988) \doi{10.1093/sf/67.2.281}. This function provides the computation of \emph{ACE} for any of the U.S. Census Bureau race or ethnicity subgroups (including Hispanic and non-Hispanic individuals).
#' 
#' The function uses the \code{\link[tidycensus]{get_acs}} function to obtain U.S. Census Bureau 5-year American Community Survey characteristics used for the computation. The yearly estimates are available for 2009 onward when ACS-5 data are available (2010 onward for \code{geo_large = 'cbsa'} and 2011 onward for \code{geo_large = 'place'}, \code{geo_large = 'csa'}, or \code{geo_large = 'metro'}) but may be available from other U.S. Census Bureau surveys. The twenty racial or ethnic subgroups (U.S. Census Bureau definitions) are:
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
#' \emph{ACE} is a measure of the degree to which racial or ethnic populations within smaller geographical units are located near the center of a larger geographical unit. \emph{ACE} can range in value from -1 to 1 and represents the spatial distribution of racial or ethnic populations within smaller geographical units compared to the distribution of land area around the center of a larger geographical unit. Positive values indicate a tendency for racial or ethnic populations to reside close to the center of a larger geographical unit, while negative values indicate a tendency to live in outlying areas. A score of 0 means that racial or ethnic populations have a uniform distribution throughout a larger geographical unit. \emph{ACE} gives the proportion of racial or ethnic populations required to change residence to achieve a uniform distribution of population around the center of a larger geographical unit.
#'
#' Larger geographical units available include states \code{geo_large = 'state'}, counties \code{geo_large = 'county'}, census tracts \code{geo_large = 'tract'}, census-designated places \code{geo_large = 'place'}, core-based statistical areas \code{geo_large = 'cbsa'}, combined statistical areas \code{geo_large = 'csa'}, and metropolitan divisions \code{geo_large = 'metro'}. Smaller geographical units available include, counties \code{geo_small = 'county'}, census tracts \code{geo_small = 'tract'}, and census block groups \code{geo_small = 'cbg'}. If a larger geographical unit is comprised of only one smaller geographical unit (e.g., a U.S county contains only one census tract), then the \emph{ACE} value returned is NA. If the larger geographical unit is census-designated places \code{geo_large = 'place'}, core-based statistical areas \code{geo_large = 'cbsa'}, combined statistical areas \code{geo_large = 'csa'}, or metropolitan divisions \code{geo_large = 'metro'}, only the smaller geographical units completely within a larger geographical unit are considered in the \emph{V} computation (see internal \code{\link[sf]{st_within}} function for more information) and recommend specifying all states within which the interested larger geographical unit are located using the internal \code{state} argument to ensure all appropriate smaller geographical units are included in the \emph{ACE} computation.
#' 
#' \emph{Important consideration}: The original metric used the location of the central business district (CBD) to compute the metric, but the U.S. Census Bureau has not defined CBDs for U.S. cities since the 1982 Census of Retail Trade. Therefore, this function uses the the centroids of each larger geographical unit as the 'centre', but may not represent the current CBD. 
#' 
#' @return An object of class 'list'. This is a named list with the following components:
#'
#' \describe{
#' \item{\code{ace}}{An object of class 'tbl' for the GEOID, name, and \emph{ACE} at specified larger census geographies.}
#' \item{\code{ace_data}}{An object of class 'tbl' for the raw census values at specified smaller census geographies.}
#' \item{\code{missing}}{An object of class 'tbl' of the count and proportion of missingness for each census variable used to compute \emph{ACE}.}
#' }
#'
#' @seealso \code{\link[tidycensus]{get_acs}} for additional arguments for geographic extent selection (i.e., \code{state} and \code{county}).
#' @seealso Relative Centralization: \code{\link{duncan_duncan}}
#' 
#' @references Duncan, OD, Cuzzort, RP, & Duncan, B (1961) \emph{Statistical Geography: Problems in Analyzing Area Data}. Free Press. LC:60007089
#' @references Massey, DS, & Denton, NA (1988) The Dimensions of Residential Segregation. \emph{Social Forces}, 67(1):281-315. \doi{10.1093/sf/67.2.281}
#' 
#' @import dplyr
#' @importFrom sf st_centroid st_distance st_drop_geometry st_transform st_within
#' @importFrom stats complete.cases
#' @importFrom stringr str_trim
#' @importFrom tidycensus get_acs
#' @importFrom tidyr pivot_longer separate
#' @importFrom tigris combined_statistical_areas core_based_statistical_areas counties metro_divisions places states tracts
#' @importFrom units drop_units set_units
#' @importFrom utils stack
#' @export
#'
#' @examples
#' \dontrun{
#' # Wrapped in \dontrun{} because these examples require a Census API key.
#'
#'   # Absolute Centralization (a measure of centralization)
#'   ## of Black populations
#'   ## in census tracts within counties of Georgia, U.S.A. (2020)
#'   duncan_cuzzort(
#'     geo_large = 'county',
#'     geo_small = 'tract',
#'     state = 'GA',
#'     year = 2020,
#'     subgroup = c('NHoLB', 'HoLB')
#'    )
#'
#' }
#'
duncan_cuzzort <- function(geo_large = 'county',
                           geo_small = 'tract',
                           year = 2020,
                           subgroup,
                           crs = 'ESRI:102008',
                           omit_NAs = TRUE,
                           quiet = FALSE,
                           ...) {
  
  # Check arguments
  match.arg(geo_large, choices = c('state', 'county', 'tract', 'place', 'cbsa', 'csa', 'metro'))
  match.arg(geo_small, choices = c('county', 'tract', 'cbg', 'block group'))
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
  out_names <- c(names(selected_vars), 'ALAND') # save for output
  in_subgroup <- paste0(subgroup, 'E')
  
  # Acquire ACE variables and sf geometries
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
  if (geo_small == 'cbg' | geo_small == 'block group') {
    out_dat <- out_dat %>%
      tidyr::separate(NAME.y, into = c('cbg', 'tract', 'county', 'state'), sep = ',') %>%
      dplyr::mutate(
        tract = gsub('[^0-9\\.]', '', tract),
        cbg = gsub('[^0-9\\.]', '', cbg)
      )
  }
  
  # Grouping IDs for ACE computation
  if (geo_large == 'state') {
    lgeom <- suppressMessages(suppressWarnings(tigris::states(year = year)))
    out_dat <- out_dat %>%
      dplyr::mutate(
        oid = STATEFP,
        state = stringr::str_trim(state)
      )
  }
  if (geo_large == 'tract') {
    lgeom <- suppressMessages(suppressWarnings(tigris::tracts(
      year = year, state = unique(out_dat$state)
    )))
    out_dat <- out_dat %>%
      dplyr::mutate(
        oid = paste0(STATEFP, COUNTYFP, TRACTCE),
        state = stringr::str_trim(state),
        county = stringr::str_trim(county)
      )
  }
  if (geo_large == 'county') {
    lgeom <- suppressMessages(suppressWarnings(tigris::counties(
      year = year, state = unique(out_dat$state)
    )))
    out_dat <- out_dat %>%
      dplyr::mutate(
        oid = paste0(STATEFP, COUNTYFP),
        state = stringr::str_trim(state),
        county = stringr::str_trim(county)
      )
  }
  if (geo_large == 'place') {
    stopifnot(is.numeric(year), year >= 2011) # Places only available 2011 onward
    lgeom <- suppressMessages(suppressWarnings(tigris::places(
      year = year, state = unique(out_dat$state)
    )))
    wlgeom <- sf::st_within(out_dat, lgeom)
    out_dat <- out_dat %>%
      dplyr::mutate(
        oid = lapply(wlgeom, function(x) { 
          tmp <- lgeom[x, 4] %>% sf::st_drop_geometry()
          lapply(tmp, function(x) { if (length(x) == 0) NA else x })
        }) %>% 
          unlist(),
        place = lapply(wlgeom, function(x) { 
          tmp <- lgeom[x, 5] %>% sf::st_drop_geometry()
          lapply(tmp, function(x) { if (length(x) == 0) NA else x })
        }) %>% 
          unlist()
      )
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
      )
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
      )
  }
  if (geo_large == 'metro') {
    stopifnot(is.numeric(year), year >= 2011) # Metropolitan Divisions only available 2011 onward
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
      )
  }
  
  # Count of racial or ethnic subgroup populations
  ## Count of racial or ethnic comparison subgroup population
  if (length(in_subgroup) == 1) {
    out_dat <- out_dat %>%
      dplyr::mutate(subgroup = as.data.frame(.)[, in_subgroup])
  } else {
    out_dat <- out_dat %>%
      dplyr::mutate(subgroup = rowSums(as.data.frame(.)[, in_subgroup]))
  }
  
  # Compute ACE
  ## From Duncan, Cuzzort, & Duncan (1961; LC:60007089)
  ## ACE = \left ( \sum_{i=2}^{n}X_{i-1}A_{i} \right ) - \left ( \sum_{i=1}^{n}X_{i}A_{i-1} \right )
  ## Where for i smaller geographical units:
  ## X_{i} is the cumulative proportion of the subgroup population through smaller geographical unit i
  ## A_{i} is the cumulative proportion of land area through smaller geographical unit i
  ## when smaller geographical units are ordered by increasing distance 
  ## from the center of a larger geographical unit
  
  ## Compute
  out_tmp <- out_dat %>%
    .[.$oid != 'NANA', ] %>%
    split(., f = list(.$oid)) %>%
    lapply(., FUN = ace_fun, lgeom = lgeom, crs = crs, omit_NAs = omit_NAs) %>%
    utils::stack(.) %>%
    dplyr::mutate(
      ACE = values,
      oid = ind
    ) %>%
    dplyr::select(ACE, oid) %>%
    sf::st_drop_geometry()
  
  # Warning for missingness of census characteristics
  missingYN <- out_dat[, c('TotalPopE', in_subgroup, 'ALAND')] %>% 
    sf::st_drop_geometry()
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
  out <- out_dat %>%
    sf::st_drop_geometry() %>%
    dplyr::left_join(out_tmp, by = dplyr::join_by(oid))
  if (geo_large == 'state') {
    out <- out %>%
      dplyr::select(oid, state, ACE) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, state, ACE)
  }
  if (geo_large == 'county') {
    out <- out %>%
      dplyr::select(oid, state, county, ACE) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, state, county, ACE)
  }
  if (geo_large == 'tract') {
    out <- out %>%
      dplyr::select(oid, state, county, tract, ACE) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, state, county, tract, ACE)
  }
  if (geo_large == 'place') {
    out <- out %>%
      dplyr::select(oid, place, ACE) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, place, ACE)
  }
  if (geo_large == 'cbsa') {
    out <- out %>%
      dplyr::select(oid, cbsa, ACE) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, cbsa, ACE)
  }
  if (geo_large == 'csa') {
    out <- out %>%
      dplyr::select(oid, csa, ACE) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, csa, ACE)
  }
  if (geo_large == 'metro') {
    out <- out %>%
      dplyr::select(oid, metro, ACE) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, metro, ACE)
  }
  
  out <- out %>%
    .[.$GEOID != 'NANA', ] %>%
    dplyr::filter(!is.na(GEOID)) %>%
    dplyr::distinct(GEOID, .keep_all = TRUE) %>%
    dplyr::arrange(GEOID) %>%
    dplyr::as_tibble()
  
  out_dat <- out_dat %>%
    dplyr::arrange(GEOID) %>%
    dplyr::as_tibble()
  
  out <- list(ace = out, ace_data = out_dat, missing = missingYN)
  
  return(out)
}
