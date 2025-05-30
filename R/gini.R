#' Gini Index based on Gini (1921)
#'
#' Compute the aspatial racial or ethnic Gini Index and retrieve the aspatial income Gini Index 
#'
#' @param geo_large Character string specifying the larger geographical unit of the data. The default is counties \code{geo_large = 'county'}.
#' @param geo_small Character string specifying the smaller geographical unit of the data. The default is census tracts \code{geo_small = 'tract'}.
#' @param year Numeric. The year to compute the estimate. The default is 2020, and the years 2009 onward are currently available.
#' @param subgroup Character string specifying the racial or ethnic subgroup(s). See Details for available choices.
#' @param omit_NAs Logical. If FALSE, will compute index for a larger geographical unit only if all of its smaller geographical units have values. The default is TRUE.
#' @param quiet Logical. If TRUE, will display messages about potential missing census information
#' @param ... Arguments passed to \code{\link[tidycensus]{get_acs}} to select state, county, and other arguments for census characteristics
#'
#' @details This function will retrieve the aspatial Gini Index (\emph{G}) of U.S. census tracts or counties for a specified geographical extent (e.g., the entire U.S. or a single state) based on Gini (1921) \doi{10.2307/2223319} for income inequality (at smaller geographical units) and race or ethnicity inequality (at larger geographical units).
#'
#' The function uses the \code{\link[tidycensus]{get_acs}} function to obtain U.S. Census Bureau 5-year American Community Survey estimates of \emph{G} for the geospatial computation. The yearly estimates are available for 2009 onward when ACS-5 data are available (2010 onward for \code{geo_large = 'cbsa'} and 2011 onward for \code{geo_large = 'place'}, \code{geo_large = 'csa'}, or \code{geo_large = 'metro'}) but are available from other U.S. Census Bureau surveys. The function will retrieve the provided income inequality metric (\strong{B19083}) and the twenty racial or ethnic subgroups (U.S. Census Bureau definitions) are:
#' \itemize{
#'  \item \strong{B03002_002}: not Hispanic or Latino \code{'NHoL'}
#'  \item \strong{B03002_003}: not Hispanic or Latino, white alone\code{'NHoLW'}
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
#' According to the U.S. Census Bureau \url{https://www.census.gov/topics/income-poverty/income-inequality/about/metrics/gini-index.html}: 'The Gini Index is a summary measure of income inequality. The Gini coefficient incorporates the detailed shares data into a single statistic, which summarizes the dispersion of income across the entire income distribution. The Gini coefficient ranges from 0, indicating perfect equality (where everyone receives an equal share), to 1, perfect inequality (where only one recipient or group of recipients receives all the income). The Gini Index is based on the difference between the Lorenz curve (the observed cumulative income distribution) and the notion of a perfectly equal income distribution.' For racial or ethnic inequality, *G* is a summary measure of racial or ethnic unevenness or the mean absolute difference between a selected subgroup proportions weighted across all pairs of geographic units, expressed as a proportion of the maximum weighted difference. 
#'
#' Larger geographical units available include states \code{geo_large = 'state'}, counties \code{geo_large = 'county'}, census tracts \code{geo_large = 'tract'}, census-designated places \code{geo_large = 'place'}, core-based statistical areas \code{geo_large = 'cbsa'}, combined statistical areas \code{geo_large = 'csa'}, and metropolitan divisions \code{geo_large = 'metro'}. Smaller geographical units available include, counties \code{geo_small = 'county'}, census tracts \code{geo_small = 'tract'}, and census block groups \code{geo_small = 'cbg'}. If a larger geographical unit is comprised of only one smaller geographical unit (e.g., a U.S county contains only one census tract), then the \emph{V} value returned is NA. If the larger geographical unit is census-designated places \code{geo_large = 'place'}, core-based statistical areas \code{geo_large = 'cbsa'}, combined statistical areas \code{geo_large = 'csa'}, or metropolitan divisions \code{geo_large = 'metro'}, only the smaller geographical units completely within a larger geographical unit are considered in the \emph{V} computation (see internal \code{\link[sf]{st_within}} function for more information) and recommend specifying all states within which the interested larger geographical unit are located using the internal \code{state} argument to ensure all appropriate smaller geographical units are included in the \emph{V} computation.
#'
#' @return An object of class 'list'. This is a named list with the following components:
#'
#' \describe{
#' \item{\code{g}}{An object of class 'tbl' for the GEOID, name, and \emph{G_re} metrics of specified census geographies.}
#' \item{\code{g_data}}{An object of class 'tbl' for the raw census values at specified smaller census geographies including \emph{G_inc}.}
#' \item{\code{missing}}{An object of class 'tbl' of the count and proportion of missingness for \emph{G_inc} and each census variable used to compute \emph{G_re}.}
#' }
#' 
#' @seealso \code{\link[tidycensus]{get_acs}} for additional arguments for geographic extent selection (i.e., \code{state} and \code{county}).
#' @seealso Other one-group evenness indices: \code{\link{atkinson}}, \code{\link{james_taeuber}}, \code{\link{sudano}}, \code{\link{theil}}
#' @seealso Between groups dissimilarity indices: \code{\link{duncan}}
#' 
#' @references Gini, C (1921) Measurement of Inequality of Incomes. \emph{The Economic Journal}, 31(121):124-126. \doi{10.2307/2223319}
#' @references Duncan, OD, & Duncan, B (1955) Residential Distribution and Occupational Stratification. \emph{American Journal of Sociology}, 60(5):493-503. \doi{10.2307/2088328}
#' @references Massey, DS, & Denton, NA (1988) The Dimensions of Residential Segregation. \emph{Social Forces}, 67(1):281-315. \doi{10.1093/sf/67.2.281}
#' 
#' @import dplyr
#' @importFrom sf st_drop_geometry st_within
#' @importFrom stats complete.cases
#' @importFrom stringr str_trim
#' @importFrom tidycensus get_acs
#' @importFrom tidyr pivot_longer separate
#' @importFrom tigris combined_statistical_areas core_based_statistical_areas metro_divisions places
#' @importFrom utils stack
#' @export
#' 
#' @examples
#' \dontrun{
#' # Wrapped in \dontrun{} because these examples require a Census API key.
#'
#'   # Gini Index (a metric of evenness) 
#'   ## of Black populations
#'   ## in census tracts of Georgia, U.S.A. (2020)
#'   gini(
#'     geo_large = 'county',
#'     geo_small = 'tract', 
#'     state = 'GA',
#'     year = 2020, 
#'     subgroup = c('NHoLB', 'HoLB')
#'    )
#'    
#' }
#'
gini <- function(geo_large = 'county',
                 geo_small = 'tract',
                 year = 2020,
                 subgroup,
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
    HoLTReSOR = 'B03002_021',
    G_inc = 'B19083_001'
  )
  
  selected_vars <- vars[c('G_inc', 'TotalPop', subgroup)]
  out_names <- names(selected_vars) # save for output
  in_subgroup <- paste0(subgroup, 'E')
  
  # Acquire Gvariables and sf geometries
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
        tract = gsub('[^0-9\\.]', '', tract), cbg = gsub('[^0-9\\.]', '', cbg)
      )
  } 
  
  # Grouping IDs for R computation
  if (geo_large == 'state') {
    out_dat <- out_dat %>%
      dplyr::mutate(oid = STATEFP, state = stringr::str_trim(state)) %>% 
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
  if (geo_large == 'place') {
    stopifnot(is.numeric(year), year >= 2011) # Places only available 2011 onward
    lgeom <- suppressMessages(suppressWarnings(tigris::places(
      year = year, state = unique(out_dat$state))
    ))
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
  
  # Compute G for race or ethnicity inequality
  ## From Gini (1921) https://doi.org/10.2307/2223319
  ## G = \sum_{n}^{i=1}\sum_{n}^{j=1}\left [ t_{i}t_{j}\left| p_{i}-p_{j}\right| /2T^{2}P(1-P)\right ]
  ## Where:
  ## t_{i} is the total population of area i
  ## t_{j} is the total population of area j
  ## p_{i} is the proportion of the subgroup population of area i
  ## p_{j} is the proportion of the subgroup population of area j
  ## T is the total population of all smaller geographical units
  ## P is the proportion of the subgroup population of all smaller geographical units
  
  ## Compute
  out_tmp <- out_dat %>%
    .[.$oid != 'NANA', ] %>%
    split(., f = list(.$oid)) %>%
    lapply(., FUN = g_fun, omit_NAs = omit_NAs) %>%
    utils::stack(.) %>%
    dplyr::mutate(G_re = values, oid = ind) %>%
    dplyr::select(G_re, oid)
  
  # Warning for missingness of census characteristics
  missingYN <- out_dat[, c('G_incE', 'TotalPopE', in_subgroup)]
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
    dplyr::left_join(out_tmp, by = dplyr::join_by(oid))
  if (geo_large == 'state') {
    out <- out %>%
      dplyr::select(oid, state, G_re) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, state, G_re)
  }
  if (geo_large == 'county') {
    out <- out %>%
      dplyr::select(oid, state, county, G_re) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, state, county, G_re)
  }
  if (geo_large == 'tract') {
    out <- out %>%
      dplyr::select(oid, state, county, tract, G_re) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, state, county, tract,G_re)
  }
  if (geo_large == 'place') {
    out <- out %>%
      dplyr::select(oid, place, G_re) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, place, G_re)
  }
  if (geo_large == 'cbsa') {
    out <- out %>%
      dplyr::select(oid, cbsa, G_re) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, cbsa, G_re)
  }
  if (geo_large == 'csa') {
    out <- out %>%
      dplyr::select(oid, csa, G_re) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, csa, G_re)
  }
  if (geo_large == 'metro') {
    out <- out %>%
      dplyr::select(oid, metro, G_re) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, metro, G_re)
  }
  
  out <- out %>%
    .[.$GEOID != 'NANA', ] %>%
    dplyr::filter(!is.na(GEOID)) %>%
    dplyr::distinct(GEOID, .keep_all = TRUE) %>%
    dplyr::arrange(GEOID) %>%
    dplyr::as_tibble()
  
  out_dat <- out_dat %>%
    dplyr::rename(G_inc = G_incE) %>%
    dplyr::arrange(GEOID) %>%
    dplyr::as_tibble()
  
  out <- list(g = out, g_data = out_dat, missing = missingYN)
  
  return(out)
}
