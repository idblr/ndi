#' Local Exposure and Isolation metric based on Bemanian & Beyer (2017)
#'
#' Compute the aspatial Local Exposure and Isolation (Bemanian & Beyer) metric of a selected racial/ethnic subgroup(s) and U.S. geographies.
#'
#' @param geo_large Character string specifying the larger geographical unit of the data. The default is counties \code{geo_large = 'county'}.
#' @param geo_small Character string specifying the smaller geographical unit of the data. The default is census tracts \code{geo_large = 'tract'}.
#' @param year Numeric. The year to compute the estimate. The default is 2020, and the years 2009 onward are currently available.
#' @param subgroup Character string specifying the racial/ethnic subgroup(s). See Details for available choices.
#' @param subgroup_ixn Character string specifying the racial/ethnic subgroup(s) as the interaction population. If the same as \code{subgroup}, will compute the simple isolation of the group. See Details for available choices.
#' @param omit_NAs Logical. If FALSE, will compute index for a larger geographical unit only if all of its smaller geographical units have values. The default is TRUE.
#' @param quiet Logical. If TRUE, will display messages about potential missing census information. The default is FALSE.
#' @param ... Arguments passed to \code{\link[tidycensus]{get_acs}} to select state, county, and other arguments for census characteristics
#'
#' @details This function will compute the aspatial Local Exposure and Isolation (\emph{LEx/Is}) metric of selected racial/ethnic subgroups and U.S. geographies for a specified geographical extent (e.g., the entire U.S. or a single state) based on Bemanian & Beyer (2017) \doi{10.1158/1055-9965.EPI-16-0926}. This function provides the computation of \emph{LEx/Is} for any of the U.S. Census Bureau race/ethnicity subgroups (including Hispanic and non-Hispanic individuals).
#'
#' The function uses the \code{\link[tidycensus]{get_acs}} function to obtain U.S. Census Bureau 5-year American Community Survey characteristics used for the aspatial computation. The yearly estimates are available for 2009 onward when ACS-5 data are available (2010 onward for \code{geo_large = 'cbsa'} and 2011 onward for \code{geo_large = 'csa'} or \code{geo_large = 'metro'}) but may be available from other U.S. Census Bureau surveys. The twenty racial/ethnic subgroups (U.S. Census Bureau definitions) are:
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
#' \emph{LEx/Is} is a measure of the probability that two individuals living within a specific smaller geography (e.g., census tract) of either different (i.e., exposure) or the same (i.e., isolation) racial/ethnic subgroup(s) will interact, assuming that individuals within a smaller geography are randomly mixed. \emph{LEx/Is} is standardized with a logit transformation and centered against an expected case that all races/ethnicities are evenly distributed across a larger geography. (Note: will adjust data by 0.025 if probabilities are zero, one, or undefined. The output will include a warning if adjusted. See \code{\link[car]{logit}} for additional details.)
#'
#' \emph{LEx/Is} can range from negative infinity to infinity. If \emph{LEx/Is} is zero then the estimated probability of the interaction between two people of the given subgroup(s) within a smaller geography is equal to the expected probability if the subgroup(s) were perfectly mixed in the larger geography. If \emph{LEx/Is} is greater than zero then the interaction is more likely to occur within the smaller geography than in the larger geography, and if \emph{LEx/Is} is less than zero then the interaction is less likely to occur within the smaller geography than in the larger geography. Note: the exponentiation of each \emph{LEx/Is} metric results in the odds ratio of the specific exposure or isolation of interest in a smaller geography relative to the larger geography.
#'
#' Larger geographies available include state \code{geo_large = 'state'}, county \code{geo_large = 'county'}, census tract \code{geo_large = 'tract'}, Core Based Statistical Area \code{geo_large = 'cbsa'}, Combined Statistical Area \code{geo_large = 'csa'}, and Metropolitan Division \code{geo_large = 'metro'} levels. Smaller geographies available include, county \code{geo_small = 'county'}, census tract \code{geo_small = 'tract'}, and census block group \code{geo_small = 'block group'} levels. If a larger geographical area is comprised of only one smaller geographical area (e.g., a U.S county contains only one census tract), then the \emph{LEx/Is} value returned is NA. If the larger geographical unit is Combined Based Statistical Areas \code{geo_large = 'csa'} or Core Based Statistical Areas \code{geo_large = 'cbsa'}, only the smaller geographical units completely within a larger geographical unit are considered in the \emph{LEx/Is} computation (see internal \code{\link[sf]{st_within}} function for more information) and recommend specifying all states within which the interested larger geographical unit are located using the internal \code{state} argument to ensure all appropriate smaller geographical units are included in the \emph{LEx/Is} computation.
#'
#' @return An object of class 'list'. This is a named list with the following components:
#'
#' \describe{
#' \item{\code{lexis}}{An object of class 'tbl' for the GEOID, name, and \emph{LEx/Is} at specified smaller census geographies.}
#' \item{\code{lexis_data}}{An object of class 'tbl' for the raw census values at specified smaller census geographies.}
#' \item{\code{missing}}{An object of class 'tbl' of the count and proportion of missingness for each census variable used to compute \emph{LEx/Is}.}
#' }
#'
#' @import dplyr
#' @importFrom car logit
#' @importFrom sf st_drop_geometry st_within
#' @importFrom stats complete.cases
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
#'   # Isolation of non-Hispanic Black vs. non-Hispanic white populations
#'   ## of census tracts within Georgia, U.S.A., counties (2020)
#'   bemanian_beyer(
#'     geo_large = 'county',
#'     geo_small = 'tract',
#'     state = 'GA',
#'     year = 2020,
#'     subgroup = 'NHoLB',
#'     subgroup_ixn = 'NHoLW'
#'    )
#'
#' }
#'
bemanian_beyer <- function(geo_large = 'county',
                           geo_small = 'tract',
                           year = 2020,
                           subgroup,
                           subgroup_ixn,
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
    match.arg(
      subgroup_ixn,
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
    
    selected_vars <- vars[c('TotalPop', subgroup, subgroup_ixn)]
    out_names <- names(selected_vars) # save for output
    in_subgroup <- paste0(subgroup, 'E')
    in_subgroup_ixn <- paste0(subgroup_ixn, 'E')
    
    # Acquire LEx/Is variables and sf geometries
    lexis_data <- suppressMessages(suppressWarnings(
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
      lexis_data <- lexis_data %>% 
        tidyr::separate(NAME.y, into = c('county', 'state'), sep = ',')
    }
    if (geo_small == 'tract') {
      lexis_data <- lexis_data %>% 
        tidyr::separate(NAME.y, into = c('tract', 'county', 'state'), sep = ',') %>%
        dplyr::mutate(tract = gsub('[^0-9\\.]', '', tract))
    }
    if (geo_small == 'block group') {
      lexis_data <- lexis_data %>% 
        tidyr::separate(NAME.y, into = c('block.group', 'tract', 'county', 'state'), sep = ',') %>%
        dplyr::mutate(
          tract = gsub('[^0-9\\.]', '', tract),
          block.group = gsub('[^0-9\\.]', '', block.group)
        )
    }
    
    # Grouping IDs for LEx/Is computation
    if (geo_large == 'state') {
      lexis_data <- lexis_data %>%
        dplyr::mutate(
          oid = STATEFP,
          state = stringr::str_trim(state)
        ) %>% 
        sf::st_drop_geometry()
    }
    if (geo_large == 'county') {
      lexis_data <- lexis_data %>%
        dplyr::mutate(
          oid = paste(STATEFP, COUNTYFP, sep = ''),
          state = stringr::str_trim(state),
          county = stringr::str_trim(county)
        ) %>% 
        sf::st_drop_geometry()
    }
    if (geo_large == 'tract') {
      lexis_data <- lexis_data %>%
        dplyr::mutate(
          oid = paste(STATEFP, COUNTYFP, TRACTCE, sep = ''),
          state = stringr::str_trim(state),
          county = stringr::str_trim(county)
        ) %>% 
        sf::st_drop_geometry()
    }
    if (geo_large == 'cbsa') {
      stopifnot(is.numeric(year), year >= 2010) # CBSAs only available 2010 onward
      dat_cbsa <- suppressMessages(suppressWarnings(tigris::core_based_statistical_areas(year = year)))
      win_cbsa <- sf::st_within(lexis_data, dat_cbsa)
      lexis_data <- lexis_data %>%
        dplyr::mutate(
          oid = lapply(win_cbsa, function(x) { 
            tmp <- dat_cbsa[x, 3] %>% sf::st_drop_geometry()
            lapply(tmp, function(x) { if (length(x) == 0) NA else x })
          }) %>% 
            unlist(),
          cbsa = lapply(win_cbsa, function(x) { 
            tmp <- dat_cbsa[x, 4] %>% sf::st_drop_geometry()
            lapply(tmp, function(x) { if (length(x) == 0) NA else x })
          }) %>% 
            unlist()
        ) %>% 
        sf::st_drop_geometry()
    }
    if (geo_large == 'csa') {
      stopifnot(is.numeric(year), year >= 2011) # CSAs only available 2011 onward
      dat_csa <- suppressMessages(suppressWarnings(tigris::combined_statistical_areas(year = year)))
      win_csa <- sf::st_within(lexis_data, dat_csa)
      lexis_data <- lexis_data %>%
        dplyr::mutate(
          oid = lapply(win_csa, function(x) { 
            tmp <- dat_csa[x, 2] %>% sf::st_drop_geometry()
            lapply(tmp, function(x) { if (length(x) == 0) NA else x })
          }) %>% 
            unlist(),
          csa = lapply(win_csa, function(x) { 
            tmp <- dat_csa[x, 3] %>% sf::st_drop_geometry()
            lapply(tmp, function(x) { if (length(x) == 0) NA else x })
          }) %>% 
            unlist()
        ) %>% 
        sf::st_drop_geometry()
    }
    if (geo_large == 'metro') {
      stopifnot(is.numeric(year), year >= 2011) # Metro Divisions only available 2011 onward
      dat_metro <- suppressMessages(suppressWarnings(tigris::metro_divisions(year = year)))
      win_metro <- sf::st_within(lexis_data, dat_metro)
      lexis_data <- lexis_data %>%
        dplyr::mutate(
          oid = lapply(win_metro, function(x) { 
            tmp <- dat_metro[x, 4] %>% sf::st_drop_geometry()
            lapply(tmp, function(x) { if (length(x) == 0) NA else x })
          }) %>% 
            unlist(),
          metro = lapply(win_metro, function(x) { 
            tmp <- dat_metro[x, 5] %>% sf::st_drop_geometry()
            lapply(tmp, function(x) { if (length(x) == 0) NA else x })
          }) %>% 
            unlist()
        ) %>% 
        sf::st_drop_geometry()
    }
    
    # Count of racial/ethnic subgroup populations
    ## Count of racial/ethnic comparison subgroup population
    if (length(in_subgroup) == 1) {
      lexis_data <- lexis_data %>%
        dplyr::mutate(subgroup = .[, in_subgroup])
    } else {
      lexis_data <- lexis_data %>%
        dplyr::mutate(subgroup = rowSums(.[, in_subgroup]))
    }
    ## Count of racial/ethnic interaction subgroup population
    if (length(in_subgroup_ixn) == 1) {
      lexis_data <- lexis_data %>%
        dplyr::mutate(subgroup_ixn = .[, in_subgroup_ixn])
    } else {
      lexis_data <- lexis_data %>%
        dplyr::mutate(subgroup_ixn = rowSums(.[, in_subgroup_ixn]))
    }
    
    # Compute LEx/Is
    ## From Bemanian & Beyer (2017) https://doi.org/10.1158/1055-9965.EPI-16-0926
    ## E^*_{m,n}(i) = log\left(\frac{p_{im} \times p_{in}}{1 - p_{im} \times p_{in}}\right) - log\left(\frac{P_{m} \times P_{n}}{1 - P_{m} \times P_{n}}\right)
    ## Where for smaller geographical unit i:
    ## p_{im} denotes the number of subgroup population m in smaller geographical unit i
    ## p_{in} denotes the number of subgroup population n in smaller geographical unit i
    ## P_{m} denotes the number of subgroup population m in larger geographical unit within which the smaller geographic unit i is located
    ## P_{n} denotes the number of subgroup population n in larger geographical unit within which the smaller geographic unit i is located
    ## If m \ne n, then computes the exposure of members of subgroup populations m and n
    ## If m = n, then computes the simple isolation experienced by members of subgroup population m
    
    ## Compute
    LExIstmp <- lexis_data %>%
      split(., f = list(lexis_data$oid)) %>%
      lapply(., FUN = lexis_fun, omit_NAs = omit_NAs) %>%
      do.call('rbind', .)
    
    # Warning for missingness of census characteristics
    missingYN <- lexis_data[, c('TotalPopE', in_subgroup, in_subgroup_ixn)]
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
    lexis <- lexis_data %>%
      dplyr::left_join(LExIstmp, by = dplyr::join_by(GEOID))
    
    if (geo_small == 'state') {
      lexis <- lexis %>%
        dplyr::select(GEOID, state, LExIs)
    }
    if (geo_small == 'county') {
      lexis <- lexis %>%
        dplyr::select(GEOID, state, county, LExIs)
    }
    if (geo_small == 'tract') {
      lexis <- lexis %>%
        dplyr::select(GEOID, state, county, tract, LExIs)
    }
    if (geo_small == 'block group') {
      lexis <- lexis %>%
        dplyr::select(GEOID, state, county, tract, block.group, LExIs)
    }
    if (geo_large == 'cbsa') {
      lexis <- lexis_data %>%
        dplyr::select(GEOID, cbsa) %>%
        dplyr::left_join(lexis, ., by = dplyr::join_by(GEOID)) %>%
        dplyr::relocate(cbsa, .after = county)
    }
    if (geo_large == 'csa') {
      lexis <- lexis_data %>%
        dplyr::select(GEOID, csa) %>%
        dplyr::left_join(lexis, ., by = dplyr::join_by(GEOID)) %>%
        dplyr::relocate(csa, .after = county)
    }
    if (geo_large == 'metro') {
      lexis <- lexis_data %>%
        dplyr::select(GEOID, metro) %>%
        dplyr::left_join(lexis, ., by = dplyr::join_by(GEOID)) %>%
        dplyr::relocate(metro, .after = county)
    }
    
    lexis <- lexis %>%
      unique(.) %>%
      .[.$GEOID != 'NANA',] %>%
      dplyr::arrange(GEOID) %>%
      dplyr::as_tibble()
    
    lexis_data <- lexis_data %>%
      dplyr::arrange(GEOID) %>%
      dplyr::as_tibble()
    
    out <- list(lexis = lexis, lexis_data = lexis_data, missing = missingYN)
    
    return(out)
  }
