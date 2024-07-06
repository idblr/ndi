#' Delta based on Hoover (1941) and Duncan et al. (1961) 
#' 
#' Compute the aspatial Delta (Hoover) of a selected racial/ethnic subgroup(s) and U.S. geographies.
#'
#' @param geo_large Character string specifying the larger geographical unit of the data. The default is counties \code{geo_large = 'county'}.
#' @param geo_small Character string specifying the smaller geographical unit of the data. The default is census tracts \code{geo_large = 'tract'}.
#' @param year Numeric. The year to compute the estimate. The default is 2020, and the years 2009 onward are currently available.
#' @param subgroup Character string specifying the racial/ethnic subgroup(s). See Details for available choices.
#' @param omit_NAs Logical. If FALSE, will compute index for a larger geographical unit only if all of its smaller geographical units have values. The default is TRUE.
#' @param quiet Logical. If TRUE, will display messages about potential missing census information. The default is FALSE.
#' @param ... Arguments passed to \code{\link[tidycensus]{get_acs}} to select state, county, and other arguments for census characteristics
#'
#' @details This function will compute the aspatial Delta (DEL) of selected racial/ethnic subgroups and U.S. geographies for a specified geographical extent (e.g., the entire U.S. or a single state) based on Hoover (1941) \doi{10.1017/S0022050700052980} and Duncan, Cuzzort, and Duncan (1961; LC:60007089). This function provides the computation of DEL for any of the U.S. Census Bureau race/ethnicity subgroups (including Hispanic and non-Hispanic individuals).
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
#' DEL is a measure of the proportion of members of one subgroup(s) residing in geographic units with above average density of members of the subgroup(s). The index provides the proportion of a subgroup population that would have to move across geographic units to achieve a uniform density. DEL can range in value from 0 to 1.
#' 
#' Larger geographies available include state \code{geo_large = 'state'}, county \code{geo_large = 'county'}, and census tract \code{geo_large = 'tract'} levels. Smaller geographies available include, county \code{geo_small = 'county'}, census tract \code{geo_small = 'tract'}, and census block group \code{geo_small = 'block group'} levels. If a larger geographical area is comprised of only one smaller geographical area (e.g., a U.S county contains only one census tract), then the DEL value returned is NA.
#' 
#' @return An object of class 'list'. This is a named list with the following components:
#' 
#' \describe{
#' \item{\code{del}}{An object of class 'tbl' for the GEOID, name, and DEL at specified larger census geographies.}
#' \item{\code{del_data}}{An object of class 'tbl' for the raw census values at specified smaller census geographies.}
#' \item{\code{missing}}{An object of class 'tbl' of the count and proportion of missingness for each census variable used to compute DEL.}
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
#'   # Delta (a measure of concentration) of non-Hispanic Black vs. non-Hispanic white populations
#'   ## of census tracts within Georgia, U.S.A., counties (2020)
#'   hoover(
#'     geo_large = 'county',
#'     geo_small = 'tract',
#'     state = 'GA',
#'     year = 2020,
#'     subgroup = 'NHoLB'
#'    )
#'   
#' }
#' 
hoover <- function(geo_large = 'county',
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
  
  selected_vars <- vars[subgroup]
  out_names <- c(names(selected_vars), 'ALAND') # save for output
  in_subgroup <- paste(subgroup, 'E', sep = '')
  
  # Acquire DEL variables and sf geometries
  del_data <- suppressMessages(suppressWarnings(
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
    del_data <- del_data %>%
      sf::st_drop_geometry() %>%
      tidyr::separate(NAME.y, into = c('county', 'state'), sep = ',')
  }
  if (geo_small == 'tract') {
    del_data <- del_data %>%
      sf::st_drop_geometry() %>%
      tidyr::separate(NAME.y, into = c('tract', 'county', 'state'), sep = ',') %>%
      dplyr::mutate(tract = gsub('[^0-9\\.]', '', tract))
  } 
  if (geo_small == 'block group') {
    del_data <- del_data %>%
      sf::st_drop_geometry() %>%
      tidyr::separate(NAME.y, into = c('block.group', 'tract', 'county', 'state'), sep = ',') %>%
      dplyr::mutate(
        tract = gsub('[^0-9\\.]', '', tract), block.group = gsub('[^0-9\\.]', '', block.group)
      )
  } 
  
  # Grouping IDs for DEL computation
  if (geo_large == 'tract') {
    del_data <- del_data %>%
      dplyr::mutate(
        oid = paste(.$STATEFP, .$COUNTYFP, .$TRACTCE, sep = ''),
        state = stringr::str_trim(state),
        county = stringr::str_trim(county)
      )
  }
  if (geo_large == 'county') {
    del_data <- del_data %>%
      dplyr::mutate(
        oid = paste(.$STATEFP, .$COUNTYFP, sep = ''),
        state = stringr::str_trim(state),
        county = stringr::str_trim(county)
      )
  }
  if (geo_large == 'state') {
    del_data <- del_data %>%
      dplyr::mutate(
        oid = .$STATEFP,
        state = stringr::str_trim(state)
      )
  }
  
  # Count of racial/ethnic subgroup populations
  ## Count of racial/ethnic comparison subgroup population
  if (length(in_subgroup) == 1) {
    del_data <- del_data %>%
      dplyr::mutate(subgroup = .[ , in_subgroup])
  } else {
    del_data <- del_data %>%
      dplyr::mutate(subgroup = rowSums(.[ , in_subgroup]))
  }
  
  # Compute DEL
  ## From Hoover (1961) https://10.1017/S0022050700052980
  ## 0.5\sum_{i=1}^{n}\left|\frac{x_{i}}{X}-\frac{a_{i}}{A}\right|
  ## Where for k geographical units i:
  ## X denotes the total number of subgroup population in study (reference) area
  ## x_{i} denotes the number of subgroup population X in geographical unit i
  ## A denotes the total land area in study (reference) area (sum of all a_{i}
  ## a_{i} denotes the land area of geographical unit i
  
  ## Compute
  DELtmp <- del_data %>%
    split(., f = list(del_data$oid)) %>%
    lapply(., FUN = del_fun, omit_NAs = omit_NAs) %>%
    utils::stack(.) %>%
    dplyr::mutate(DEL = values, oid = ind) %>%
    dplyr::select(DEL, oid)
  
  # Warning for missingness of census characteristics
  missingYN <- del_data[ , c(in_subgroup, 'ALAND')]
  names(missingYN) <- out_names
  missingYN <- missingYN %>%
    tidyr::pivot_longer(cols = dplyr::everything(), names_to = 'variable', values_to = 'val') %>%
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
    del <- del_data %>%
      dplyr::left_join(DELtmp, by = dplyr::join_by(oid)) %>%
      dplyr::select(oid, state, DEL) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, state, DEL) %>%
      .[.$GEOID != 'NANA', ]
  }
  if (geo_large == 'county') {
    del <- del_data %>%
      dplyr::left_join(DELtmp, by = dplyr::join_by(oid)) %>%
      dplyr::select(oid, state, county, DEL) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, state, county, DEL) %>%
      .[.$GEOID != 'NANA', ]
  }
  if (geo_large == 'tract') {
    del <- del_data %>%
      dplyr::left_join(DELtmp, by = dplyr::join_by(oid)) %>%
      dplyr::select(oid, state, county, tract, DEL) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, state, county, tract, DEL) %>%
      .[.$GEOID != 'NANA', ]
  }
  
  del <- del %>%
    dplyr::arrange(GEOID) %>%
    dplyr::as_tibble()
  
  del_data <- del_data %>%
    dplyr::arrange(GEOID) %>%
    dplyr::as_tibble() 
  
  out <- list(del = del, del_data = del_data, missing = missingYN)
  
  return(out)
}
