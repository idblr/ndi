#' Racial Isolation Index based on Anthopolos et al. (2011)
#'
#' Compute the spatial Racial Isolation Index (Anthopolos) of selected subgroup(s).
#'
#' @param geo Character string specifying the geography of the data either counties \code{geo = 'county'}, census tracts \code{geo = 'tract'} (the default), or census block groups \code{geo = 'cbg'}.
#' @param year Numeric. The year to compute the estimate. The default is 2020, and the years 2009 onward are currently available.
#' @param subgroup Character string specifying the racial or ethnic subgroup(s). See Details for available choices.
#' @param quiet Logical. If TRUE, will display messages about potential missing census information. The default is FALSE.
#' @param ... Arguments passed to \code{\link[tidycensus]{get_acs}} to select state, county, and other arguments for census characteristics
#'
#' @details This function will compute the spatial Racial Isolation Index (\emph{RI}) of U.S. census tracts or counties for a specified geographical extent (e.g., the entire U.S. or a single state) based on Anthopolos et al. (2011) \doi{10.1016/j.sste.2011.06.002} who originally designed the metric for the racial isolation of non-Hispanic Black individuals. This function provides the computation of \emph{RI} for any of the U.S. Census Bureau race or ethnicity subgroups (including Hispanic and non-Hispanic individuals).
#'
#' The function uses the \code{\link[tidycensus]{get_acs}} function to obtain U.S. Census Bureau 5-year American Community Survey characteristics used for the geospatial computation. The yearly estimates are available for 2009 onward when ACS-5 data are available but are available from other U.S. Census Bureau surveys. The twenty racial or ethnic subgroups (U.S. Census Bureau definitions) are:
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
#' Use the internal \code{state} and \code{county} arguments within the \code{\link[tidycensus]{get_acs}} function to specify geographic extent of the data output. NOTE: Current version does not correct for edge effects (e.g., census geographies along the specified spatial extent border, coastline, or U.S.-Mexico / U.S.-Canada border) may have few neighboring census geographies, and \emph{RI} values in these census geographies may be unstable. A stop-gap solution for the former source of edge effect is to compute the \emph{RI} for neighboring census geographies (i.e., the states bordering a study area of interest) and then use the estimates of the study area of interest.
#'
#' A census geography (and its neighbors) that has nearly all of its population who identify with the specified race or ethnicity subgroup(s) (e.g., non-Hispanic or Latino, Black or African American alone) will have an \emph{RI} value close to 1. In contrast, a census geography (and its neighbors) that has nearly none of its population who identify with the specified race or ethnicity subgroup(s) (e.g., not non-Hispanic or Latino, Black or African American alone) will have an \emph{RI} value close to 0.
#'
#' @return An object of class 'list'. This is a named list with the following components:
#'
#' \describe{
#' \item{\code{ri}}{An object of class 'tbl' for the GEOID, name, \emph{RI}, and raw census values of specified census geographies.}
#' \item{\code{missing}}{An object of class 'tbl' of the count and proportion of missingness for each census variable used to compute \emph{RI}.}
#' }
#'
#' @import dplyr
#' @importFrom Matrix sparseMatrix
#' @importFrom sf st_drop_geometry st_geometry st_intersects
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
#'   # Racial Isolation Index of Black populations
#'   ## of census tracts within Georgia, U.S.A., counties (2020)
#'   anthopolos(
#'     geo = 'tract',
#'     state = 'GA',
#'     year = 2020,
#'     subgroup = c('NHoLB', 'HoLB')
#'    )
#'
#' }
#'
anthopolos <- function(geo = 'tract',
                       year = 2020,
                       subgroup,
                       quiet = FALSE,
                       ...) {
  
  # Check arguments
  match.arg(geo, choices = c('county', 'tract', 'cbg', 'block group'))
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
  prefix <- 'subgroup'
  suffix <- seq(1:length(subgroup))
  names(selected_vars) <- c('TotalPop', paste0(prefix, suffix))
  in_names <- paste0(names(selected_vars), 'E')
  
  # Acquire RI variables and sf geometries
  out_dat <- suppressMessages(suppressWarnings(
    tidycensus::get_acs(
      geography = geo,
      year = year,
      output = 'wide',
      variables = selected_vars,
      geometry = TRUE,
      ...
    )
  ))
  
  # Format output
  if (geo == 'county') {
    out_dat <- out_dat %>%
      tidyr::separate(NAME, into = c('county', 'state'), sep = ',')
  }
  if (geo == 'tract') {
    out_dat <- out_dat %>%
      tidyr::separate(NAME, into = c('tract', 'county', 'state'), sep = ',') %>%
      dplyr::mutate(tract = gsub('[^0-9\\.]', '', tract))
  } 
  if (geo == 'cbg' | geo == 'block group') {
    out_dat <- out_dat %>%
      tidyr::separate(NAME, into = c('cbg', 'tract', 'county', 'state'), sep = ',') %>%
      dplyr::mutate(
        tract = gsub('[^0-9\\.]', '', tract), cbg = gsub('[^0-9\\.]', '', cbg)
      )
  } 
  
  out_dat <- out_dat %>%
    dplyr::mutate(subgroup = rowSums(sf::st_drop_geometry(out_dat[, in_names[-1]])))
  
  # Compute RI
  ## From Anthopolos et al. (2011) https://doi.org/10.1016/j.sste.2011.06.002
  ## RI_{im} = (Sigma_{j∈∂_{i}} w_{ij} * T_{jm}) / (Sigma_{j∈∂_{i}} w_{ij} * T_{j})
  ## Where:
  ## ∂_{i} denotes the set of index units i and its neighbors
  ## Given M mutually exclusive racial or ethnic subgroups, m indexes the subgroups of M
  ## T_{i} denotes the total population in region i (TotalPop)
  ## T_{im} denotes the population of the selected subgroup(s) (subgroup1, ...)
  ## w_{ij} denotes a nXn first-order adjacency matrix, where n is the number of census geometries in the study area
  ### and the entries of w_{ij} are set to 1 if a boundary is shared by region i and region j and zero otherwise
  ### Entries of the main diagonal (since i∈∂_{i}, w_{ij} = w_{ii} when j = i) of w_{ij} are set to 1.5
  ### such that the weight of the index unit, i, is larger than the weights assigned to adjacent tracts
  
  ## Geospatial adjacency matrix (w_ij)
  tmp <- out_dat %>%
    sf::st_geometry() %>%
    sf::st_intersects(sparse = TRUE)
  names(tmp) <- as.character(seq_len(nrow(out_dat)))
  tmp_L <- length(tmp)
  tmp_counts <- unlist(Map(length, tmp))
  tmp_i <- rep(1:tmp_L, tmp_counts)
  tmp_j <- unlist(tmp)
  w_ij <- Matrix::sparseMatrix(
    i = tmp_i,
    j = tmp_j,
    x = 1,
    dims = c(tmp_L, tmp_L)
  )
  diag(w_ij) <- 1.5
  
  ## Compute
  out_dat <- out_dat %>%
    sf::st_drop_geometry() # drop geometries (can join back later)
  out_tmp <- list()
  for (i in 1:dim(w_ij)[1]) {
    out_tmp[[i]] <- sum(as.matrix(w_ij[i,]) * out_dat[, 'subgroup']) /
      sum(as.matrix(w_ij[i,]) * out_dat[, 'TotalPopE'])
  }
  out_dat$RI <- unlist(out_tmp)
  
  # Warning for missingness of census characteristics
  missingYN <- out_dat[, in_names]
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
  if (geo == 'county') {
    out <- out_dat %>%
      dplyr::select(c(GEOID, state, county, RI, dplyr::all_of(in_names)))
    names(out) <- c('GEOID', 'state', 'county', 'RI', out_names)
  }
  if (geo == 'tract') {
    out <- out_dat %>%
      dplyr::select(c(GEOID, state, county, tract, RI, dplyr::all_of(in_names)))
    names(out) <- c('GEOID', 'state', 'county', 'tract', 'RI', out_names)
  }
  if (geo == 'cbg' | geo == 'block group') {
    out <- out_dat %>%
      dplyr::select(c(GEOID, state, county, tract, cbg, RI, dplyr::all_of(in_names)))
    names(out) <- c('GEOID', 'state', 'county', 'tract', 'cbg', 'RI', out_names)
  }

  out <- out %>%
    dplyr::mutate(
      state = stringr::str_trim(state),
      county = stringr::str_trim(county)
    ) %>%
    dplyr::arrange(GEOID) %>%
    dplyr::as_tibble()
  
  out <- list(ri = out, missing = missingYN)
  
  return(out)
}
