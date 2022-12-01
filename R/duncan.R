#' Dissimilarity Index based on Duncan & Duncan (1955) 
#' 
#' Compute the Dissimilarity Index (Duncan) values for a selected subgroup(s) and U.S. geographies
#'
#' @param geo_large Character string specifying the larger geographical unit of the data. The default is counties \code{geo_large = "county"}.
#' @param geo_small Character string specifying the smaller geographical unit of the data. The default is census tracts \code{geo_large = "tract"}.
#' @param year Numeric. The year to compute the estimate. The default is 2020, and the years 2009 onward are currently available.
#' @param subgroup Character string specifying the racial/ethnic subgroup(s) as the comparison population. See Details for available choices.
#' @param subgroup_ref Character string specifying the racial/ethnic subgroup(s) as the reference population. See Details for available choices.
#' @param quiet Logical. If TRUE, will display messages about potential missing census information. The default is FALSE.
#' @param ... Arguments passed to \code{\link[tidycensus]{get_acs}} to select state, county, and other arguments for census characteristics
#'
#' @details This function will compute the Dissimilarity Index (DI) of selected racial/ethnic subgroups and U.S. geographies for a specified geographical extent (e.g., the entire U.S. or a single state) based on Duncan & Duncan (1955) \doi{10.2307/2088328}. This function provides the computation of DI for any of the U.S. Census Bureau race/ethnicity subgroups (including Hispanic and non-Hispanic individuals).
#' 
#' The function uses the \code{\link[tidycensus]{get_acs}} function to obtain U.S. Census Bureau 5-year American Community Survey characteristics used for the geospatial computation. The yearly estimates are available for 2009 onward when ACS-5 data are available but are available from other U.S. Census Bureau surveys. The twenty racial/ethnic subgroups (U.S. Census Bureau definitions) are:
#' \itemize{
#'  \item{B03002_002: }{not Hispanic or Latino "NHoL"}
#'  \item{B03002_003: }{not Hispanic or Latino, white alone "NHoLW"}
#'  \item{B03002_004: }{not Hispanic or Latino, Black or African American alone "NHoLB"}
#'  \item{B03002_005: }{not Hispanic or Latino, American Indian and Alaska Native alone "NHoLAIAN"}
#'  \item{B03002_006: }{not Hispanic or Latino, Asian alone "NHoLA"}
#'  \item{B03002_007: }{not Hispanic or Latino, Native Hawaiian and Other Pacific Islander alone "NHoLNHOPI"}
#'  \item{B03002_008: }{not Hispanic or Latino, Some other race alone "NHoLSOR"}
#'  \item{B03002_009: }{not Hispanic or Latino, Two or more races "NHoLTOMR"}
#'  \item{B03002_010: }{not Hispanic or Latino, Two races including Some other race "NHoLTRiSOR"}
#'  \item{B03002_011: }{not Hispanic or Latino, Two races excluding Some other race, and three or more races "NHoLTReSOR"}
#'  \item{B03002_012: }{Hispanic or Latino "HoL"}
#'  \item{B03002_013: }{Hispanic or Latino, white alone "HoLW"}
#'  \item{B03002_014: }{Hispanic or Latino, Black or African American alone "HoLB"}
#'  \item{B03002_015: }{Hispanic or Latino, American Indian and Alaska Native alone "HoLAIAN"}
#'  \item{B03002_016: }{Hispanic or Latino, Asian alone "HoLA"}
#'  \item{B03002_017: }{Hispanic or Latino, Native Hawaiian and Other Pacific Islander alone "HoLNHOPI"}
#'  \item{B03002_018: }{Hispanic or Latino, Some other race alone "HoLSOR"}
#'  \item{B03002_019: }{Hispanic or Latino, Two or more races "HoLTOMR"}
#'  \item{B03002_020: }{Hispanic or Latino, Two races including Some other race "HoLTRiSOR"}
#'  \item{B03002_021: }{Hispanic or Latino, Two races excluding Some other race, and three or more races "HoLTReSOR"}
#' }
#' 
#' Use the internal \code{state} and \code{county} arguments within the \code{\link[tidycensus]{get_acs}} function to specify geographic extent of the data output.
#' 
#' DI is a measure of the evenness of racial/ethnic residential segregation when comparing smaller geographical areas to larger ones within which the smaller geographical areas are located. The DI metric can range in value from 0 to 1 and represents the proportion of racial/ethnic subgroup members that would have to change their area of residence to achieve an even distribution within the larger geographical area under conditions of maximum segregation.
#' 
#' Larger geographies available include state \code{geo_large = "state"}, county \code{geo_large = "county"}, and census tract \code{geo_large = "tract"} levels. Smaller geographies available include, county \code{geo_small = "county"}, census tract \code{geo_small = "tract"}, and census block group \code{geo_small = "block group"} levels. If a larger geographical area is comprised of only one smaller geographical area (e.g., a U.S county contains only one census tract), then the DI value returned is NA.
#' 
#' @return An object of class 'list'. This is a named list with the following components:
#' 
#' \describe{
#' \item{\code{di}}{An object of class 'tbl' for the GEOID, name, and DI metric at specified larger census geographies.}
#' \item{\code{di_data}}{An object of class 'tbl' for the raw census values at specified smaller census geographies.}
#' \item{\code{missing}}{An object of class 'tbl' of the count and proportion of missingness for each census variable used to compute the DI}
#' }
#' 
#' @import dplyr
#' @importFrom sf st_drop_geometry
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
#'   # Dissimilarity Index of non-Hispanic Black vs. non-Hispanic white populations
#'   ## of census tracts within counties (2020)
#'   duncan(geo_large = "county", geo_small = "tract", state = "GA", year = 2020,
#'          subgroup = "NHoLB", subgroup_ref = "NHoLW")
#'   
#' }
#' 
duncan <- function(geo_large = "county", geo_small = "tract", year = 2020, subgroup, subgroup_ref, quiet = FALSE, ...) {
  
 # Check arguments
  match.arg(geo_large, choices = c("state", "county", "tract"))
  match.arg(geo_small, choices = c("county", "tract", "block group"))
  stopifnot(is.numeric(year), year >= 2009) # all variables available 2009 onward
  match.arg(subgroup, several.ok = TRUE,
            choices = c("NHoL", "NHoLW", "NHoLB", "NHoLAIAN", "NHoLA", "NHoLNHOPI",
                        "NHoLSOR", "NHoLTOMR", "NHoLTRiSOR", "NHoLTReSOR",
                        "HoL", "HoLW", "HoLB", "HoLAIAN", "HoLA", "HoLNHOPI",
                        "HoLSOR", "HoLTOMR", "HoLTRiSOR", "HoLTReSOR"))
  match.arg(subgroup_ref, several.ok = TRUE,
            choices = c("NHoL", "NHoLW", "NHoLB", "NHoLAIAN", "NHoLA", "NHoLNHOPI",
                        "NHoLSOR", "NHoLTOMR", "NHoLTRiSOR", "NHoLTReSOR",
                        "HoL", "HoLW", "HoLB", "HoLAIAN", "HoLA", "HoLNHOPI",
                        "HoLSOR", "HoLTOMR", "HoLTRiSOR", "HoLTReSOR"))
  
  # Select census variables
  vars <- c(NHoL = "B03002_002",
            NHoLW = "B03002_003",
            NHoLB = "B03002_004",
            NHoLAIAN = "B03002_005",
            NHoLA = "B03002_006",
            NHoLNHOPI = "B03002_007",
            NHoLSOR = "B03002_008",
            NHoLTOMR = "B03002_009",
            NHoLTRiSOR = "B03002_010",
            NHoLTReSOR = "B03002_011",
            HoL = "B03002_012",
            HoLW = "B03002_013",
            HoLB = "B03002_014",
            HoLAIAN = "B03002_015",
            HoLA = "B03002_016",
            HoLNHOPI = "B03002_017",
            HoLSOR = "B03002_018",
            HoLTOMR = "B03002_019",
            HoLTRiSOR = "B03002_020",
            HoLTReSOR = "B03002_021")
  
  selected_vars <- vars[c(subgroup, subgroup_ref)]
  out_names <- names(selected_vars) # save for output
  in_subgroup <- paste(subgroup, "E", sep = "")
  in_subgroup_ref <- paste(subgroup_ref, "E", sep = "")
  
  # Acquire DI variables and sf geometries
  di_data <- suppressMessages(suppressWarnings(tidycensus::get_acs(geography = geo_small,
                                                                   year = year,
                                                                   output = "wide",
                                                                   variables = selected_vars,
                                                                   geometry = TRUE,
                                                                   keep_geo_vars = TRUE, ...)))
  
  # Format output
  if (geo_small == "county") {
    di_data <- sf::st_drop_geometry(di_data) %>%
      tidyr::separate(NAME.y, into = c("county", "state"), sep = ",")
  }
  if (geo_small == "tract") {
    di_data <- sf::st_drop_geometry(di_data) %>%
      tidyr::separate(NAME.y, into = c("tract", "county", "state"), sep = ",") %>%
      dplyr::mutate(tract = gsub("[^0-9\\.]", "", tract))
  } 
  if (geo_small == "block group") {
    di_data <- sf::st_drop_geometry(di_data) %>%
      tidyr::separate(NAME.y, into = c("block.group", "tract", "county", "state"), sep = ",") %>%
      dplyr::mutate(tract = gsub("[^0-9\\.]", "", tract),
                    block.group = gsub("[^0-9\\.]", "", block.group))
  } 
  
  # Grouping IDs for DI computation
  if (geo_large == "tract") {
    di_vars <- di_data %>%
      dplyr::mutate(oid = paste(.$STATEFP, .$COUNTYFP, .$TRACTCE, sep = ""),
                    state = stringr::str_trim(state),
                    county = stringr::str_trim(county))
  }
  if (geo_large == "county") {
    di_vars <- di_data %>%
      dplyr::mutate(oid = paste(.$STATEFP, .$COUNTYFP, sep = ""),
                    state = stringr::str_trim(state),
                    county = stringr::str_trim(county))
  }
  if (geo_large == "state") {
    di_vars <- di_data %>%
      dplyr::mutate(oid = .$STATEFP,
                    state = stringr::str_trim(state))
  }

  # Count of racial/ethnic subgroup populations
  ## Count of racial/ethnic comparison subgroup population
  if (length(in_subgroup == 1)) {
    di_vars <- di_vars %>%
      dplyr::mutate(subgroup = .[ , in_subgroup])
  } else {
    di_vars <- di_vars %>%
      dplyr::mutate(subgroup = rowSums(.[ , in_subgroup]))
  }
  ## Count of racial/ethnic reference subgroup population
  if (length(in_subgroup_ref == 1)) {
    di_vars <- di_vars %>%
      dplyr::mutate(subgroup_ref = .[ , in_subgroup_ref])
  } else {
    di_vars <- di_vars %>%
      dplyr::mutate(subgroup_ref = rowSums(.[ , in_subgroup_ref]))
  }

  # Compute DI
  ## From Duncan & Duncan (1955) https://doi.org/10.2307/2088328
  ## D_{jt} = 1/2 \sum_{i=1}^{k} | \frac{x_{ijt}}{X_{jt}}-\frac{y_{ijt}}{Y_{jt}}|
  ## Where for k smaller geographies:
  ## D_{jt} denotes the DI of larger geography j at time t
  ## x_{ijt} denotes the racial/ethnic subgroup population of smaller geography i within larger geography j at time t
  ## X_{jt} denotes the racial/ethnic subgroup population of larger geography j at time t
  ## y_{ijt} denotes the racial/ethnic referent subgroup population of smaller geography i within larger geography j at time t
  ## Y_{jt} denotes the racial/ethnic referent subgroup population of larger geography j at time t

  ## Compute
  DItmp <- di_vars %>%
    split(., f = list(di_vars$oid)) %>%
    lapply(., FUN = di_fun) %>%
    utils::stack(.) %>%
    dplyr::mutate(DI = values,
                  oid = ind) %>%
    dplyr::select(DI, oid)

  # Warning for missingness of census characteristics
  missingYN <- di_vars[ , c(in_subgroup, in_subgroup_ref)]
  names(missingYN) <- out_names
  missingYN <- missingYN %>%
    tidyr::pivot_longer(cols = dplyr::everything(),
                        names_to = "variable",
                        values_to = "val") %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(total = dplyr::n(),
                     n_missing = sum(is.na(val)),
                     percent_missing = paste0(round(mean(is.na(val)) * 100, 2), " %"))

  if (quiet == FALSE) {
    # Warning for missing census data
    if (sum(missingYN$n_missing) > 0) {
      message("Warning: Missing census data")
    } else {
      returnValue(missingYN)
    }
  }

  # Format output
  if (geo_large == "state") {
    di <- merge(di_vars, DItmp) %>%
      dplyr::select(oid, state, DI) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, state, DI) %>%
      .[.$GEOID != "NANA", ]
  }
  if (geo_large == "county") {
    di <- merge(di_vars, DItmp) %>%
      dplyr::select(oid, state, county, DI) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, state, county, DI) %>%
      .[.$GEOID != "NANA", ]
  }
  if (geo_large == "tract") {
    di <- merge(di_vars, DItmp) %>%
      dplyr::select(oid, state, county, tract, DI) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, state, county, tract, DI) %>%
      .[.$GEOID != "NANA", ]
  }

  di <- di %>%
    dplyr::arrange(GEOID) %>%
    dplyr::as_tibble()
  
  di_data <- di_data %>%
    dplyr::arrange(GEOID) %>%
    dplyr::as_tibble() 

  out <- list(di = di,
              di_data = di_data,
              missing = missingYN)
  
  return(out)
}
