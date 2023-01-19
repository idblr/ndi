#' Correlation Ratio based on Bell (1954) and White (1986) 
#' 
#' Compute the aspatial Correlation Ratio (White) of a selected racial/ethnic subgroup(s) and U.S. geographies.
#'
#' @param geo_large Character string specifying the larger geographical unit of the data. The default is counties \code{geo_large = "county"}.
#' @param geo_small Character string specifying the smaller geographical unit of the data. The default is census tracts \code{geo_large = "tract"}.
#' @param year Numeric. The year to compute the estimate. The default is 2020, and the years 2009 onward are currently available.
#' @param subgroup Character string specifying the racial/ethnic subgroup(s). See Details for available choices.
#' @param omit_NAs Logical. If FALSE, will compute index for a larger geographical unit only if all of its smaller geographical units have values. The default is TRUE.
#' @param quiet Logical. If TRUE, will display messages about potential missing census information. The default is FALSE.
#' @param ... Arguments passed to \code{\link[tidycensus]{get_acs}} to select state, county, and other arguments for census characteristics
#'
#' @details This function will compute the aspatial Correlation Ratio (V or \eqn{Eta^{2}}{Eta^2}) of selected racial/ethnic subgroups and U.S. geographies for a specified geographical extent (e.g., the entire U.S. or a single state) based on Bell (1954) \doi{10.2307/2574118} and White (1986) \doi{10.2307/3644339}. This function provides the computation of V for any of the U.S. Census Bureau race/ethnicity subgroups (including Hispanic and non-Hispanic individuals).
#' 
#' The function uses the \code{\link[tidycensus]{get_acs}} function to obtain U.S. Census Bureau 5-year American Community Survey characteristics used for the aspatial computation. The yearly estimates are available for 2009 onward when ACS-5 data are available but are available from other U.S. Census Bureau surveys. The twenty racial/ethnic subgroups (U.S. Census Bureau definitions) are:
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
#' V removes the asymmetry from the Isolation Index (Bell) by controlling for the effect of population composition. The Isolation Index (Bell) is some measure of the probability that a member of one subgroup(s) will meet or interact with a member of another subgroup(s) with higher values signifying higher probability of interaction (less isolation). V can range in value from 0 to 1.
#' 
#' Larger geographies available include state \code{geo_large = "state"}, county \code{geo_large = "county"}, and census tract \code{geo_large = "tract"} levels. Smaller geographies available include, county \code{geo_small = "county"}, census tract \code{geo_small = "tract"}, and census block group \code{geo_small = "block group"} levels. If a larger geographical area is comprised of only one smaller geographical area (e.g., a U.S county contains only one census tract), then the DI value returned is NA.
#' 
#' @return An object of class 'list'. This is a named list with the following components:
#' 
#' \describe{
#' \item{\code{v}}{An object of class 'tbl' for the GEOID, name, and V at specified larger census geographies.}
#' \item{\code{v_data}}{An object of class 'tbl' for the raw census values at specified smaller census geographies.}
#' \item{\code{missing}}{An object of class 'tbl' of the count and proportion of missingness for each census variable used to compute V}
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
#'   # Isolation of non-Hispanic Black and non-Hispanic white populations
#'   ## of census tracts within Georgia, U.S.A., counties (2020)
#'   white(geo_large = "county", geo_small = "tract", state = "GA", year = 2020,
#'         subgroup = "NHoLB")
#'   
#' }
#' 
white <- function(geo_large = "county", geo_small = "tract", year = 2020, subgroup, omit_NAs = TRUE, quiet = FALSE, ...) {
  
  # Check arguments
  match.arg(geo_large, choices = c("state", "county", "tract"))
  match.arg(geo_small, choices = c("county", "tract", "block group"))
  stopifnot(is.numeric(year), year >= 2009) # all variables available 2009 onward
  match.arg(subgroup, several.ok = TRUE,
            choices = c("NHoL", "NHoLW", "NHoLB", "NHoLAIAN", "NHoLA", "NHoLNHOPI",
                        "NHoLSOR", "NHoLTOMR", "NHoLTRiSOR", "NHoLTReSOR",
                        "HoL", "HoLW", "HoLB", "HoLAIAN", "HoLA", "HoLNHOPI",
                        "HoLSOR", "HoLTOMR", "HoLTRiSOR", "HoLTReSOR"))
  
  # Select census variables
  vars <- c(TotalPop = "B03002_001",
            NHoL = "B03002_002",
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
  
  selected_vars <- vars[c("TotalPop", subgroup)]
  out_names <- names(selected_vars) # save for output
  in_subgroup <- paste(subgroup, "E", sep = "")
  
  # Acquire V variables and sf geometries
  v_data <- suppressMessages(suppressWarnings(tidycensus::get_acs(geography = geo_small,
                                                                  year = year, 
                                                                  output = "wide",
                                                                  variables = selected_vars, 
                                                                  geometry = TRUE,
                                                                  keep_geo_vars = TRUE, ...)))
  
  # Format output
  if (geo_small == "county") {
    v_data <- sf::st_drop_geometry(v_data) %>%
      tidyr::separate(NAME.y, into = c("county", "state"), sep = ",")
  }
  if (geo_small == "tract") {
    v_data <- sf::st_drop_geometry(v_data) %>%
      tidyr::separate(NAME.y, into = c("tract", "county", "state"), sep = ",") %>%
      dplyr::mutate(tract = gsub("[^0-9\\.]", "", tract))
  } 
  if (geo_small == "block group") {
    v_data <- sf::st_drop_geometry(v_data) %>%
      tidyr::separate(NAME.y, into = c("block.group", "tract", "county", "state"), sep = ",") %>%
      dplyr::mutate(tract = gsub("[^0-9\\.]", "", tract),
                    block.group = gsub("[^0-9\\.]", "", block.group))
  } 
  
  # Grouping IDs for R computation
  if (geo_large == "tract") {
    v_vars <- v_data %>%
      dplyr::mutate(oid = paste(.$STATEFP, .$COUNTYFP, .$TRACTCE, sep = ""),
                    state = stringr::str_trim(state),
                    county = stringr::str_trim(county))
  }
  if (geo_large == "county") {
    v_vars <- v_data %>%
      dplyr::mutate(oid = paste(.$STATEFP, .$COUNTYFP, sep = ""),
                    state = stringr::str_trim(state),
                    county = stringr::str_trim(county))
  }
  if (geo_large == "state") {
    v_vars <- v_data %>%
      dplyr::mutate(oid = .$STATEFP,
                    state = stringr::str_trim(state))
  }
  
  # Count of racial/ethnic subgroup populations
  ## Count of racial/ethnic comparison subgroup population
  if (length(in_subgroup) == 1) {
    v_vars <- v_vars %>%
      dplyr::mutate(subgroup = .[ , in_subgroup])
  } else {
    v_vars <- v_vars %>%
      dplyr::mutate(subgroup = rowSums(.[ , in_subgroup]))
  }
  
  # Compute V or \mathit{Eta}^{2}
  ## From White (1986) https://doi.org/10.2307/3644339
  ## V = \mathit{Eta}^2 = [(_{x}P_{x}^* - P) / (1 - P)]
  ## Where:
  ## _{x}P_{x}^* denotes the Isolation Index (Bell) of subgroup x
  ## P denotes the proportion of subgroup x of study (reference) area
  
  ## Compute
  Vtmp <- v_vars %>%
    split(., f = list(v_vars$oid)) %>%
    lapply(., FUN = v_fun, omit_NAs = omit_NAs) %>%
    utils::stack(.) %>%
    dplyr::mutate(V = values,
                  oid = ind) %>%
    dplyr::select(V, oid)
  
  # Warning for missingness of census characteristics
  missingYN <- v_vars[ , c("TotalPopE", in_subgroup)]
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
    }
  }
  
  # Format output
  if (geo_large == "state") {
    v <- merge(v_vars, Vtmp) %>%
      dplyr::select(oid, state, V) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, state, V) %>%
      .[.$GEOID != "NANA", ]
  }
  if (geo_large == "county") {
    v <- merge(v_vars, Vtmp) %>%
      dplyr::select(oid, state, county, V) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, state, county, V) %>%
      .[.$GEOID != "NANA", ]
  }
  if (geo_large == "tract") {
    v <- merge(v_vars, Vtmp) %>%
      dplyr::select(oid, state, county, tract, V) %>%
      unique(.) %>%
      dplyr::mutate(GEOID = oid) %>%
      dplyr::select(GEOID, state, county, tract, V) %>%
      .[.$GEOID != "NANA", ]
  }
  
  v <- v %>%
    dplyr::arrange(GEOID) %>%
    dplyr::as_tibble()
  
  v_data <- v_data %>%
    dplyr::arrange(GEOID) %>%
    dplyr::as_tibble() 
  
  out <- list(v = v,
              v_data = v_data,
              missing = missingYN)
  
  return(out)
}
