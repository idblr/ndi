#' Racial Isolation Index based on Anthopolos et al. (2011) 
#' 
#' Compute the Racial Isolation Index (Anthopolos) values for a selected subgroup(s).
#'
#' @param geo Character string specifying the geography of the data either census tracts \code{geo = "tract"} (the default) or counties \code{geo = "county"}.
#' @param year Numeric. The year to compute the estimate. The default is 2020, and the years 2009 onward are currently available.
#' @param subgroup Character string specifying the racial/ethnic subgroup(s). See Details for available choices.
#' @param quiet Logical. If TRUE, will display messages about potential missing census information. The default is FALSE.
#' @param ... Arguments passed to \code{\link[tidycensus]{get_acs}} to select state, county, and other arguments for census characteristics
#'
#' @details This function will compute the Racial Isolation Index (RI) of U.S. census tracts or counties for a specified geographical extent (e.g., the entire U.S. or a single state) based on Anthopolos et al. (2011) \doi{10.1016/j.sste.2011.06.002} who originally designed the metric for the racial isolation of non-Hispanic Black individuals. This function provides the computation of RI for any of the U.S. Census Bureau race/ethnicity subgroups (including Hispanic and non-Hispanic individuals).
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
#' Use the internal \code{state} and \code{county} arguments within the \code{\link[tidycensus]{get_acs}} function to specify geographic extent of the data output. NOTE: Current version does not correct for edge effects (e.g., census geographies along the specified spatial extent border, coastline, or U.S.-Mexico / U.S.-Canada border) may have few neighboring census geographies, and RI values in these census geographies may be unstable. A stop-gap solution for the former source of edge effect is to compute the RI for neighboring census geographies (i.e., the states bordering a study area of interest) and then use the estimates of the study area of interest.
#' 
#' A census geography (and its neighbors) that has nearly all of its population who identify with the specified race/ethnicity subgroup(s) (e.g., non-Hispanic or Latino, Black or African American alone) will have an RI value close to 1. In contrast, a census geography (and its neighbors) that is nearly none of its population who identify with the specified race/ethnicity subgroup(s) (e.g., not non-Hispanic or Latino, Black or African American alone) will have an RI value close to 0.
#' 
#' @return An object of class 'list'. This is a named list with the following components:
#' 
#' \describe{
#' \item{\code{ri}}{An object of class 'tbl' for the GEOID, name, RI, and raw census values of specified census geographies.}
#' \item{\code{missing}}{An object of class 'tbl' of the count and proportion of missingness for each census variable used to compute the RI.}
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
#'   # Tract-level metric (2020)
#'   anthopolos(geo = "tract", state = "GA", year = 2020, subgroup = c("NHoLB", "HoLB"))
#'   
#'   # County-level metric (2020)
#'   anthopolos(geo = "county", state = "GA", year = 2020, subgroup = c("NHoLB", "HoLB"))
#'   
#' }
#' 
anthopolos <- function(geo = "tract", year = 2020, subgroup, quiet = FALSE, ...) {
  
  # Check arguments
  match.arg(geo, choices = c("county", "tract"))
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
  prefix <- "subgroup"
  suffix <- seq(1:length(subgroup))
  names(selected_vars) <- c("TotalPop", paste(prefix, suffix, sep = ""))
  in_names <- paste(names(selected_vars), "E", sep = "")
  
  # Acquire RI variables and sf geometries
  ri_vars <- suppressMessages(suppressWarnings(tidycensus::get_acs(geography = geo,
                                                                   year = year, 
                                                                   output = "wide",
                                                                   variables = selected_vars, 
                                                                   geometry = TRUE, ...)))
  
  
  if (geo == "tract") {
    ri_vars <- ri_vars %>%
      tidyr::separate(NAME, into = c("tract", "county", "state"), sep = ",") %>%
      dplyr::mutate(tract = gsub("[^0-9\\.]","", tract))
  } else {
    ri_vars <- ri_vars %>% tidyr::separate(NAME, into = c("county", "state"), sep = ",") 
  }
  
  ri_vars <- ri_vars %>% 
    dplyr::mutate(subgroup = rowSums(sf::st_drop_geometry(ri_vars[ , in_names[-1]])))
  
  # Compute RI
  ## From Anthopolos et al. (2011) https://doi.org/10.1016/j.sste.2011.06.002
  ## RI_{im} = (Sigma_{j∈∂_{i}} w_{ij} * T_{jm}) / (Sigma_{j∈∂_{i}} w_{ij} * T_{j})
  ## Where:
  ## ∂_{i} denotes the set of index units i and its neighbors
  ## Given M mutually exclusive racial/ethnic subgroups, m indexes the subgroups of M
  ## T_{i} denotes the total population in region i (TotalPop)
  ## T_{im} denotes the population of the selected subgroup(s) (subgroup1, ...)
  ## w_{ij} denotes a nXn first-order adjacency matrix, where n is the number of census geometries in the study area
  ### and the entries of w_{ij} are set to 1 if a boundary is shared by region i and region j and zero otherwise
  ### Entries of the main diagonal (since i∈∂_{i}, w_{ij} = w_{ii} when j = i) of w_{ij} are set to 1.5
  ### such that the weight of the index unit, i, is larger than the weights assigned to adjacent tracts
  
  ## Geospatial adjacency matrix (wij)
  tmp <- sf::st_intersects(sf::st_geometry(ri_vars), sparse = TRUE)
  names(tmp) <- as.character(seq_len(nrow(ri_vars)))
  tmpL <- length(tmp)
  tmpcounts <- unlist(Map(length, tmp))
  tmpi <- rep(1:tmpL, tmpcounts)
  tmpj <- unlist(tmp)
  wij <- Matrix::sparseMatrix(i = tmpi, j = tmpj, x = 1, dims = c(tmpL, tmpL))
  diag(wij) <- 1.5
  
  ## Compute
  ri_vars <- sf::st_drop_geometry(ri_vars) # drop geometries (can join back later)
  RIim <- list()
  for (i in 1:dim(wij)[1]){
    RIim[[i]] <- sum(as.matrix(wij[i, ])*ri_vars[ , "subgroup"]) / sum(as.matrix(wij[i, ])*ri_vars[, "TotalPopE"])
  }
  ri_vars$RI <- unlist(RIim)
  
  # Warning for missingness of census characteristics
  missingYN <- ri_vars[ , in_names]
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
    if (nrow(missingYN) != 0) {
      message("Warning: Missing census data")
    } else {
      returnValue(missingYN)
    }
  }
  
  # Format output
  if (geo == "tract") {
    ri <- ri_vars %>%
      dplyr::select(c("GEOID",
                      "state",
                      "county",
                      "tract",
                      "RI",
                      in_names))
    names(ri) <- c("GEOID", "state", "county", "tract", "RI", out_names)
  } else {
    ri <- ri_vars %>%
      dplyr::select(c("GEOID",
                      "state",
                      "county",
                      "RI",
                      in_names))
    names(ri) <- c("GEOID", "state", "county", "RI", out_names)
  }
  
  ri <- ri %>%
    dplyr::mutate(state = stringr::str_trim(state),
                  county = stringr::str_trim(county)) %>%
    dplyr::arrange(GEOID) %>%
    dplyr::as_tibble() 
  
  out <- list(ri = ri,
              missing = missingYN)
  
  return(out)
}
