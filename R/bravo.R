#' Educational Isolation Index based on Bravo et al. (2021)
#' 
#' Compute the Educational Isolation Index (Bravo) values for selected educational attainment category(ies).
#'
#' @param geo Character string specifying the geography of the data either census tracts \code{geo = "tract"} (the default) or counties \code{geo = "county"}.
#' @param year Numeric. The year to compute the estimate. The default is 2020, and the years 2009 onward are currently available.
#' @param subgroup Character string specifying the educational attainment category(ies). See Details for available choices.
#' @param quiet Logical. If TRUE, will display messages about potential missing census information. The default is FALSE.
#' @param ... Arguments passed to \code{\link[tidycensus]{get_acs}} to select state, county, and other arguments for census characteristics
#'
#' @details This function will compute the Educational Isolation Index (EI) of U.S. census tracts or counties for a specified geographical extent (e.g., the entire U.S. or a single state) based on Bravo et al. (2021) \doi{10.3390/ijerph18179384} who originally designed the metric for the educational isolation of individual without a college degree. This function provides the computation of EI for any of the U.S. Census Bureau educational attainment levels.
#' 
#' The function uses the \code{\link[tidycensus]{get_acs}} function to obtain U.S. Census Bureau 5-year American Community Survey characteristics used for the geospatial computation. The yearly estimates are available for 2009 onward when ACS-5 data are available but are available from other U.S. Census Bureau surveys. The five educational attainment levels (U.S. Census Bureau definitions) are:
#' \itemize{
#'  \item{B06009_002: }{Less than high school graduate "LtHS"}
#'  \item{B06009_003: }{High school graduate (includes equivalency) "HSGiE"}
#'  \item{B06009_004: }{Some college or associate's degree "SCoAD"}
#'  \item{B06009_005: }{Bachelor's degree "BD"}
#'  \item{B06009_006: }{Graduate or professional degree "GoPD"}
#' }
#' 
#' Use the internal \code{state} and \code{county} arguments within the \code{\link[tidycensus]{get_acs}} function to specify geographic extent of the data output. NOTE: Current version does not correct for edge effects (e.g., census geographies along the specified spatial extent border, coastline, or U.S.-Mexico / U.S.-Canada border) may have few neighboring census geographies, and EI values in these census geographies may be unstable. A stop-gap solution for the former source of edge effect is to compute the EI for neighboring census geographies (i.e., the states bordering a study area of interest) and then use the estimates of the study area of interest.
#' 
#' A census geography (and its neighbors) that has nearly all of its population with the specified educational attainment category (e.g., a Bachelor's degree or more) will have an EI value close to 1. In contrast, a census geography (and its neighbors) that is nearly none of its population with the specified educational attainment category (e.g., less than a Bachelor's degree) will have an EI value close to 0.
#' 
#' @return An object of class 'list'. This is a named list with the following components:
#' 
#' \describe{
#' \item{\code{ei}}{An object of class 'tbl' for the GEOID, name, EI, and raw census values of specified census geographies.}
#' \item{\code{missing}}{An object of class 'tbl' of the count and proportion of missingness for each census variable used to compute the EI.}
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
#'   bravo(geo = "tract", state = "GA", year = 2020, subgroup = c("LtHS", "HSGiE"))
#'   
#'   # County-level metric (2020)
#'   bravo(geo = "county", state = "GA", year = 2020, subgroup = c("LtHS", "HSGiE"))
#'   
#' }
#' 
bravo <- function(geo = "tract", year = 2020, subgroup, quiet = FALSE, ...) {
  
  # Check arguments
  match.arg(geo, choices = c("county", "tract"))
  stopifnot(is.numeric(year), year >= 2009) # all variables available 2009 onward
  match.arg(subgroup, several.ok = TRUE,
            choices = c("LtHS", "HSGiE", "SCoAD", "BD", "GoPD"))
  
  # Select census variables
  vars <- c(TotalPop = "B06009_001",
            LtHS = "B06009_002",
            HSGiE = "B06009_003",
            SCoAD = "B06009_004",
            BD = "B06009_005",
            GoPD = "B06009_006")
  
  selected_vars <- vars[c("TotalPop", subgroup)]
  out_names <- names(selected_vars) # save for output
  prefix <- "subgroup"
  suffix <- seq(1:length(subgroup))
  names(selected_vars) <- c("TotalPop", paste(prefix, suffix, sep = ""))
  in_names <- paste(names(selected_vars), "E", sep = "")
  
  # Acquire EI variables and sf geometries
  ei_vars <- suppressMessages(suppressWarnings(tidycensus::get_acs(geography = geo,
                                                                   year = year, 
                                                                   output = "wide",
                                                                   variables = selected_vars, 
                                                                   geometry = TRUE, ...)))
  
  
  if (geo == "tract") {
    ei_vars <- ei_vars %>%
      tidyr::separate(NAME, into = c("tract", "county", "state"), sep = ",") %>%
      dplyr::mutate(tract = gsub("[^0-9\\.]","", tract))
  } else {
    ei_vars <- ei_vars %>% tidyr::separate(NAME, into = c("county", "state"), sep = ",") 
  }
  
  ei_vars <- ei_vars %>% 
    dplyr::mutate(subgroup = rowSums(sf::st_drop_geometry(ei_vars[ , in_names[-1]])))
  
  # Compute EI
  ## From Bravo et al. (2021) https://doi.org/10.3390/ijerph18179384
  ## EI_{im} = (Sigma_{j∈∂_{i}} w_{ij} * T_{jm}) / (Sigma_{j∈∂_{i}} w_{ij} * T_{j})
  ## Where:
  ## ∂_{i} denotes the set of index units i and its neighbors
  ## Given M mutually exclusive subgroups of educational attainment categories, m indexes the subgroups of M
  ## T_{i} denotes the total population in region i (TotalPop)
  ## T_{im} denotes the population of the selected subgroup(s) (subgroup1, ...)
  ## w_{ij} denotes a nXn first-order adjacency matrix, where n is the number of census geometries in the study area
  ### and the entries of w_{ij} are set to 1 if a boundary is shared by region i and region j and zero otherwise
  ### Entries of the main diagonal (since i∈∂_{i}, w_{ij} = w_{ii} when j = i) of w_{ij} are set to 1.5
  ### such that the weight of the index unit, i, is larger than the weights assigned to adjacent tracts
  
  ## Geospatial adjacency matrix (wij)
  tmp <- sf::st_intersects(sf::st_geometry(ei_vars), sparse = TRUE)
  names(tmp) <- as.character(seq_len(nrow(ei_vars)))
  tmpL <- length(tmp)
  tmpcounts <- unlist(Map(length, tmp))
  tmpi <- rep(1:tmpL, tmpcounts)
  tmpj <- unlist(tmp)
  wij <- Matrix::sparseMatrix(i = tmpi, j = tmpj, x = 1, dims = c(tmpL, tmpL))
  diag(wij) <- 1.5
  
  ## Compute
  ei_vars <- sf::st_drop_geometry(ei_vars) # drop geometries (can join back later)
  EIim <- list()
  for (i in 1:dim(wij)[1]){
    EIim[[i]] <- sum(as.matrix(wij[i, ])*ei_vars[ , "subgroup"]) / sum(as.matrix(wij[i, ])*ei_vars[, "TotalPopE"])
  }
  ei_vars$EI <- unlist(EIim)
  
  # Warning for missingness of census characteristics
  missingYN <- ei_vars[ , in_names]
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
    ei <- ei_vars %>%
      dplyr::select(c("GEOID",
                      "state",
                      "county",
                      "tract",
                      "EI",
                      in_names))
    names(ei) <- c("GEOID", "state", "county", "tract", "EI", out_names)
  } else {
    ei <- ei_vars %>%
      dplyr::select(c("GEOID",
                      "state",
                      "county",
                      "EI",
                      in_names))
    names(ei) <- c("GEOID", "state", "county", "EI", out_names)
  }
  
  ei <- ei %>%
    dplyr::mutate(state = stringr::str_trim(state),
                  county = stringr::str_trim(county)) %>%
    dplyr::arrange(GEOID) %>%
    dplyr::as_tibble() 
  
  out <- list(ei = ei,
              missing = missingYN)
  
  return(out)
}
