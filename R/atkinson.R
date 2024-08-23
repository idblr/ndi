#' Atkinson Index based on Atkinson (1970)
#'
#' Compute the aspatial Atkinson Index of income or selected racial/ethnic subgroup(s) and U.S. geographies.
#'
#' @param geo_large Character string specifying the larger geographical unit of the data. The default is counties \code{geo_large = 'county'}.
#' @param geo_small Character string specifying the smaller geographical unit of the data. The default is census tracts \code{geo_large = 'tract'}.
#' @param year Numeric. The year to compute the estimate. The default is 2020, and the years 2009 onward are currently available.
#' @param subgroup Character string specifying the income or racial/ethnic subgroup(s) as the comparison population. See Details for available choices.
#' @param epsilon Numerical. Shape parameter that denotes the aversion to inequality. Value must be between 0 and 1.0 (the default is 0.5).
#' @param holder Logical. If TRUE, will compute index using the Hölder mean. If FALSE, will not compute with the Hölder mean. The default is FALSE.
#' @param omit_NAs Logical. If FALSE, will compute index for a larger geographical unit only if all of its smaller geographical units have values. The default is TRUE.
#' @param quiet Logical. If TRUE, will display messages about potential missing census information. The default is FALSE.
#' @param ... Arguments passed to \code{\link[tidycensus]{get_acs}} to select state, county, and other arguments for census characteristics
#'
#' @details This function will compute the aspatial Atkinson Index (\emph{A}) of income or selected racial/ethnic subgroups and U.S. geographies for a specified geographical extent (e.g., the entire U.S. or a single state) based on Atkinson (1970) \doi{10.1016/0022-0531(70)90039-6}. This function provides the computation of \emph{A} for median household income and any of the U.S. Census Bureau race/ethnicity subgroups (including Hispanic and non-Hispanic individuals).
#'
#' The function uses the \code{\link[tidycensus]{get_acs}} function to obtain U.S. Census Bureau 5-year American Community Survey characteristics used for the aspatial computation. The yearly estimates are available for 2009 onward when ACS-5 data are available (2010 onward for \code{geo_large = 'cbsa'} and 2011 onward for \code{geo_large = 'csa'} or \code{geo_large = 'metro'}) but may be available from other U.S. Census Bureau surveys. When \code{subgroup = 'MedHHInc'}, the metric will be computed for median household income ('B19013_001'). The twenty racial/ethnic subgroups (U.S. Census Bureau definitions) are:
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
#' \emph{A} is a measure of the evenness of residential inequality (e.g., racial/ethnic segregation) when comparing smaller geographical areas to larger ones within which the smaller geographical areas are located. \emph{A} can range in value from 0 to 1 with smaller values indicating lower levels of inequality (e.g., less segregation).
#'
#' The \code{epsilon} argument that determines how to weight the increments to inequality contributed by different proportions of the Lorenz curve. A user must explicitly decide how heavily to weight smaller geographical units at different points on the Lorenz curve (i.e., whether the index should take greater account of differences among areas of over- or under-representation). The \code{epsilon} argument must have values between 0 and 1.0. For \code{0 <= epsilon < 0.5} or less 'inequality-averse,' smaller geographical units with a subgroup proportion smaller than the subgroup proportion of the larger geographical unit contribute more to inequality ('over-representation'). For \code{0.5 < epsilon <= 1.0} or more 'inequality-averse,' smaller geographical units with a subgroup proportion larger than the subgroup proportion of the larger geographical unit contribute more to inequality ('under-representation'). If \code{epsilon = 0.5} (the default), units of over- and under-representation contribute equally to the index. See Section 2.3 of Saint-Jacques et al. (2020) \doi{10.48550/arXiv.2002.05819} for one method to select \code{epsilon}.
#'
#' Larger geographies available include state \code{geo_large = 'state'}, county \code{geo_large = 'county'}, census tract \code{geo_large = 'tract'}, Core Based Statistical Area \code{geo_large = 'cbsa'}, Combined Statistical Area \code{geo_large = 'csa'}, and Metropolitan Division \code{geo_large = 'metro'} levels. Smaller geographies available include, county \code{geo_small = 'county'}, census tract \code{geo_small = 'tract'}, and census block group \code{geo_small = 'block group'} levels. If a larger geographical area is comprised of only one smaller geographical area (e.g., a U.S county contains only one census tract), then the \emph{A} value returned is NA. If the larger geographical unit is Combined Based Statistical Areas \code{geo_large = 'csa'} or Core Based Statistical Areas \code{geo_large = 'cbsa'}, only the smaller geographical units completely within a larger geographical unit are considered in the \emph{A} computation (see internal \code{\link[sf]{st_within}} function for more information) and recommend specifying all states within which the interested larger geographical unit are located using the internal \code{state} argument to ensure all appropriate smaller geographical units are included in the \emph{A} computation.
#' 
#' @return An object of class 'list'. This is a named list with the following components:
#'
#' \describe{
#' \item{\code{a}}{An object of class 'tbl' for the GEOID, name, and \emph{A} at specified larger census geographies.}
#' \item{\code{a_data}}{An object of class 'tbl' for the raw census values at specified smaller census geographies.}
#' \item{\code{missing}}{An object of class 'tbl' of the count and proportion of missingness for each census variable used to compute \emph{A}.}
#' }
#'
#' @import dplyr
#' @importFrom sf st_drop_geometry st_within
#' @importFrom stats na.omit
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
#'   # Atkinson Index (a measure of the evenness) of Black populations
#'   ## of census tracts within counties within Georgia, U.S.A., counties (2020)
#'   atkinson(
#'     geo_large = 'county',
#'     geo_small = 'tract',
#'     state = 'GA',
#'     year = 2020,
#'     subgroup = c('NHoLB', 'HoLB')
#'    )
#'
#' }
#'
atkinson <- function(geo_large = 'county',
                     geo_small = 'tract',
                     year = 2020,
                     subgroup,
                     epsilon = 0.5,
                     holder = FALSE,
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
        'HoLTReSOR',
        'MedHHInc'
      )
    )
    stopifnot(is.numeric(epsilon), epsilon >= 0 , epsilon <= 1) # values between 0 and 1
    
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
    in_subgroup <- paste0(subgroup, 'E')
    
    # Acquire A variables and sf geometries
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
    
    # Grouping IDs for A computation
    if (geo_large == 'state') {
      out_dat <- out_dat %>%
        dplyr::mutate(
          oid = STATEFP,
          state = stringr::str_trim(state)
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
    if (geo_large == 'tract') {
      out_dat <- out_dat %>%
        dplyr::mutate(
          oid = paste0(STATEFP, COUNTYFP, TRACTCE),
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
      stopifnot(is.numeric(year), year >= 2011) # Metro Divisions only available 2011 onward
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
      stopifnot(is.numeric(year), year >= 2011) # CSAs only available 2011 onward
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

    # Count of racial/ethnic subgroup populations
    ## Count of racial/ethnic subgroup population
    if (length(in_subgroup) == 1) {
      out_dat <- out_dat %>%
        dplyr::mutate(subgroup = .[, in_subgroup])
    } else {
      out_dat <- out_dat %>%
        dplyr::mutate(subgroup = rowSums(.[, in_subgroup]))
    }
    
    # Compute A
    ## From Atkinson (1970) https://doi.org/10.1016/0022-0531(70)90039-6
    ## A_{\epsilon}(x_{1},...,x_{n}) = \begin{Bmatrix}
    ## 1 - (\frac{1}{n}\sum_{i=1}^{n}x_{i}^{1-\epsilon})^{1/(1-\epsilon)}/(\frac{1}{n}\sum_{i=1}^{n}x_{i}) & \mathrm{if\:} \epsilon \neq 1 \\
    ## 1 - (\prod_{i=1}^{n}x_{i})^{1/n}/(\frac{1}{n}\sum_{i=1}^{n}x_{i}) & \mathrm{if\:} \epsilon = 1 \\
    ## \end{Bmatrix}
    ## Where the Atkinson index (A) is defined for a population subgroup count (x) of a given smaller geographical unit (i) for n smaller geographical units
    ## and an inequality-aversion parameter (epsilon)
    ## If denoting the Hölder mean (based on `Atkinson()` function in 'DescTools' package) by
    ## M_{p}(x_{1},...,x_{n}) = \begin{Bmatrix}
    ## (\frac{1}{n}\sum_{i=1}^{n}x_{i}^{p})^{1/p} & \mathrm{if\:} p \neq 0 \\
    ## (\prod_{i=1}^{n}x_{i})^{1/n} & \mathrm{if\:} p = 0 \\
    ## \end{Bmatrix}
    ## then A is
    ## A_{\epsilon}(x_{1},...,x_{n}) = 1 - \frac{M_{1-\epsilon}(x_{1},...,x_{n})}{M_{1}(x_{1},...,x_{n})}
    
    ## Compute
    out_tmp <- out_dat %>%
      split(., f = list(out_dat$oid)) %>%
      lapply(., FUN = a_fun, epsilon = epsilon, omit_NAs = omit_NAs, holder = holder) %>%
      utils::stack(.) %>%
      dplyr::mutate(
        A = values,
        oid = ind
      ) %>%
      dplyr::select(A, oid)
    
    # Warning for missingness of census characteristics
    missingYN <- as.data.frame(out_dat[, c('TotalPopE', in_subgroup)])
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
        dplyr::select(oid, state, A) %>%
        unique(.) %>%
        dplyr::mutate(GEOID = oid) %>%
        dplyr::select(GEOID, state, A) %>%
        .[.$GEOID != 'NANA',]
    }
    if (geo_large == 'county') {
      out <- out_dat %>%
        dplyr::left_join(out_tmp, by = dplyr::join_by(oid)) %>%
        dplyr::select(oid, state, county, A) %>%
        unique(.) %>%
        dplyr::mutate(GEOID = oid) %>%
        dplyr::select(GEOID, state, county, A) %>%
        .[.$GEOID != 'NANA',]
    }
    if (geo_large == 'tract') {
      out <- out_dat %>%
        dplyr::left_join(out_tmp, by = dplyr::join_by(oid)) %>%
        dplyr::select(oid, state, county, tract, A) %>%
        unique(.) %>%
        dplyr::mutate(GEOID = oid) %>%
        dplyr::select(GEOID, state, county, tract, A) %>%
        .[.$GEOID != 'NANA',]
    }
    if (geo_large == 'cbsa') {
      out <- out_dat %>%
        dplyr::left_join(out_tmp, by = dplyr::join_by(oid)) %>%
        dplyr::select(oid, cbsa, A) %>%
        unique(.) %>%
        dplyr::mutate(GEOID = oid) %>%
        dplyr::select(GEOID, cbsa, A) %>%
        .[.$GEOID != 'NANA', ] %>%
        dplyr::distinct(GEOID, .keep_all = TRUE) %>%
        dplyr::filter(stats::complete.cases(.))
    }
    if (geo_large == 'csa') {
      out <- out_dat %>%
        dplyr::left_join(out_tmp, by = dplyr::join_by(oid)) %>%
        dplyr::select(oid, csa, A) %>%
        unique(.) %>%
        dplyr::mutate(GEOID = oid) %>%
        dplyr::select(GEOID, csa, A) %>%
        .[.$GEOID != 'NANA', ] %>%
        dplyr::distinct(GEOID, .keep_all = TRUE) %>%
        dplyr::filter(stats::complete.cases(.))
    }
    if (geo_large == 'metro') {
      out <- out_dat %>%
        dplyr::left_join(out_tmp, by = dplyr::join_by(oid)) %>%
        dplyr::select(oid, metro, A) %>%
        unique(.) %>%
        dplyr::mutate(GEOID = oid) %>%
        dplyr::select(GEOID, metro, A) %>%
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
    
    out <- list(a = out, a_data = out_dat, missing = missingYN)
    
    return(out)
  }
