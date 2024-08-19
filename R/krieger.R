#' Index of Concentration at the Extremes based on Feldman et al. (2015) and Krieger et al. (2016)
#'
#' Compute the aspatial Index of Concentration at the Extremes (Krieger).
#'
#' @param geo Character string specifying the geography of the data either census tracts \code{geo = 'tract'} (the default) or counties \code{geo = 'county'}.
#' @param year Numeric. The year to compute the estimate. The default is 2020, and the years 2009 onward are currently available.
#' @param quiet Logical. If TRUE, will display messages about potential missing census information. The default is FALSE.
#' @param ... Arguments passed to \code{\link[tidycensus]{get_acs}} to select state, county, and other arguments for census characteristics
#'
#' @details This function will compute three aspatial Index of Concentration at the Extremes (*ICE*) of U.S. census tracts or counties for a specified geographical extent (e.g., entire U.S. or a single state) based on Feldman et al. (2015) \doi{10.1136/jech-2015-205728} and Krieger et al. (2016) \doi{10.2105/AJPH.2015.302955}. The authors expanded the metric designed by Massey in a chapter of Booth & Crouter (2001) \doi{10.4324/9781410600141} who initially designed the metric for residential segregation. This function computes five *ICE* metrics:
#'
#' \itemize{
#' \item **Income**: 80th income percentile vs. 20th income percentile
#' \item **Education**: less than high school vs. four-year college degree or more
#' \item **Race/Ethnicity**: white non-Hispanic vs. black non-Hispanic
#' \item **Income and race/ethnicity combined**: white non-Hispanic in 80th income percentile vs. black alone (including Hispanic) in 20th income percentile
#' \item **Income and race/ethnicity combined**: white non-Hispanic in 80th income percentile vs. white non-Hispanic in 20th income percentile
#' }
#'
#' The function uses the \code{\link[tidycensus]{get_acs}} function to obtain U.S. Census Bureau 5-year American Community Survey characteristics used for the geospatial computation. The yearly estimates are available for 2009 onward when ACS-5 data are available but are available from other U.S. Census Bureau surveys. The ACS-5 groups used in the computation of the five *ICE* metrics are:
#' \itemize{
#'  \item **B03002**: HISPANIC OR LATINO ORIGIN BY RACE
#'  \item **B15002**: SEX BY EDUCATIONAL ATTAINMENT FOR THE POPULATION 25 YEARS AND OVER
#'  \item **B19001**: HOUSEHOLD INCOME IN THE PAST 12 MONTHS (IN 20XX INFLATION-ADJUSTED DOLLARS)
#'  \item **B19001B**: HOUSEHOLD INCOME IN THE PAST 12 MONTHS (IN 20XX INFLATION-ADJUSTED DOLLARS) (BLACK OR AFRICAN AMERICAN ALONE HOUSEHOLDER)
#'  \item **B19001H**: HOUSEHOLD INCOME IN THE PAST 12 MONTHS (IN 20XX INFLATION-ADJUSTED DOLLARS) (WHITE ALONE, NOT HISPANIC OR LATINO HOUSEHOLDER)
#' }
#'
#' Use the internal \code{state} and \code{county} arguments within the \code{\link[tidycensus]{get_acs}} function to specify geographic extent of the data output.
#'
#' *ICE* metrics can range in value from -1 (most deprived) to 1 (most privileged). A value of 0 can thus represent two possibilities: (1) none of the residents are in the most privileged or most deprived categories, or (2) an equal number of persons are in the most privileged and most deprived categories, and in both cases indicates that the area is not dominated by extreme concentrations of either of the two groups.
#'
#' @return An object of class 'list'. This is a named list with the following components:
#'
#' \describe{
#' \item{\code{ice}}{An object of class 'tbl' for the GEOID, name, *ICE* metrics, and raw census values of specified census geographies.}
#' \item{\code{missing}}{An object of class 'tbl' of the count and proportion of missingness for each census variable used to compute the *ICE* metrics.}
#' }
#'
#' @import dplyr
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
#'   krieger(geo = 'tract', state = 'GA', year = 2020)
#'
#'   # County-level metric (2020)
#'   krieger(geo = 'county', state = 'GA', year = 2020)
#'
#' }
#'
krieger <- function(geo = 'tract',
                    year = 2020,
                    quiet = FALSE,
                    ...) {
  
    # Check arguments
    match.arg(geo, choices = c('county', 'tract'))
    stopifnot(is.numeric(year), year >= 2009) # all variables available 2009 onward
    
    # Select census variables
    vars <- c(
      TotalPopi = 'B19001_001',
      TotalPopedu = 'B15002_001',
      TotalPopre = 'B03002_001',
      U10i = 'B19001_002',
      B1015i = 'B19001_003',
      B1520i = 'B19001_004',
      B2025i = 'B19001_005',
      B2530i = 'B19001_006',
      B100125i = 'B19001_014',
      B125150i = 'B19001_015',
      B150200i = 'B19001_016',
      O200i = 'B19001_017',
      O25MNSC = 'B15002_003',
      O25FNSC = 'B15002_020',
      O25MNt4G = 'B15002_004',
      O25FNt4G = 'B15002_021',
      O25M5t6G = 'B15002_005',
      O25F5t6G = 'B15002_022',
      O25M7t8G = 'B15002_006',
      O25F7t8G = 'B15002_023',
      O25M9G = 'B15002_007',
      O25F9G = 'B15002_024',
      O25M10G = 'B15002_008',
      O25F10G = 'B15002_025',
      O25M11G = 'B15002_009',
      O25F11G = 'B15002_026',
      O25M12GND = 'B15002_010',
      O25F12GND = 'B15002_027',
      O25MBD = 'B15002_015',
      O25FBD = 'B15002_032',
      O25MMD = 'B15002_016',
      O25FMD = 'B15002_033',
      O25MPSD = 'B15002_017',
      O25FPSD = 'B15002_034',
      O25MDD = 'B15002_018',
      O25FDD = 'B15002_035',
      NHoLW = 'B03002_003',
      NHoLB = 'B03002_004',
      U10nhw = 'B19001H_002',
      B1015nhw = 'B19001H_003',
      B1520nhw = 'B19001H_004',
      B2025nhw = 'B19001H_005',
      B2530nhw = 'B19001H_006',
      B100125nhw = 'B19001H_014',
      B125150nhw = 'B19001H_015',
      B150200nhw = 'B19001H_016',
      O200nhw = 'B19001H_017',
      U10bih = 'B19001B_002',
      B1015bih = 'B19001B_003',
      B1520bih = 'B19001B_004',
      B2025bih = 'B19001B_005',
      B2530bih = 'B19001B_006'
    )
    
    # Acquire ICE variables
    ice_data <- suppressMessages(suppressWarnings(
      tidycensus::get_acs(
        geography = geo,
        year = year,
        output = 'wide',
        variables = vars,
        ...
      )
    ))
      
    
    if (geo == 'tract') {
      ice_data <- ice_data %>%
        tidyr::separate(NAME, into = c('tract', 'county', 'state'), sep = ',') %>%
        dplyr::mutate(tract = gsub('[^0-9\\.]', '', tract))
    } else {
      ice_data <- ice_data %>% 
        tidyr::separate(NAME, into = c('county', 'state'), sep = ',')
    }
    
    ice_data <- ice_data %>%
      dplyr::mutate(
        TotalPop_inc = TotalPopiE,
        TotalPop_edu = TotalPopeduE,
        TotalPop_re = TotalPopreE,
        U10i = U10iE,
        B1015i = B1015iE,
        B1520i = B1520iE,
        B2025i = B2025iE,
        B2530i = B2530iE,
        B100125i = B100125iE,
        B125150i = B125150iE,
        B150200i = B150200iE,
        O200i = O200iE,
        O25MNSC = O25MNSCE,
        O25FNSC = O25FNSCE,
        O25MNt4G = O25MNt4GE,
        O25FNt4G = O25FNt4GE,
        O25M5t6G = O25M5t6GE,
        O25F5t6G = O25F5t6GE,
        O25M7t8G = O25M7t8GE,
        O25F7t8G = O25F7t8GE,
        O25M9G = O25M9GE,
        O25F9G = O25F9GE,
        O25M10G = O25M10GE,
        O25F10G = O25F10GE,
        O25M11G = O25M11GE,
        O25F11G = O25F11GE,
        O25M12GND = O25M12GNDE,
        O25F12GND = O25F12GNDE,
        O25MBD = O25MBDE,
        O25FBD = O25FBDE,
        O25MMD = O25MMDE,
        O25FMD = O25FMDE,
        O25MPSD = O25MPSDE,
        O25FPSD = O25FPSDE,
        O25MDD = O25MDDE,
        O25FDD = O25FDDE,
        NHoLW = NHoLWE,
        NHoLB = NHoLBE,
        U10nhw = U10nhwE,
        B1015nhw = B1015nhwE,
        B1520nhw = B1520nhwE,
        B2025nhw = B2025nhwE,
        B2530nhw = B2530nhwE,
        B100125nhw = B100125nhwE,
        B125150nhw = B125150nhwE,
        B150200nhw = B150200nhwE,
        O200nhw = O200nhwE,
        U10bih = U10bihE,
        B1015bih = B1015bihE,
        B1520bih = B1520bihE,
        B2025bih = B2025bihE,
        B2530bih = B2530bihE
      )
    
    # Sum educational attainment categories
    # A_{edu} = Less than high school / 12 year / GED
    # P_{edu} = Four-year college degree or more
    ice_data <- ice_data %>%
      dplyr::mutate(
        A_edu = O25MBD + O25FBD + O25MMD + O25FMD + O25MPSD + O25FPSD + O25MDD + O25FDD,
        P_edu = O25MNSC + O25FNSC + O25MNt4G + O25FNt4G + O25M5t6G + O25F5t6G + O25M7t8G + 
          O25F7t8G + O25M9G + O25F9G + O25M10G + O25F10G + O25M11G + O25F11G + O25M12GND + 
          O25F12GND
      )
    
    # Sum income percentile counts
    ## A_income (A_{inc}) is the 80th income percentile
    ## P_income (P_{inc}) is the 20th income percentile
    ## Add 'Total, $25,000 to $34,999' for years 2016 and after
    ## Remove 'Total, $100,000 to $124,999' for years 2016 and after
    ## According to U.S. Census Bureau Table A-4a
    ## 'Selected Measures of Household Income Dispersion: 1967 to 2020'
    if (year < 2016) {
      ice_data <- ice_data %>%
        dplyr::mutate(
          A_inc = B100125i + B125150i + B150200i + O200i,
          P_inc = U10i + B1015i + B1520i + B2025i,
          A_wbinc = B100125nhw + B125150nhw + B150200nhw + O200nhw,
          P_wbinc = U10bih + B1015bih + B1520bih + B2025bih,
          A_wpcinc = B100125nhw + B125150nhw + B150200nhw + O200nhw,
          P_wpcinc = U10nhw + B1015nhw + B1520nhw + B2025nhw
        )
    } else {
      ice_data <- ice_data %>%
        dplyr::mutate(
          A_inc = B125150i + B150200i + O200i,
          P_inc = U10i + B1015i + B1520i + B2025i + B2530i,
          A_wbinc = B125150nhw + B150200nhw + O200nhw,
          P_wbinc = U10bih + B1015bih + B1520bih + B2025bih + B2530bih,
          A_wpcinc = B125150nhw + B150200nhw + O200nhw,
          P_wpcinc = U10nhw + B1015nhw + B1520nhw + B2025nhw + B2530nhw
        )
    }
    
    # Compute ICEs
    ## From Kreiger et al. (2016) https://doi.org/10.2105%2FAJPH.2015.302955
    ## ICE_{i} = (A_{i} - P_{i}) / T_{i}
    ## Where:
    ## A_{i} denotes the count within the lowest extreme (e.g., households in 20th income percentile)
    ## P_{i} denotes the count within the highest extreme (e.g., households in 80th income percentile)
    ## T_{i} denotes the total population in region i (TotalPop)
    
    ice_data <- ice_data %>%
      dplyr::mutate(
        ICE_inc = (A_inc - P_inc) / TotalPop_inc,
        ICE_edu = (A_edu - P_edu) / TotalPop_edu,
        ICE_rewb = (NHoLW - NHoLB) / TotalPop_re,
        ICE_wbinc = (A_wbinc - P_wbinc) / TotalPop_inc,
        ICE_wpcinc = (A_wpcinc - P_wpcinc) / TotalPop_inc
      )
    
    # Warning for missingness of census characteristics
    missingYN <- ice_data %>%
      dplyr::select(
        U10i,
        B1015i,
        B1520i,
        B2025i,
        B2530i,
        B100125i,
        B125150i,
        B150200i,
        O200i,
        O25MNSC,
        O25FNSC,
        O25MNt4G,
        O25FNt4G,
        O25M5t6G,
        O25F5t6G,
        O25M7t8G,
        O25F7t8G,
        O25M9G,
        O25F9G,
        O25M10G,
        O25F10G,
        O25M11G,
        O25F11G,
        O25M12GND,
        O25F12GND,
        O25MBD,
        O25FBD,
        O25MMD,
        O25FMD,
        O25MPSD,
        O25FPSD,
        O25MDD,
        O25FDD,
        NHoLW,
        NHoLB,
        U10nhw,
        B1015nhw,
        B1520nhw,
        B2025nhw,
        B2530nhw,
        B100125nhw,
        B125150nhw,
        B150200nhw,
        O200nhw,
        U10bih,
        B1015bih,
        B1520bih,
        B2025bih,
        B2530bih,
        TotalPop_inc,
        TotalPop_edu,
        TotalPop_re
      ) %>%
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
    if (geo == 'tract') {
      ice <- ice_data %>%
        dplyr::select(
          GEOID,
          state,
          county,
          tract,
          ICE_inc,
          ICE_edu,
          ICE_rewb,
          ICE_wbinc,
          ICE_wpcinc,
          U10i,
          B1015i,
          B1520i,
          B2025i,
          B2530i,
          B100125i,
          B125150i,
          B150200i,
          O200i,
          O25MNSC,
          O25FNSC,
          O25MNt4G,
          O25FNt4G,
          O25M5t6G,
          O25F5t6G,
          O25M7t8G,
          O25F7t8G,
          O25M9G,
          O25F9G,
          O25M10G,
          O25F10G,
          O25M11G,
          O25F11G,
          O25M12GND,
          O25F12GND,
          O25MBD,
          O25FBD,
          O25MMD,
          O25FMD,
          O25MPSD,
          O25FPSD,
          O25MDD,
          O25FDD,
          NHoLW,
          NHoLB,
          U10nhw,
          B1015nhw,
          B1520nhw,
          B2025nhw,
          B2530nhw,
          B100125nhw,
          B125150nhw,
          B150200nhw,
          O200nhw,
          U10bih,
          B1015bih,
          B1520bih,
          B2025bih,
          B2530bih,
          TotalPop_inc,
          TotalPop_edu,
          TotalPop_re
        )
    } else {
      ice <- ice_data %>%
        dplyr::select(
          GEOID,
          state,
          county,
          ICE_inc,
          ICE_edu,
          ICE_rewb,
          ICE_wbinc,
          ICE_wpcinc,
          U10i,
          B1015i,
          B1520i,
          B2025i,
          B2530i,
          B100125i,
          B125150i,
          B150200i,
          O200i,
          O25MNSC,
          O25FNSC,
          O25MNt4G,
          O25FNt4G,
          O25M5t6G,
          O25F5t6G,
          O25M7t8G,
          O25F7t8G,
          O25M9G,
          O25F9G,
          O25M10G,
          O25F10G,
          O25M11G,
          O25F11G,
          O25M12GND,
          O25F12GND,
          O25MBD,
          O25FBD,
          O25MMD,
          O25FMD,
          O25MPSD,
          O25FPSD,
          O25MDD,
          O25FDD,
          NHoLW,
          NHoLB,
          U10nhw,
          B1015nhw,
          B1520nhw,
          B2025nhw,
          B2530nhw,
          B100125nhw,
          B125150nhw,
          B150200nhw,
          O200nhw,
          U10bih,
          B1015bih,
          B1520bih,
          B2025bih,
          B2530bih,
          TotalPop_inc,
          TotalPop_edu,
          TotalPop_re
        )
    }
    
    ice <- ice %>%
      dplyr::mutate(
        state = stringr::str_trim(state),
        county = stringr::str_trim(county)
      ) %>%
      dplyr::arrange(GEOID) %>%
      dplyr::as_tibble()
    
    out <- list(ice = ice, missing = missingYN)
    
    return(out)
  }
