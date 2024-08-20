#' Neighborhood Deprivation Index based on Messer et al. (2006) 
#' 
#' Compute the aspatial Neighborhood Deprivation Index (Messer).
#'
#' @param geo Character string specifying the geography of the data either census tracts \code{geo = 'tract'} (the default) or counties \code{geo = 'county'}.
#' @param year Numeric. The year to compute the estimate. The default is 2020, and the years 2010 onward are currently available.
#' @param imp Logical. If TRUE, will impute missing census characteristics within the internal \code{\link[psych]{principal}}. If FALSE (the default), will not impute. 
#' @param quiet Logical. If TRUE, will display messages about potential missing census information and the proportion of variance explained by principal component analysis. The default is FALSE.
#' @param round_output Logical. If TRUE, will round the output of raw census and \emph{NDI} values from the \code{\link[tidycensus]{get_acs}} at one and four significant digits, respectively. The default is FALSE.
#' @param df Optional. Pass a pre-formatted \code{'dataframe'} or \code{'tibble'} with the desired variables through the function. Bypasses the data obtained by \code{\link[tidycensus]{get_acs}}. The default is NULL. See Details below.
#' @param ... Arguments passed to \code{\link[tidycensus]{get_acs}} to select state, county, and other arguments for census characteristics
#'
#' @details This function will compute the aspatial Neighborhood Deprivation Index (\emph{NDI}) of U.S. census tracts or counties for a specified geographical referent (e.g., US-standardized) based on Messer et al. (2006) \doi{10.1007/s11524-006-9094-x}.
#' 
#' The function uses the \code{\link[tidycensus]{get_acs}} function to obtain U.S. Census Bureau 5-year American Community Survey characteristics used for computation involving a principal component analysis with the \code{\link[psych]{principal}} function. The yearly estimates are available for 2010 and after when all census characteristics became available. The eight characteristics are:
#' \itemize{
#'  \item \strong{OCC (C24030)}: percent males in management, science, and arts occupation
#'  \item \strong{CWD (B25014)}: percent of crowded housing
#'  \item \strong{POV (B17017)}: percent of households in poverty
#'  \item \strong{FHH (B25115)}: percent of female headed households with dependents
#'  \item \strong{PUB (B19058)}: percent of households on public assistance
#'  \item \strong{U30 (B19001)}: percent of households earning <$30,000 per year
#'  \item \strong{EDU (B06009)}: percent earning less than a high school education
#'  \item \strong{EMP (B23025)}: percent unemployed (2011 onward)
#'  \item \strong{EMP (B23001)}: percent unemployed (2010 only)
#' }
#' 
#' Use the internal \code{state} and \code{county} arguments within the \code{\link[tidycensus]{get_acs}} function to specify the referent for standardizing the \emph{NDI} (Messer) values. For example, if all U.S. states are specified for the \code{state} argument, then the output would be a U.S.-standardized index.
#' 
#' The continuous \emph{NDI} (Messer) values are z-transformed, i.e., 'standardized,' and the categorical \emph{NDI} (Messer) values are quartiles of the standardized continuous \emph{NDI} (Messer) values. 
#' 
#' Check if the proportion of variance explained by the first principal component is high (more than 0.5).
#' 
#' Users can bypass \code{\link[tidycensus]{get_acs}} by specifying a pre-formatted data frame or tibble using the \code{df} argument. This function will compute an index using the first component of a principal component analysis (PCA) with a Varimax rotation (the default for \code{\link[psych]{principal}}) and only one factor (note: PCA set-up not unspecified in Messer et al. (2006)). The recommended structure of the data frame or tibble is an ID (e.g., GEOID) in the first feature (column), followed by the variables of interest (in any order) and no additional information (e.g., omit state or county names from the \code{df} argument input).
#' 
#' @return An object of class 'list'. This is a named list with the following components:
#' 
#' \describe{
#' \item{\code{ndi}}{An object of class 'tbl' for the GEOID, name, \emph{NDI} (standardized), \emph{NDI} (quartile), and raw census values of specified census geographies.}
#' \item{\code{pca}}{An object of class 'principal', returns the output of \code{\link[psych]{principal}} used to compute the \emph{NDI} values.}
#' \item{\code{missing}}{An object of class 'tbl' of the count and proportion of missingness for each census variable used to compute \emph{NDI}.}
#' }
#' 
#' @import dplyr
#' @importFrom psych principal
#' @importFrom stats quantile
#' @importFrom stringr str_trim
#' @importFrom tidycensus get_acs
#' @importFrom tidyr pivot_longer separate
#' @export
#' 
#' @seealso \code{\link[tidycensus]{get_acs}} for additional arguments for geographic referent selection (i.e., \code{state} and \code{county}).
#'
#' @examples
#' 
#' messer(df = DCtracts2020[ , c(1, 3:10)])
#' 
#' \dontrun{
#' # Wrapped in \dontrun{} because these examples require a Census API key.
#'   
#'   # Tract-level metric (2020)
#'   messer(geo = 'tract', state = 'GA', year = 2020)
#'
#'   # Impute NDI for tracts (2020) with missing census information (median values)
#'   messer(state = 'tract', 'GA', year = 2020, imp = TRUE)
#'   
#'   # County-level metric (2020)
#'   messer(geo = 'county', state = 'GA', year = 2020)
#'   
#' }
#' 
messer <- function(geo = 'tract',
                   year = 2020,
                   imp = FALSE,
                   quiet = FALSE,
                   round_output = FALSE,
                   df = NULL,
                   ...) {
  
  # Check arguments
  if (!is.null(df) &
      !inherits(df, c('tbl_df', 'tbl', 'data.frame'))) {
    stop("df' must be class 'data.frame' or 'tbl'")
  }
  
  if (is.null(df)) {
    # Check additional arguments
    match.arg(geo, choices = c('county', 'tract'))
    stopifnot(is.numeric(year), year >= 2010) # all variables available 2010 onward
    
    # Select census variables
    vars <-
      c(
        PctMenMgmtBusScArti_num1 = 'C24030_018',
        PctMenMgmtBusScArti_num2 = 'C24030_019',
        PctMenMgmtBusScArti_den = 'C24030_002',
        PctCrwdHH_num1 = 'B25014_005',
        PctCrwdHH_num2 = 'B25014_006',
        PctCrwdHH_num3 = 'B25014_007',
        PctCrwdHH_num4 = 'B25014_011',
        PctCrwdHH_num5 = 'B25014_012',
        PctCrwdHH_num6 = 'B25014_013',
        PctCrwdHH_den = 'B25014_001',
        PctHHPov_num = 'B17017_002',
        PctHHPov_den = 'B17017_001',
        PctFemHeadKids_num1 = 'B25115_012',
        PctFemHeadKids_num2 = 'B25115_025',
        PctFemHeadKids_den = 'B25115_001',
        PctPubAsst_num = 'B19058_002',
        PctPubAsst_den = 'B19058_001',
        PctHHUnder30K_num1 = 'B19001_002',
        PctHHUnder30K_num2 = 'B19001_003',
        PctHHUnder30K_num3 = 'B19001_004',
        PctHHUnder30K_num4 = 'B19001_005',
        PctHHUnder30K_num5 = 'B19001_006',
        PctHHUnder30K_den = 'B19001_001',
        PctEducLessThanHS_num = 'B06009_002',
        PctEducLessThanHS_den = 'B06009_001',
        PctUnemp_num = 'B23025_005',
        PctUnemp_den = 'B23025_003'
      )
    
    if (year == 2010) {
      # Select census variables
      vars <- c(
        vars[-c(26, 27)],
        PctUnemp_den = 'B23001_001',
        PctUnemp_1619M = 'B23001_008',
        PctUnemp_2021M = 'B23001_015',
        PctUnemp_2224M = 'B23001_022',
        PctUnemp_2529M = 'B23001_029',
        PctUnemp_3034M = 'B23001_036',
        PctUnemp_3544M = 'B23001_043',
        PctUnemp_4554M = 'B23001_050',
        PctUnemp_5559M = 'B23001_057',
        PctUnemp_6061M = 'B23001_064',
        PctUnemp_6264M = 'B23001_071',
        PctUnemp_6569M = 'B23001_076',
        PctUnemp_7074M = 'B23001_081',
        PctUnemp_75upM = 'B23001_086',
        PctUnemp_1619F = 'B23001_094',
        PctUnemp_2021F = 'B23001_101',
        PctUnemp_2224F = 'B23001_108',
        PctUnemp_2529F = 'B23001_115',
        PctUnemp_3034F = 'B23001_122',
        PctUnemp_3544F = 'B23001_129',
        PctUnemp_4554F = 'B23001_136',
        PctUnemp_5559F = 'B23001_143',
        PctUnemp_6061F = 'B23001_150',
        PctUnemp_6264F = 'B23001_157',
        PctUnemp_6569F = 'B23001_162',
        PctUnemp_7074F = 'B23001_167',
        PctUnemp_75upF = 'B23001_172'
      )
      
      # Acquire NDI variables
      ndi_data <- suppressMessages(suppressWarnings(
        tidycensus::get_acs(
          geography = geo,
          year = year,
          output = 'wide',
          variables = vars,
          ...
        )
      ))
      
      if (geo == 'tract') {
        ndi_data <- ndi_data %>%
          tidyr::separate(NAME, into = c('tract', 'county', 'state'), sep = ',') %>%
          dplyr::mutate(tract = gsub('[^0-9\\.]', '', tract))
      } else {
        ndi_data <-
          ndi_data %>% tidyr::separate(NAME, into = c('county', 'state'), sep = ',')
      }
      
      ndi_data <- ndi_data %>%
        dplyr::mutate(
          OCC = (PctMenMgmtBusScArti_num1E + PctMenMgmtBusScArti_num2E) / PctMenMgmtBusScArti_denE,
          CWD = (
            PctCrwdHH_num1E + PctCrwdHH_num2E + PctCrwdHH_num3E + PctCrwdHH_num4E + 
              PctCrwdHH_num5E + PctCrwdHH_num6E
          ) / PctCrwdHH_denE,
          POV = PctHHPov_numE / PctHHPov_denE,
          FHH = (PctFemHeadKids_num1E + PctFemHeadKids_num2E) / PctFemHeadKids_denE,
          PUB = PctPubAsst_numE / PctPubAsst_denE,
          U30 = (
            PctHHUnder30K_num1E + PctHHUnder30K_num2E + PctHHUnder30K_num3E + PctHHUnder30K_num4E + 
              PctHHUnder30K_num5E
          ) / PctHHUnder30K_denE,
          EDU = PctEducLessThanHS_numE / PctEducLessThanHS_denE,
          EMP = (
            PctUnemp_1619ME + PctUnemp_2021ME +
              PctUnemp_2224ME + PctUnemp_2529ME +
              PctUnemp_4554ME + PctUnemp_5559ME +
              PctUnemp_6061ME + PctUnemp_6264ME +
              PctUnemp_6569ME + PctUnemp_7074ME +
              PctUnemp_75upME + PctUnemp_1619FE +
              PctUnemp_2021FE + PctUnemp_2224FE +
              PctUnemp_2529FE + PctUnemp_4554FE +
              PctUnemp_5559FE + PctUnemp_6061FE +
              PctUnemp_6264FE + PctUnemp_6569FE +
              PctUnemp_7074FE + PctUnemp_75upME
          ) / PctUnemp_denE
        )
    } else {
      # Acquire NDI variables
      ndi_data <- suppressMessages(suppressWarnings(
        tidycensus::get_acs(
          geography = geo,
          year = year,
          output = 'wide',
          variables = vars,
          ...
        )
      ))
        
      if (geo == 'tract') {
        ndi_data <- ndi_data %>%
          tidyr::separate(NAME, into = c('tract', 'county', 'state'), sep = ',') %>%
          dplyr::mutate(tract = gsub('[^0-9\\.]', '', tract))
      } else {
        ndi_data <-
          ndi_data %>% tidyr::separate(NAME, into = c('county', 'state'), sep = ',')
      }
      
      ndi_data <- ndi_data %>%
        dplyr::mutate(
          OCC = (PctMenMgmtBusScArti_num1E + PctMenMgmtBusScArti_num2E) / PctMenMgmtBusScArti_denE,
          CWD = (
            PctCrwdHH_num1E + PctCrwdHH_num2E + PctCrwdHH_num3E + PctCrwdHH_num4E + 
              PctCrwdHH_num5E + PctCrwdHH_num6E
          ) / PctCrwdHH_denE,
          POV = PctHHPov_numE / PctHHPov_denE,
          FHH = (PctFemHeadKids_num1E + PctFemHeadKids_num2E) / PctFemHeadKids_denE,
          PUB = PctPubAsst_numE / PctPubAsst_denE,
          U30 = (
            PctHHUnder30K_num1E + PctHHUnder30K_num2E + PctHHUnder30K_num3E + PctHHUnder30K_num4E + 
              PctHHUnder30K_num5E
          ) / PctHHUnder30K_denE,
          EDU = PctEducLessThanHS_numE / PctEducLessThanHS_denE,
          EMP = PctUnemp_numE / PctUnemp_denE
        )
    }
    
    # Generate NDI
    ndi_data_pca <- ndi_data %>%
      dplyr::select(OCC, CWD, POV, FHH, PUB, U30, EDU, EMP)
  } else {
    # If inputing pre-formatted data:
    ndi_data <- dplyr::as_tibble(df)
    # omit the first feature (column) typically an ID (e.g., GEOID or FIPS)
    ndi_data_pca <- df[,-1] 
  }
  
  # Replace infinite values as zero (typically because denominator is zero)
  ndi_data_pca <- do.call(
    data.frame,
    lapply(ndi_data_pca, function(x) replace(x, is.infinite(x), 0))
  )
  
  # Run principal component analysis
  pca <- psych::principal(
    ndi_data_pca,
    nfactors = 1,
    n.obs = nrow(ndi_data_pca),
    covar = FALSE,
    scores = TRUE,
    missing = imp
  )
  
  # Warning for missingness of census characteristics
  missingYN <- ndi_data_pca %>%
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
    
    # Warning for proportion of variance explained by PC1
    if (pca$Vaccounted[2] < 0.50) {
      message('Warning: The proportion of variance explained by PC1 is less than 0.50.')
    }
  }
  
  # NDI quartiles
  NDIQuart <- data.frame(PC1 = pca$scores) %>%
    dplyr::mutate(
      NDI = PC1 / pca$value[1] ^ 2,
      NDIQuart = cut(
        NDI,
        breaks = stats::quantile(NDI, probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE),
        labels = c(
          '1-Least deprivation',
          '2-BelowAvg deprivation',
          '3-AboveAvg deprivation',
          '4-Most deprivation'
        ),
        include.lowest = TRUE
      ),
      NDIQuart = factor(
        replace(as.character(NDIQuart), is.na(NDIQuart), '9-NDI not avail'),
        c(levels(NDIQuart), '9-NDI not avail')
      )
    ) %>%
    dplyr::select(NDI, NDIQuart)
  
  if (is.null(df)) {
    # Format output
    if (round_output == TRUE) {
      ndi <- cbind(ndi_data, NDIQuart) %>%
        dplyr::mutate(
          OCC = round(OCC, digits = 1),
          CWD = round(CWD, digits = 1),
          POV = round(POV, digits = 1),
          FHH = round(FHH, digits = 1),
          PUB = round(PUB, digits = 1),
          U30 = round(U30, digits = 1),
          EDU = round(EDU, digits = 1),
          EMP = round(EMP, digits = 1),
          NDI = round(NDI, digits = 4)
        )
    } else {
      ndi <- cbind(ndi_data, NDIQuart)
    }
    
    if (geo == 'tract') {
      ndi <- ndi %>%
        dplyr::select(
          GEOID,
          state,
          county,
          tract,
          NDI,
          NDIQuart,
          OCC,
          CWD,
          POV,
          FHH,
          PUB,
          U30,
          EDU,
          EMP
        )
    } else {
      ndi <- ndi %>%
        dplyr::select(
          GEOID,
          state,
          county,
          NDI,
          NDIQuart,
          OCC,
          CWD,
          POV,
          FHH,
          PUB,
          U30,
          EDU,
          EMP
        )
    }
    
    ndi <- ndi %>%
      dplyr::mutate(
        state = stringr::str_trim(state),
        county = stringr::str_trim(county)
      ) %>%
      dplyr::arrange(GEOID) %>%
      dplyr::as_tibble()
    
  } else {
    ndi <- cbind(df[, 1], NDIQuart, df[, 2:ncol(df)])
    ndi <- dplyr::as_tibble(ndi[order(ndi[, 1]),])
  }
  
  out <- list(ndi = ndi, pca = pca, missing = missingYN)
  
  return(out)
}
