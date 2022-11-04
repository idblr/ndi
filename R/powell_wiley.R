#' Neighborhood Deprivation Index based on Andrews et al. (2020) and Slotman et al. (2022)
#' 
#' Compute the Neighborhood Deprivation Index (Powell-Wiley) values.
#'
#' @param geo Character string specifying the geography of the data either census tracts \code{geo = "tract"} (the default) or counties \code{geo = "county"}.
#' @param year Numeric. The year to compute the estimate. The default is 2020, and the years 2010 onward are currently available.
#' @param imp Logical. If TRUE, will impute missing census characteristics within the internal \code{\link[psych]{principal}} using median values of variables. If FALSE (the default), will not impute. 
#' @param quiet Logical. If TRUE, will display messages about potential missing census information, standardized Cronbach's alpha, and proportion of variance explained by principal component analysis. The default is FALSE.
#' @param round_output Logical. If TRUE, will round the output of raw census and NDI values from the \code{\link[tidycensus]{get_acs}} at one and four significant digits, respectively. The default is FALSE.
#' @param df Optional. Pass a pre-formatted \code{'dataframe'} or \code{'tibble'} with the desired variables through the function. Bypasses the data obtained by \code{\link[tidycensus]{get_acs}}. The default is NULL. See Details below.
#' @param ... Arguments passed to \code{\link[tidycensus]{get_acs}} to select state, county, and other arguments for census characteristics
#'
#' @details This function will compute the Neighborhood Deprivation Index (NDI) of U.S. census tracts or counties for a specified geographical referent (e.g., US-standardized) based on Andrews et al. (2020) \doi{10.1080/17445647.2020.1750066} and Slotman et al. (2022) \doi{10.1016/j.dib.2022.108002}.
#' 
#' The function uses the \code{\link[tidycensus]{get_acs}} function to obtain U.S. Census Bureau 5-year American Community Survey characteristics used for computation involving a factor analysis with the \code{\link[psych]{principal}} function. The yearly estimates are available in 2010 and after when all census characteristics became available. The thirteen characteristics chosen by Roux and Mair (2010) \doi{10.1111/j.1749-6632.2009.05333.x} are:
#' \itemize{
#'  \item{MedHHInc (5B19013): }{median household income (dollars)}
#'  \item{PctRecvIDR (B19054): }{percent of households receiving dividends, interest, or rental income}
#'  \item{PctPubAsst (B19058): }{percent of households receiving public assistance}
#'  \item{MedHomeVal (B25077): }{median home value (dollars)}
#'  \item{PctMgmtBusScArti (C24060): }{percent in a management, business, science, or arts occupation}
#'  \item{PctFemHeadKids (B11005): }{percent of households that are female headed with any children under 18 years}
#'  \item{PctOwnerOcc (DP04): }{percent of housing units that are owner occupied}
#'  \item{PctNoPhone (DP04): }{percent of households without a telephone}
#'  \item{PctNComPlm (DP04): }{percent of households without complete plumbing facilities}
#'  \item{PctEducHSPlus (S1501): }{percent with a high school degree or higher (population 25 years and over)}
#'  \item{PctEducBchPlus (S1501): }{percent with a college degree or higher (population 25 years and over)}
#'  \item{PctFamBelowPov (S1702): }{percent of families with incomes below the poverty level}
#'  \item{PctUnempl (S2301): }{percent unemployed}
#' }
#' 
#' Use the internal \code{state} and \code{county} arguments within the \code{\link[tidycensus]{get_acs}} function to specify the referent for standardizing the NDI (Powell-Wiley) values. For example, if all U.S. states are specified for the \code{state} argument, then the output would be a U.S.-standardized index. Please note: the NDI (Powell-Wiley) values will not exactly match (but will highly correlate with) those found in Andrews et al. (2020) \doi{10.1080/17445647.2020.1750066} and Slotman et al. (2022) \doi{10.1016/j.dib.2022.108002} because the two studies used a different statistical platform (i.e., SPSS and SAS, respectively) that intrinsically calculate the principal component analysis differently from R.
#' 
#' The categorical NDI (Powell-Wiley) values are population-weighted quintiles of the continuous NDI (Powell-Wiley) values. 
#' 
#' Check if the proportion of variance explained by the first principal component is high (more than 0.5).
#' 
#' Users can bypass \code{\link[tidycensus]{get_acs}} by specifying a pre-formatted data frame or tibble using the \code{df} argument. This function will compute an index using the first component of a principal component analysis (PCA) with a Promax (oblique) rotation and a minimum Eigenvalue of 1, omitting variables with absolute loading score < 0.4. The recommended structure of the data frame or tibble is an ID (e.g., GEOID) in the first feature (column), an estimate of the total population in the second feature (column), followed by the variables of interest (in any order) and no additional information (e.g., omit state or county names from the \code{df} argument input).
#' 
#' @return An object of class 'list'. This is a named list with the following components:
#' 
#' \describe{
#' \item{\code{ndi}}{An object of class 'tbl' for the GEOID, name, NDI continuous, NDI quintiles, and raw census values of specified census geographies.}
#' \item{\code{pca}}{An object of class 'principal', returns the output of \code{\link[psych]{principal}} used to compute the NDI values.}
#' \item{\code{missing}}{An object of class 'tbl' of the count and proportion of missingness for each census variable used to compute the NDI.}
#' \item{\code{cronbach}}{An object of class 'character' or 'numeric' for the results of the Cronbach's alpha calculation. If only one factor is computed, a message is returned. If more than one factor is computed, Cronbach's alpha is calculated and should check that it is >0.7 for respectable internal consistency between factors.}
#' }
#' 
#' @import dplyr 
#' @importFrom MASS ginv
#' @importFrom psych alpha principal
#' @importFrom stats complete.cases cor cov2cor loadings median promax quantile sd
#' @importFrom stringr str_trim
#' @importFrom tidycensus get_acs
#' @importFrom tidyr pivot_longer separate
#' @export
#' 
#' @seealso \code{\link[tidycensus]{get_acs}} for additional arguments for geographic referent selection (i.e., \code{state} and \code{county}).
#'
#' @examples
#' 
#' powell_wiley(df = DCtracts2020[ , -c(3:10)])
#' 
#' \dontrun{
#' # Wrapped in \dontrun{} because these examples require a Census API key.
#'   
#'   # Tract-level metric (2020)
#'   powell_wiley(geo = "tract", state = "GA", year = 2020)
#'
#'   # Impute NDI for tracts (2020) with missing census information (median values)
#'   powell_wiley(state = "tract", "GA", year = 2020, imp = TRUE)
#'   
#'   # County-level metric (2020)
#'   powell_wiley(geo = "county", state = "GA", year = 2020)
#'   
#' }
#' 
powell_wiley <- function(geo = "tract", year = 2020, imp = FALSE, quiet = FALSE, round_output = FALSE, df = NULL, ...) {
  
  # Check arguments
  if (!is.null(df) & !inherits(df, c("tbl_df", "tbl", "data.frame"))) { stop("'df' must be class 'data.frame' or 'tbl'") }
  
  if (is.null(df)) {
    
    # Check additional arguments
    match.arg(geo, choices = c("county", "tract"))
    stopifnot(is.numeric(year), year >= 2010) # all variables available 2010 onward
    
    # Select census variables
    vars <- c(MedHHInc = "B19013_001",
              PctRecvIDR_num = "B19054_002", PctRecvIDR_den = "B19054_001",
              PctPubAsst_num = "B19058_002", PctPubAsst_den = "B19058_001",
              MedHomeVal = "B25077_001",
              PctMgmtBusScArti_num = "C24060_002",  PctMgmtBusScArti_den = "C24060_001", 
              PctFemHeadKids_num1 = "B11005_007", PctFemHeadKids_num2 = "B11005_010",
              PctFemHeadKids_den = "B11005_001",
              PctOwnerOcc = "DP04_0046P",
              PctNoPhone = "DP04_0075P",
              PctNComPlmb = "DP04_0073P",
              PctEduc_num25upHS = "S1501_C01_009",
              PctEduc_num25upSC = "S1501_C01_010",
              PctEduc_num25upAD = "S1501_C01_011",
              PctEduc_num25upBD = "S1501_C01_012",
              PctEduc_num25upGD = "S1501_C01_013",
              PctEduc_den25up = "S1501_C01_006",
              PctFamBelowPov = "S1702_C02_001",
              PctUnempl = "S2301_C04_001",
              TotalPopulation = "B01001_001")
    
    # Updated census variable definition(s)
    if (year < 2015){ vars <- c(vars[-13], PctNoPhone = "DP04_0074P") }
    
    # Acquire NDI variables
    ndi_vars <- suppressMessages(suppressWarnings(tidycensus::get_acs(geography = geo,
                                                                      year = year,
                                                                      output = "wide",
                                                                      variables = vars, ...)))
    
    if (geo == "tract") {
      ndi_vars <- ndi_vars %>%
        tidyr::separate(NAME, into = c("tract", "county", "state"), sep = ",") %>%
        dplyr::mutate(tract = gsub("[^0-9\\.]","", tract))
    } else {
      ndi_vars <- ndi_vars %>% tidyr::separate(NAME, into = c("county", "state"), sep = ",") 
    }
    
    ndi_vars <- ndi_vars %>%
      dplyr::mutate(MedHHInc = MedHHIncE,
                    PctRecvIDR = PctRecvIDR_numE / PctRecvIDR_denE * 100,
                    PctPubAsst = PctPubAsst_numE / PctPubAsst_denE * 100,
                    MedHomeVal = MedHomeValE,
                    PctMgmtBusScArti = PctMgmtBusScArti_numE / PctMgmtBusScArti_denE * 100,
                    PctFemHeadKids = (PctFemHeadKids_num1E + PctFemHeadKids_num2E) / PctFemHeadKids_denE * 100,
                    PctOwnerOcc = PctOwnerOccE,
                    PctNoPhone = PctNoPhoneE,
                    PctNComPlmb = PctNComPlmbE,
                    PctEducHSPlus = (PctEduc_num25upHSE + PctEduc_num25upSCE + PctEduc_num25upADE +
                                       PctEduc_num25upBDE + PctEduc_num25upGDE) / PctEduc_den25upE * 100,
                    PctEducBchPlus = (PctEduc_num25upBDE + PctEduc_num25upGDE) / PctEduc_den25upE * 100,
                    PctFamBelowPov = PctFamBelowPovE,
                    PctUnempl = PctUnemplE,
                    TotalPop = TotalPopulationE) %>%
      # Log transform median household income and median home value
      # Reverse code percentages so that higher values represent more deprivation 
      # Round percentages to 1 decimal place
      dplyr::mutate(logMedHHInc = log(MedHHInc),
                    logMedHomeVal = log(MedHomeVal),
                    PctNoIDR = 100 - PctRecvIDR,
                    PctWorkClass = 100 - PctMgmtBusScArti,
                    PctNotOwnerOcc = 100 - PctOwnerOcc,
                    PctEducLTHS = 100 - PctEducHSPlus,
                    PctEducLTBch = 100 - PctEducBchPlus) %>%
      # Z-standardize the percentages
      dplyr::mutate(PctNoIDRZ = scale(PctNoIDR),
                    PctPubAsstZ = scale(PctPubAsst),
                    PctWorkClassZ = scale(PctWorkClass),
                    PctFemHeadKidsZ = scale(PctFemHeadKids),
                    PctNotOwnerOccZ = scale(PctNotOwnerOcc),
                    PctNoPhoneZ = scale(PctNoPhone),
                    PctNComPlmbZ = scale(PctNComPlmb),
                    PctEducLTHSZ = scale(PctEducLTHS),
                    PctEducLTBchZ = scale(PctEducLTBch),
                    PctFamBelowPovZ = scale(PctFamBelowPov),
                    PctUnemplZ = scale(PctUnempl))
    
    # generate NDI
    ndi_vars_pca <- ndi_vars %>%
      dplyr::select(logMedHHInc, PctNoIDRZ, PctPubAsstZ, logMedHomeVal, PctWorkClassZ,
                    PctFemHeadKidsZ, PctNotOwnerOccZ, PctNoPhoneZ, PctNComPlmbZ, PctEducLTHSZ,
                    PctEducLTBchZ, PctFamBelowPovZ, PctUnemplZ)
  } else {
    # If inputing pre-formatted data: 
    colnames(df)[1:2] <- c("GEOID", "TotalPop") # rename first and second features (columns) with name to match above 
    ndi_vars <- dplyr::as_tibble(df)
    ndi_vars_pca <- ndi_vars[ , -c(1:2)] # omits the first two features (columns) typically an ID (e.g., GEOID or FIPS) and TotalPop
  }
  # Run a factor analysis using Promax (oblique) rotation and a minimum Eigenvalue of 1
  nfa <- eigen(stats::cor(ndi_vars_pca, use = "complete.obs"))
  nfa <- sum(nfa$values > 1) # count of factors with a minimum Eigenvalue of 1
  fit <- psych::principal(ndi_vars_pca, 
                          nfactors = nfa,
                          rotate = "none")
  fit_rotate <- stats::promax(stats::loadings(fit), m = 3)
  
  # Calculate the factors using only variables with an absolute loading score > 0.4 for the first factor
  ## If number of factors > 2, use structure matrix, else pattern matrix
  if (nfa > 1) {
    P_mat <- matrix(stats::loadings(fit_rotate), nrow = 13, ncol = nfa)
    
    # Structure matrix (S_mat) from under-the-hood of the psych::principal() function
    rot.mat <- fit_rotate$rotmat # rotation matrix
    ui <- solve(rot.mat)
    Phi <- cov2cor(ui %*% t(ui)) # interfactor correlation
    S_mat <- P_mat %*% Phi # pattern matrix multiplied by interfactor correlation
    
  } else {
    P_mat <- matrix(fit_rotate, nrow = 13, ncol = 1)
    Phi <- 1
    S_mat <- P_mat
  }
  
  ## Variable correlation matrix (R_mat)
  R_mat <- as.matrix(cor(ndi_vars_pca[complete.cases(ndi_vars_pca), ]))
  
  ## standardized score coefficients or weight matrix (B_mat)
  B_mat <- solve(R_mat, S_mat)
  
  # Additional PCA Information
  fit_rotate$rotation <- "promax"
  fit_rotate$Phi <- Phi
  fit_rotate$Structure <- S_mat
  
  if (nfa > 1) {
    fit_rotate$communality <- rowSums(P_mat^2)
  } else {
    fit_rotate$communality <- P_mat^2
  }
  fit_rotate$uniqueness <- diag(R_mat) - fit_rotate$communality
  
  if (nfa > 1) {
    vx <- colSums(P_mat^2)
  } else {
    vx <- sum(P_mat^2)
  }
  
  vtotal <- sum(fit_rotate$communality + fit_rotate$uniqueness)
  vx <- diag(Phi %*% t(P_mat) %*% P_mat)
  names(vx) <- colnames(loadings)
  varex <- rbind(`SS loadings` = vx)
  varex <- rbind(varex, `Proportion Var` = vx/vtotal)
  if (nfa > 1) {
    varex <- rbind(varex, `Cumulative Var` = cumsum(vx/vtotal))
    varex <- rbind(varex, `Proportion Explained` = vx/sum(vx))
    varex <- rbind(varex, `Cumulative Proportion` = cumsum(vx/sum(vx)))
  }
  fit_rotate$Vaccounted <- varex
  
  if (imp == TRUE) {
    ndi_vars_scrs <- as.matrix(ndi_vars_pca)
    miss <- which(is.na(ndi_vars_scrs), arr.ind = TRUE)
    item.med <- apply(ndi_vars_scrs, 2, stats::median, na.rm = TRUE)
    ndi_vars_scrs[miss] <- item.med[miss[, 2]]
  } else {
    ndi_vars_scrs <- ndi_vars_pca
  }
  
  scrs <- as.matrix(scale(ndi_vars_scrs[complete.cases(ndi_vars_scrs), abs(S_mat[ , 1]) > 0.4 ])) %*% B_mat[abs(S_mat[ , 1]) > 0.4, 1]
  
  ndi_vars_NA <- ndi_vars[complete.cases(ndi_vars_scrs), ]
  ndi_vars_NA$NDI <- c(scrs)
  
  ndi_vars_NDI <- dplyr::left_join(ndi_vars[ , c("GEOID", "TotalPop")], ndi_vars_NA[ , c("GEOID", "NDI")], by = "GEOID", all.x = TRUE)
  
  # Calculate Cronbach's alpha correlation coefficient among the factors and verify values are above 0.7. 
  if (nfa == 1) { 
    crnbch <- "Only one factor with minimum Eigenvalue of 1. Cannot calculate Cronbach's alpha."
  } else {
    cronbach <- suppressMessages(psych::alpha(ndi_vars_pca[ , abs(S_mat[ , 1]) > 0.4 ], check.keys = TRUE, na.rm = TRUE, warnings = FALSE))
    crnbch <- cronbach$total$std.alpha
  }
  
  # Warning for missingness of census characteristics
  missingYN <- ndi_vars_pca %>%
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
    
    # Warning for Cronbach's alpha < 0.7
    if (cronbach$total$std.alpha < 0.7) {
      message("Warning: Cronbach's alpha correlation coefficient among the factors is less than 0.7.")
    }
    
    # Warning for proportion of variance explained by FA1
    if (fit_rotate$Vaccounted[2] < 0.50) {
      message("Warning: The proportion of variance explained by PC1 is less than 0.50.")
    }
  }
  
  # NDI quintiles weighted by tract population
  NDIQuint <- ndi_vars_NDI %>%
    dplyr::mutate(NDIQuint = cut(NDI*log(TotalPop),
                                 breaks = stats::quantile(NDI*log(TotalPop),
                                                          probs = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                                          na.rm = TRUE),
                                 labels = c("1-Least deprivation", "2-BelowAvg deprivation",
                                            "3-Average deprivation","4-AboveAvg deprivation",
                                            "5-Most deprivation"),
                                 include.lowest = TRUE),
                  NDIQuint = factor(replace(as.character(NDIQuint), 
                                            is.na(NDIQuint) | is.infinite(NDIQuint),
                                            "9-NDI not avail"),
                                    c(levels(NDIQuint), "9-NDI not avail"))) %>%
    dplyr::select(NDI, NDIQuint)
  
  if (is.null(df)) {
    # Format output
    if (round_output == TRUE) {
      ndi <- cbind(ndi_vars, NDIQuint) %>%
        dplyr::mutate(PctRecvIDR = round(PctRecvIDR, digits = 1),
                      PctPubAsst = round(PctPubAsst, digits = 1),
                      PctMgmtBusScArti = round(PctMgmtBusScArti, digits = 1),
                      PctFemHeadKids = round(PctFemHeadKids, digits = 1),
                      PctOwnerOcc = round(PctOwnerOcc, digits = 1),
                      PctNoPhone = round(PctNoPhone, digits = 1),
                      PctNComPlmb = round(PctNComPlmb, digits = 1),
                      PctEducHSPlus = round(PctEducHSPlus, digits = 1),
                      PctEducBchPlus = round(PctEducBchPlus, digits = 1),
                      PctFamBelowPov = round(PctFamBelowPov, digits = 1),
                      PctUnempl = round(PctUnempl, digits = 1))
    } else {
      ndi <- cbind(ndi_vars, NDIQuint)
    }
    
    if (geo == "tract") {
      ndi <- ndi %>%
        dplyr::select(GEOID, 
                      state,
                      county,
                      tract,
                      NDI, NDIQuint, 
                      MedHHInc, PctRecvIDR, PctPubAsst, MedHomeVal, PctMgmtBusScArti,
                      PctFemHeadKids,PctOwnerOcc, PctNoPhone, PctNComPlmb, PctEducHSPlus,
                      PctEducBchPlus, PctFamBelowPov, PctUnempl, TotalPop) 
    } else {
      ndi <- ndi %>%
        dplyr::select(GEOID,
                      state,
                      county,
                      NDI, NDIQuint,
                      MedHHInc, PctRecvIDR, PctPubAsst, MedHomeVal, PctMgmtBusScArti,
                      PctFemHeadKids,PctOwnerOcc, PctNoPhone, PctNComPlmb, PctEducHSPlus,
                      PctEducBchPlus, PctFamBelowPov, PctUnempl, TotalPop) 
    }
    
    ndi <- ndi %>%
      dplyr::mutate(state = stringr::str_trim(state),
                    county = stringr::str_trim(county)) %>%
      dplyr::arrange(GEOID) %>%
      dplyr::as_tibble()
    
  } else {
    ndi <- cbind(df[ , 1], NDIQuint, df[ , 2:ncol(df)])
    ndi <- dplyr::as_tibble(ndi[order(ndi[ , 1]), ])
  }
  
  out <- list(ndi = ndi,
              pca = fit_rotate,
              missing = missingYN,
              cronbach = crnbch)
  
  return(out)
}
