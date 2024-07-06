#' Formatted U.S. Census American Community Survey 5-year estimate data for DC census tracts (2020) from the 'tidycensus' package
#'
#' A sample dataset containing information about U.S. Census American Community Survey 5-year estimate data for the District of Columbia census tracts (2020). The data are obtained from the \code{\link[tidycensus]{get_acs}} function and formatted for the \code{\link{messer}} and \code{\link{powell_wiley}} functions input.
#'
#' @format A data frame with 206 rows and 23 variables:
#' \describe{
#'   \item{GEOID}{census tract ID}
#'   \item{TotalPop}{arcsinh-transformed CD3}
#'   \item{OCC}{percent males in management, science, and arts occupation}
#'   \item{CWD}{percent of crowded housing}
#'   \item{POV}{percent of households in poverty}
#'   \item{FHH}{percent of female headed households with dependents}
#'   \item{PUB}{percent of households on public assistance}
#'   \item{U30}{percent of households earning <$30,000 per year}
#'   \item{EDU}{percent earning less than a high school education}
#'   \item{EMP}{percent unemployed}
#'   \item{logMedHHInc}{median household income (dollars), natural log-transformed}
#'   \item{PctNoIDRZ}{percent of households receiving dividends, interest, or rental income, Z-transformed}
#'   \item{PctPubAsstZ}{percent of households receiving public assistance, Z-transformed}
#'   \item{logMedHomeVal}{median home value (dollars), natural log-transformed}
#'   \item{PctWorkClassZ}{percent in a management, business, science, or arts occupation, Z-transformed}
#'   \item{PctFemHeadKidsZ}{percent of households that are female headed with any children under 18 years, Z-transformed}
#'   \item{PctNotOwnerOccZ}{percent of housing units that are owner occupied, Z-transformed}
#'   \item{PctNoPhoneZ}{percent of households without a telephone, Z-transformed}
#'   \item{PctNComPlmbZ}{percent of households without complete plumbing facilities, Z-transformed}
#'   \item{PctEducLTHSZ}{percent with a high school degree or higher (population 25 years and over), Z-transformed}
#'   \item{PctEducLTBchZ}{percent with a college degree or higher (population 25 years and over), Z-transformed}
#'   \item{PctFamBelowPovZ}{percent of families with incomes below the poverty level, Z-transformed}
#'   \item{PctUnemplZ}{percent unemployed, Z-transformed}
#' }
#' @examples
#' head(DCtracts2020)
#'
#' @source \url{https://github.com/idblr/ndi/blob/master/README.md}
'DCtracts2020'
