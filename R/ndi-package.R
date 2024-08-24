#' The ndi Package: Neighborhood Deprivation Indices
#'
#' Computes various metrics of socio-economic deprivation and disparity in the United States based on information available from the U.S. Census Bureau.
#'
#' @details The 'ndi' package computes various metrics of socio-economic deprivation and disparity in the United States. Some metrics are considered "spatial" because they consider the values of neighboring (i.e., adjacent) census geographies in their computation, while other metrics are "aspatial" because they only consider the value within each census geography. Two types of aspatial neighborhood deprivation indices (\emph{NDI}) are available: (1) based on Messer et al. (2006) \doi{10.1007/s11524-006-9094-x} and (2) based on Andrews et al. (2020) \doi{10.1080/17445647.2020.1750066} and Slotman et al. (2022) \doi{10.1016/j.dib.2022.108002} who use variables chosen by Roux and Mair (2010) \doi{10.1111/j.1749-6632.2009.05333.x}. Both are a decomposition of multiple demographic characteristics from the U.S. Census Bureau American Community Survey 5-year estimates (ACS-5; 2006-2010 onward). Using data from the ACS-5 (2005-2009 onward), the package can also compute the (1) spatial Racial Isolation Index (\emph{RI}) based on Anthopolos et al. (2011) \doi{10.1016/j.sste.2011.06.002}, (2) spatial Educational Isolation Index (\emph{EI}) based on Bravo et al. (2021) \doi{10.3390/ijerph18179384}, (3) aspatial Index of Concentration at the Extremes (\emph{ICE}) based on Feldman et al. (2015) \doi{10.1136/jech-2015-205728} and Krieger et al. (2016) \doi{10.2105/AJPH.2015.302955}, (4) aspatial racial or ethnic Dissimilarity Index (\emph{DI}) based on Duncan & Duncan (1955) \doi{10.2307/2088328}, (5) aspatial income or racial or ethnic Atkinson Index (\emph{AI}) based on Atkinson (1970) \doi{10.1016/0022-0531(70)90039-6}, (6) aspatial racial or ethnic Isolation Index (\emph{II}) based on Shevky & Williams (1949; ISBN-13:978-0-837-15637-8) and Bell (1954) \doi{10.2307/2574118}, (7) aspatial racial or ethnic Correlation Ratio (\emph{V}) based on Bell (1954) \doi{10.2307/2574118} and White (1986) \doi{10.2307/3644339}, (8) aspatial racial or ethnic Location Quotient (\emph{LQ}) based on Merton (1939) \doi{10.2307/2084686} and Sudano et al. (2013) \doi{10.1016/j.healthplace.2012.09.015}, (9) aspatial racial or ethnic Local Exposure and Isolation (\emph{LEx/Is}) metric based on Bemanian & Beyer (2017) \doi{10.1158/1055-9965.EPI-16-0926}, (10) aspatial racial or ethnic Delta (\emph{DEL}) based on Hoover (1941) \doi{10.1017/S0022050700052980} and Duncan et al. (1961; LC:60007089), and (11) an index of spatial proximity (\emph{SP}) based on White (1986) \doi{10.2307/3644339} and Blau (1977; ISBN-13:978-0-029-03660-0). Also using data from the ACS-5 (2005-2009 onward), the package can retrieve the aspatial Gini Index (\emph{G}) based on Gini (1921) \doi{10.2307/2223319}.
#' 
#' Key content of the 'ndi' package include:\cr
#' 
#' \bold{Metrics of Socio-Economic Deprivation and Disparity}
#' 
#' \code{\link{anthopolos}} Computes the spatial Racial Isolation Index (\emph{RI}) based on Anthopolos (2011) \doi{10.1016/j.sste.2011.06.002}.
#' 
#' \code{\link{atkinson}} Computes the aspatial income or racial or ethnic Atkinson Index (\emph{A}) based on Atkinson (1970) \doi{10.1016/0022-0531(70)90039-6}.
#' 
#' \code{\link{bell}} Computes the aspatial racial or ethnic Interaction Index (\emph{xPy\*}) based on Shevky & Williams (1949; ISBN-13:978-0-837-15637-8) and Bell (1954) \doi{10.2307/2574118}.
#' 
#' \code{\link{bemanian_beyer}} Computes the aspatial racial or ethnic Local Exposure and Isolation (\emph{LEx/Is}) metric based on Bemanian & Beyer (2017) \doi{10.1158/1055-9965.EPI-16-0926}.
#' 
#' \code{\link{bravo}} Computes the spatial Educational Isolation Index (\emph{EI}) based on Bravo (2021) \doi{10.3390/ijerph18179384}.
#' 
#' \code{\link{duncan}} Computes the aspatial racial or ethnic Dissimilarity Index (\emph{D}) based on Duncan & Duncan (1955) \doi{10.2307/2088328}.
#' 
#' \code{\link{gini}} Computes the aspatial Gini Index (\emph{G}) of racial or ethnic inequality and retrieves the aspatial Gini Index (\emph{G}) of income inequality based on Gini (1921) \doi{10.2307/2223319}.
#' 
#' \code{\link{hoover}} Computes the aspatial racial or ethnic Delta (\emph{DEL}) based on Hoover (1941) \doi{doi:10.1017/S0022050700052980} and Duncan et al. (1961; LC:60007089).
#'
#' \code{\link{james_taeuber}} Computes the aspatial racial or ethnic Dissimilarity Index (\emph{D}) based on James & Taeuber (1985) \doi{10.2307/270845}.
#' 
#' \code{\link{krieger}} Computes the aspatial Index of Concentration at the Extremes based on Feldman et al. (2015) \doi{10.1136/jech-2015-205728} and Krieger et al. (2016) \doi{10.2105/AJPH.2015.302955}.
#' 
#' \code{\link{lieberson}} Computes the aspatial racial or ethnic Isolation Index (\emph{xPx\*}) based on Lieberson (1981; ISBN-13:978-1-032-53884-6) and Bell (1954) \doi{10.2307/2574118}.
#' 
#' \code{\link{messer}} Computes the aspatial Neighborhood Deprivation Index (\emph{NDI}) based on Messer et al. (2006) \doi{10.1007/s11524-006-9094-x}.
#' 
#' \code{\link{powell_wiley}} Computes the aspatial Neighborhood Deprivation Index (\emph{NDI}) based on Andrews et al. (2020) \doi{10.1080/17445647.2020.1750066} and Slotman et al. (2022) \doi{10.1016/j.dib.2022.108002} who use variables chosen by Roux and Mair (2010) \doi{10.1111/j.1749-6632.2009.05333.x}.
#' 
#' \code{\link{sudano}} Computes the aspatial racial or ethnic Location Quotient (\emph{LQ}) based on Merton (1939) \doi{10.2307/2084686} and Sudano et al. (2013) \doi{10.1016/j.healthplace.2012.09.015}.
#' 
#' \code{\link{white}} Computes the aspatial racial or ethnic Correlation Ratio (\emph{V}) based on Bell (1954) \doi{10.2307/2574118} and White (1986) \doi{10.2307/3644339}.
#' 
#' \code{\link{white_blau}} Computes an index of spatial proximity (\emph{SP}) based on White (1986) \doi{10.2307/3644339} and Blau (1977; ISBN-13:978-0-029-03660-0).
#' 
#' \bold{Pre-formatted U.S. Census Data}
#' 
#' \code{\link{DCtracts2020}} A sample dataset containing information about U.S. Census American Community Survey 5-year estimate data for the District of Columbia census tracts (2020). The data are obtained from the \code{\link[tidycensus]{get_acs}} function and formatted for the \code{\link{messer}} and \code{\link{powell_wiley}} functions input.
#' 
#' @name ndi-package
#' @aliases ndi-package ndi
#' 
#' @section Dependencies: The 'ndi' package relies heavily upon \code{\link{tidycensus}} to retrieve data from the U.S. Census Bureau American Community Survey five-year estimates and the \code{\link{psych}} for computing the neighborhood deprivation indices. The \code{\link{messer}} function builds upon code developed by Hruska et al. (2022) \doi{10.17605/OSF.IO/M2SAV} by fictionalizing, adding the percent of households earning <$30,000 per year to the NDI computation, and providing the option for computing the ACS-5 2006-2010 NDI values. There is no code companion to compute NDI included in Andrews et al. (2020) \doi{10.1080/17445647.2020.1750066} or Slotman et al. (2022) \doi{10.1016/j.dib.2022.108002}, but the package author worked directly with the Slotman et al. (2022) \doi{10.1016/j.dib.2022.108002} authors to replicate their SAS code in R. The spatial metrics RI and EI rely on the \code{\link{sf}} and \code{\link{Matrix}} packages to compute the geospatial adjacency matrix between census geographies. Internal function to calculate AI is based on \code{\link[DescTools]{Atkinson}} function. There is no code companion to compute RI, EI, DI, II, V, LQ, or LEx/Is included in Anthopolos et al. (2011) \doi{10.1016/j.sste.2011.06.002}, Bravo et al. (2021) \doi{10.3390/ijerph18179384}, Duncan & Duncan (1955) \doi{10.2307/2088328}, Bell (1954) \doi{10.2307/2574118}, White (1986) \doi{10.2307/3644339}, Sudano et al. (2013) \doi{10.1016/j.healthplace.2012.09.015}, or Bemanian & Beyer (2017) \doi{10.1158/1055-9965.EPI-16-0926}, respectively.
#' 
#' @author Ian D. Buller\cr \emph{DLH Corporation (formerly Social & Scientific Systems, Inc.), Bethesda, Maryland, USA (current); Occupational and Environmental Epidemiology Branch, Division of Cancer Epidemiology and Genetics, National Cancer Institute, National Institutes of Health, Rockville, Maryland, USA (original).} \cr
#' 
#' Maintainer: I.D.B. \email{ian.buller@@alumni.emory.edu}
#'
#' @keywords internal
'_PACKAGE'

#' @import dplyr
#' @importFrom car logit
#' @importFrom MASS ginv
#' @importFrom Matrix sparseMatrix
#' @importFrom psych alpha principal
#' @importFrom sf st_drop_geometry st_geometry st_intersects st_within
#' @importFrom stats complete.cases cor cov2cor loadings median na.omit promax quantile sd setNames
#' @importFrom stringr str_trim
#' @importFrom tidycensus get_acs
#' @importFrom tidyr pivot_longer separate
#' @importFrom tigris combined_statistical_areas core_based_statistical_areas metro_divisions
#' @importFrom units drop_units set_units
#' @importFrom utils stack
NULL
