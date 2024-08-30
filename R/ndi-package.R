#' The ndi Package: Neighborhood Deprivation Indices
#'
#' Computes various geospatial indices of socioeconomic deprivation and disparity in the United States based on information available from the U.S. Census Bureau.
#'
#' @details The 'ndi' package computes various indices of socioeconomic deprivation and disparity in the United States. Some indices are considered "spatial" because they consider the values of neighboring (i.e., adjacent) census geographies in their computation, while other indices are "aspatial" because they only consider the value within each census geography. Two types of aspatial neighborhood deprivation indices (\emph{NDI}) are available: (1) based on Messer et al. (2006) \doi{10.1007/s11524-006-9094-x} and (2) based on Andrews et al. (2020) \doi{10.1080/17445647.2020.1750066} and Slotman et al. (2022) \doi{10.1016/j.dib.2022.108002} who use variables chosen by Roux and Mair (2010) \doi{10.1111/j.1749-6632.2009.05333.x}. Both are a decomposition of multiple demographic characteristics from the U.S. Census Bureau American Community Survey 5-year estimates (ACS-5; 2006-2010 onward). Using data from the ACS-5 (2005-2009 onward), the package can also compute indices of racial or ethnic residential segregation, including but limited to those discussed in Massey & Denton (1988) \doi{10.1093/sf/67.2.281}, and additional indices of socioeconomic disparity.
#' 
#' Key content of the 'ndi' package include:\cr
#' 
#' \strong{Neighborhood Deprivation Indices}
#' 
#' \code{\link{messer}} Computes the aspatial Neighborhood Deprivation Index (\emph{NDI}) based on Messer et al. (2006) \doi{10.1007/s11524-006-9094-x}.
#' 
#' \code{\link{powell_wiley}} Computes the aspatial Neighborhood Deprivation Index (\emph{NDI}) based on Andrews et al. (2020) \doi{10.1080/17445647.2020.1750066} and Slotman et al. (2022) \doi{10.1016/j.dib.2022.108002} who use variables chosen by Roux and Mair (2010) \doi{10.1111/j.1749-6632.2009.05333.x}.
#' 
#' \strong{Indices of Racial or Ethnic Residential Segregation}
#' 
#' \code{\link{anthopolos}} Computes the spatial Racial Isolation Index (\emph{RI}) based on Anthopolos (2011) \doi{10.1016/j.sste.2011.06.002}.
#' 
#' \code{\link{atkinson}} Computes the aspatial Atkinson Index (\emph{A}) based on Atkinson (1970) \doi{10.1016/0022-0531(70)90039-6}.
#' 
#' \code{\link{bell}} Computes the aspatial Interaction Index (\emph{xPy\*}) based on Shevky & Williams (1949; ISBN-13:978-0-837-15637-8) and Bell (1954) \doi{10.2307/2574118}.
#' 
#' \code{\link{bemanian_beyer}} Computes the aspatial Local Exposure and Isolation (\emph{LEx/Is}) based on Bemanian & Beyer (2017) \doi{10.1158/1055-9965.EPI-16-0926}.
#' 
#' \code{\link{denton}} Computes the aspatial Relative Clustering (\emph{RCL}) based on Massey & Denton (1988) \doi{10.1093/sf/67.2.281}.
#'
#' \code{\link{duncan}} Computes the aspatial Dissimilarity Index (\emph{D}) based on Duncan & Duncan (1955a) \doi{10.2307/2088328}.
#' 
#' \code{\link{duncan_cuzzort}} Computes the aspatial Absolute Centralization (\emph{ACE}) based on Duncan, Cuzzort, & Duncan (1961; LC:60007089) and Massey & Denton (1988) \doi{10.1093/sf/67.2.281}.
#'
#' \code{\link{duncan_duncan}} Computes the aspatial Relative Centralization (\emph{RCE}) based on Duncan & Duncan (1955b) \doi{10.1086/221609} and Massey & Denton (1988) \doi{10.1093/sf/67.2.281}.
#'
#' \code{\link{gini}} Computes the aspatial Gini Index (\emph{G}) based on Gini (1921) \doi{10.2307/2223319}.
#' 
#' \code{\link{hoover}} Computes the aspatial Delta (\emph{DEL}) based on Hoover (1941) \doi{doi:10.1017/S0022050700052980} and Duncan, Cuzzort, & Duncan (1961; LC:60007089).
#'
#' \code{\link{james_taeuber}} Computes the aspatial Dissimilarity Index (\emph{D}) based on James & Taeuber (1985) \doi{10.2307/270845}.
#' 
#' \code{\link{krieger}} Computes the aspatial Index of Concentration at the Extremes based on Feldman et al. (2015) \doi{10.1136/jech-2015-205728} and Krieger et al. (2016) \doi{10.2105/AJPH.2015.302955}.
#' 
#' \code{\link{lieberson}} Computes the aspatial Isolation Index (\emph{xPx\*}) based on Lieberson (1981; ISBN-13:978-1-032-53884-6) and Bell (1954) \doi{10.2307/2574118}.
#' 
#' \code{\link{massey}} Computes the aspatial Absolute Clustering (\emph{ACL}) based on Massey & Denton (1988) \doi{10.1093/sf/67.2.281}.
#'
#' \code{\link{sudano}} Computes the aspatial Location Quotient (\emph{LQ}) based on Merton (1939) \doi{10.2307/2084686} and Sudano et al. (2013) \doi{10.1016/j.healthplace.2012.09.015}.
#' 
#' \code{\link{theil}} Computes the aspatial Entropy (\emph{H}) based on Theil (1972; ISBN-13:978-0-444-10378-9) and Theil & Finizza (1971) \doi{110.1080/0022250X.1971.9989795}.
#' 
#' \code{\link{white}} Computes the aspatial Correlation Ratio (\emph{V}) based on Bell (1954) \doi{10.2307/2574118} and White (1986) \doi{10.2307/3644339}.
#' 
#' \code{\link{white_blau}} Computes an index of spatial proximity (\emph{SP}) based on White (1986) \doi{10.2307/3644339} and Blau (1977; ISBN-13:978-0-029-03660-0).
#' 
#' \strong{Additional Indices of Socioeconomic Disparity}
#' 
#' \code{\link{atkinson}} Also computes the aspatial Atkinson Index (\emph{A}) of income based on Atkinson (1970) \doi{10.1016/0022-0531(70)90039-6}.
#' 
#' \code{\link{bravo}} Computes the spatial Educational Isolation Index (\emph{EI}) based on Bravo (2021) \doi{10.3390/ijerph18179384}.
#' 
#' \code{\link{gini}} Also retrieves the aspatial Gini Index (\emph{G}) of income inequality based on Gini (1921) \doi{10.2307/2223319}.
#' 
#' \strong{Pre-formatted U.S. Census Data}
#' 
#' \code{\link{DCtracts2020}} A sample dataset containing information about U.S. Census American Community Survey 5-year estimate data for the District of Columbia census tracts (2020). The data are obtained from the \code{\link[tidycensus]{get_acs}} function and formatted for the \code{\link{messer}} and \code{\link{powell_wiley}} functions input.
#' 
#' @name ndi-package
#' @aliases ndi-package ndi
#' 
#' @section Dependencies: The 'ndi' package relies heavily upon \code{\link{tidycensus}} to retrieve data from the U.S. Census Bureau American Community Survey five-year estimates and the \code{\link{psych}} for computing the neighborhood deprivation indices. The \code{\link{messer}} function builds upon code developed by Hruska et al. (2022) \doi{10.17605/OSF.IO/M2SAV} by fictionalizing, adding the percent of households earning <$30,000 per year to the \emph{NDI} computation, and providing the option for computing the ACS-5 2006-2010 \emph{NDI} values. There is no code companion to compute \emph{NDI} included in Andrews et al. (2020) \doi{10.1080/17445647.2020.1750066} or Slotman et al. (2022) \doi{10.1016/j.dib.2022.108002}, but the package author worked directly with the Slotman et al. (2022) \doi{10.1016/j.dib.2022.108002} authors to replicate their SAS code in \strong{R}. The indices of racial or ethnic residential segregation rely heavily on the \code{\link{sf}} and \code{\link{tigris}} packages to assign the smaller geographical units within larger geographical units and, occasionally, perform geospatial projection for distance-based metrics. The computation of \emph{RI} and \emph{EI} also relies on the \code{\link{Matrix}} package to compute the geospatial adjacency matrix between census geographies. Internal function to calculate \emph{AI} using the Hölder mean is based on \code{\link[DescTools]{Atkinson}} function.
#' 
#' @author Ian D. Buller\cr \emph{DLH, LLC (formerly DLH Corporation and Social & Scientific Systems, Inc.), Bethesda, Maryland, USA (current); Occupational and Environmental Epidemiology Branch, Division of Cancer Epidemiology and Genetics, National Cancer Institute, National Institutes of Health, Rockville, Maryland, USA (original).} \cr
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
#' @importFrom sf st_drop_geometry st_geometry st_intersects st_transform st_within
#' @importFrom stats complete.cases cor cov2cor loadings median na.omit promax quantile sd setNames
#' @importFrom stringr str_trim
#' @importFrom tidycensus get_acs
#' @importFrom tidyr pivot_longer separate
#' @importFrom tigris combined_statistical_areas core_based_statistical_areas counties metro_divisions places states tracts
#' @importFrom units drop_units set_units
#' @importFrom utils stack
NULL
