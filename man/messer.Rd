% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/messer.R
\name{messer}
\alias{messer}
\title{Neighborhood Deprivation Index based on Messer et al. (2006)}
\usage{
messer(
  geo = "tract",
  year = 2020,
  imp = FALSE,
  quiet = FALSE,
  round_output = FALSE,
  df = NULL,
  ...
)
}
\arguments{
\item{geo}{Character string specifying the geography of the data either census tracts \code{geo = 'tract'} (the default) or counties \code{geo = 'county'}.}

\item{year}{Numeric. The year to compute the estimate. The default is 2020, and the years 2010 onward are currently available.}

\item{imp}{Logical. If TRUE, will impute missing census characteristics within the internal \code{\link[psych]{principal}}. If FALSE (the default), will not impute.}

\item{quiet}{Logical. If TRUE, will display messages about potential missing census information and the proportion of variance explained by principal component analysis. The default is FALSE.}

\item{round_output}{Logical. If TRUE, will round the output of raw census and \emph{NDI} values from the \code{\link[tidycensus]{get_acs}} at one and four significant digits, respectively. The default is FALSE.}

\item{df}{Optional. Pass a pre-formatted \code{'dataframe'} or \code{'tibble'} with the desired variables through the function. Bypasses the data obtained by \code{\link[tidycensus]{get_acs}}. The default is NULL. See Details below.}

\item{...}{Arguments passed to \code{\link[tidycensus]{get_acs}} to select state, county, and other arguments for census characteristics}
}
\value{
An object of class 'list'. This is a named list with the following components:

\describe{
\item{\code{ndi}}{An object of class 'tbl' for the GEOID, name, \emph{NDI} (standardized), \emph{NDI} (quartile), and raw census values of specified census geographies.}
\item{\code{pca}}{An object of class 'principal', returns the output of \code{\link[psych]{principal}} used to compute the \emph{NDI} values.}
\item{\code{missing}}{An object of class 'tbl' of the count and proportion of missingness for each census variable used to compute \emph{NDI}.}
}
}
\description{
Compute the aspatial Neighborhood Deprivation Index (Messer).
}
\details{
This function will compute the aspatial Neighborhood Deprivation Index (\emph{NDI}) of U.S. census tracts or counties for a specified geographical referent (e.g., US-standardized) based on Messer et al. (2006) \doi{10.1007/s11524-006-9094-x}.

The function uses the \code{\link[tidycensus]{get_acs}} function to obtain U.S. Census Bureau 5-year American Community Survey characteristics used for computation involving a principal component analysis with the \code{\link[psych]{principal}} function. The yearly estimates are available for 2010 and after when all census characteristics became available. The eight characteristics are:
\itemize{
\item \strong{OCC (C24030)}: percent males in management, science, and arts occupation
\item \strong{CWD (B25014)}: percent of crowded housing
\item \strong{POV (B17017)}: percent of households in poverty
\item \strong{FHH (B25115)}: percent of female headed households with dependents
\item \strong{PUB (B19058)}: percent of households on public assistance
\item \strong{U30 (B19001)}: percent of households earning <$30,000 per year
\item \strong{EDU (B06009)}: percent earning less than a high school education
\item \strong{EMP (B23025)}: percent unemployed (2011 onward)
\item \strong{EMP (B23001)}: percent unemployed (2010 only)
}

Use the internal \code{state} and \code{county} arguments within the \code{\link[tidycensus]{get_acs}} function to specify the referent for standardizing the \emph{NDI} (Messer) values. For example, if all U.S. states are specified for the \code{state} argument, then the output would be a U.S.-standardized index.

The continuous \emph{NDI} (Messer) values are z-transformed, i.e., 'standardized,' and the categorical \emph{NDI} (Messer) values are quartiles of the standardized continuous \emph{NDI} (Messer) values.

Check if the proportion of variance explained by the first principal component is high (more than 0.5).

Users can bypass \code{\link[tidycensus]{get_acs}} by specifying a pre-formatted data frame or tibble using the \code{df} argument. This function will compute an index using the first component of a principal component analysis (PCA) with a Varimax rotation (the default for \code{\link[psych]{principal}}) and only one factor (note: PCA set-up not unspecified in Messer et al. (2006)). The recommended structure of the data frame or tibble is an ID (e.g., GEOID) in the first feature (column), followed by the variables of interest (in any order) and no additional information (e.g., omit state or county names from the \code{df} argument input).
}
\examples{

messer(df = DCtracts2020[ , c(1, 3:10)])

\dontrun{
# Wrapped in \dontrun{} because these examples require a Census API key.
  
  # Tract-level NDI (Messer; 2020)
  messer(geo = 'tract', state = 'GA', year = 2020)

  # Impute NDI (Messer; 2020) for tracts with missing census information (median values)
  messer(state = 'tract', state = 'GA', year = 2020, imp = TRUE)
  
}

}
\references{
Messer, LC, Laraia, BA, Kaufman, JS, Eyster, J, Holzman, C, Culhane, J, Elo, I, Burke, J, O'Campo, P (2006) The Development of a Standardized Neighborhood Deprivation Index. \emph{Journal of Urban Health}, 83(6):1041-1062. \doi{10.1007/s11524-006-9094-x}
}
\seealso{
\code{\link[tidycensus]{get_acs}} for additional arguments for geographic referent selection (i.e., \code{state} and \code{county}).

Neighborhood Deprivation Index: \code{\link{powell_wiley}}
}
