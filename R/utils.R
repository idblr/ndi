# Internal function for the Dissimilarity Index (Duncan & Duncan 1955)
## Returns NA value if only one smaller geography in a larger geography
di_fun <- function(x, omit_NAs) {
  xx <- x[ , c('subgroup', 'subgroup_ref')]
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(x) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    0.5 * sum(
      abs(
        xx$subgroup / sum(xx$subgroup, na.rm = TRUE) - 
          xx$subgroup_ref / sum(xx$subgroup_ref, na.rm = TRUE)
      ), 
      na.rm = TRUE)
  }
}

# Internal function for the Atkinson Index (Atkinson 1970)
## Returns NA value if only one smaller geography in a larger geography
## If denoting the Hölder mean
ai_fun <- function(x, epsilon, omit_NAs) {
  if (omit_NAs == TRUE) { 
    xx <- stats::na.omit(x$subgroup)
  } else {
    xx <- x$subgroup
  } 
  if (nrow(x) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    if (epsilon == 1) {
      1 - (exp(mean(log(stats::na.omit(xx)))) / mean(xx, na.rm = TRUE))
    } else {
      xxx <- (xx / mean(xx, na.rm = TRUE)) ^ (1 - epsilon)
      1 - mean(xxx, na.rm = TRUE) ^ (1 / (1 - epsilon)) 
    }
  }
}

# Internal function for the aspatial Racial Isolation Index (Bell 1954)
## Returns NA value if only one smaller geography in a larger geography
ii_fun <- function(x, omit_NAs) {
  xx <- x[ , c('TotalPopE', 'subgroup', 'subgroup_ixn')]
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(x) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    sum(
      (xx$subgroup / sum(xx$subgroup, na.rm = TRUE)) * (xx$subgroup_ixn / xx$TotalPopE),
      na.rm = TRUE
    )
  }
}

# Internal function for the aspatial Correlation Ratio (White 1986)
## Returns NA value if only one smaller geography in a larger geography
v_fun <- function(x, omit_NAs) {
  xx <- x[ , c('TotalPopE', 'subgroup')]
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(x) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    xxx <- sum(
      (xx$subgroup / sum(xx$subgroup, na.rm = TRUE)) * (xx$subgroup / xx$TotalPopE),
      na.rm = TRUE
    )
    px <- sum(xx$subgroup, na.rm = TRUE) / sum(xx$TotalPopE, na.rm = TRUE)
    (xxx - px) / (1 - px)
  }
}

# Internal function for the aspatial Location Quotient (Sudano et al. 2013)
## Returns NA value if only one smaller geography in a larger geography
lq_fun <- function(x, omit_NAs) {
  xx <- x[ , c('TotalPopE', 'subgroup', 'GEOID')]
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(x) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    p_im <- xx$subgroup / xx$TotalPopE
    if (anyNA(p_im)) { p_im[is.na(p_im), ] <- 0 }
    LQ <- p_im / (sum(xx$subgroup, na.rm = TRUE) / sum(xx$TotalPopE, na.rm = TRUE))
    df <-  data.frame(LQ = LQ, GEOID = xx$GEOID)
    return(df)
  }
}

# Internal function for the aspatial Local Exposure & Isolation (Bemanian & Beyer 2017) metric
## Returns NA value if only one smaller geography in a larger geography
lexis_fun <- function(x, omit_NAs) {
  xx <- x[ , c('TotalPopE', 'subgroup', 'subgroup_ixn', 'GEOID')]
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(x) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    p_im <- xx$subgroup / xx$TotalPopE
    if (anyNA(p_im)) { p_im[is.na(p_im), ] <- 0 }
    p_in <- xx$subgroup_ixn / xx$TotalPopE
    if (anyNA(p_in)) { p_in[is.na(p_in), ] <- 0 }
    P_m <- sum(xx$subgroup, na.rm = TRUE) / sum(xx$TotalPopE, na.rm = TRUE)
    P_n <- sum(xx$subgroup_ixn, na.rm = TRUE) / sum(xx$TotalPopE, na.rm = TRUE)
    LExIs <- car::logit(p_im * p_in) - car::logit(P_m * P_n)
    df <-  data.frame(LExIs = LExIs, GEOID = xx$GEOID)
    return(df)
  }
}

# Internal function for the aspatial Delta (Hoover 1941)
## Returns NA value if only one smaller geography in a larger geography
del_fun <- function(x, omit_NAs) {
  xx <- x[ , c('subgroup', 'ALAND')]
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(x) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    0.5 * sum(
      abs((xx$subgroup / sum(xx$subgroup, na.rm = TRUE)) - (xx$ALAND / sum(xx$ALAND, na.rm = TRUE))
      ),
      na.rm = TRUE
    )
  }
}
