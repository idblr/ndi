# Internal function for Dissimilarity Index (Duncan & Duncan 1955)
## Returns NA value if only one smaller geography in a larger geography
di_fun <- function(x, omit_NAs) {
  xx <- x[ , c("subgroup", "subgroup_ref")]
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(x) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    1/2 * sum(abs(xx$subgroup / sum(xx$subgroup, na.rm = TRUE) - xx$subgroup_ref / sum(xx$subgroup_ref, na.rm = TRUE)))
  }
}

# Internal function for Atkinson Index (Atkinson 1970)
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

# Internal function for aspatial Racial Isolation Index (Bell 1954)
## Returns NA value if only one smaller geography in a larger geography
ii_fun <- function(x, omit_NAs) {
  xx <- x[ , c("TotalPopE", "subgroup", "subgroup_ixn")]
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(x) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    sum((xx$subgroup / sum(xx$subgroup, na.rm = TRUE)) * (xx$subgroup_ixn / xx$TotalPopE))
  }
}

# Internal function for aspatial Correlation Ratio (White 1986)
## Returns NA value if only one smaller geography in a larger geography
v_fun <- function(x, omit_NAs) {
  xx <- x[ , c("TotalPopE", "subgroup")]
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(x) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    xxx <- sum((xx$subgroup / sum(xx$subgroup, na.rm = TRUE)) * (xx$subgroup / xx$TotalPopE))
    px <- sum(xx$subgroup, na.rm = TRUE) / sum(xx$TotalPopE, na.rm = TRUE)
    (xxx - px) / (1 - px)
  }
}

# Internal function for aspatial Location Quotient (Sudano et al. 2013)
## Returns NA value if only one smaller geography in a larger geography
lq_fun <- function(x, omit_NAs) {
  xx <- x[ , c("TotalPopE", "subgroup", "GEOID")]
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(x) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    LQ <- (xx$subgroup / xx$TotalPopE) / (sum(xx$subgroup, na.rm = TRUE) / sum(xx$TotalPopE, na.rm = TRUE))
    df <-  data.frame(LQ = LQ,
                      GEOID = xx$GEOID)
    return(df)
  }
}

# Internal function for aspatial Local Exposure & Isolation (Bemanian & Beyer 2017) metric
## Returns NA value if only one smaller geography in a larger geography
lexis_fun <- function(x, omit_NAs) {
  xx <- x[ , c("TotalPopE", "subgroup", "subgroup_ixn", "GEOID")]
  if (omit_NAs == TRUE) { xx <- xx[stats::complete.cases(xx), ] }
  if (nrow(x) < 2 || any(xx < 0) || any(is.na(xx))) {
    NA
  } else {
    p_im <- xx$subgroup / xx$TotalPopE
    p_in <- xx$subgroup_ixn / xx$TotalPopE
    P_m <- sum(xx$subgroup, na.rm = TRUE) / sum(xx$TotalPopE, na.rm = TRUE)
    P_n <- sum(xx$subgroup_ixn, na.rm = TRUE) / sum(xx$TotalPopE, na.rm = TRUE)
    
    # Logit transformation and center with assumption of even mixing across large geography
    # (log((p_im * p_in)/(1 - p_im * xx$p_in)) -  log((P_m * P_n)/(1 - P_m * P_n)))
    LExIs <- car::logit(p_im * p_in) - car::logit(P_m * P_n)
    df <-  data.frame(LExIs = LExIs,
                      GEOID = xx$GEOID)
    return(df)
  }
}
