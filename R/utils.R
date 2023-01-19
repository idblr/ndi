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
