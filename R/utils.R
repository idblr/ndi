# Internal function for Dissimilarity Index (James & Taeuber)
## Returns NA value if only one smaller geography in a larger geography
di_fun <- function(x) {
  if (nrow(x) < 2) {
    NA
  } else {
    1/2*sum(abs(x$subgroup/sum(x$subgroup, na.rm = TRUE) - x$subgroup_ref/sum(x$subgroup_ref, na.rm = TRUE)))
  }
}
