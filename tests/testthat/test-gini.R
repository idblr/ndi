context("gini")

###################
# gini testthat #
###################

test_that("gini throws error with invalid arguments", {
  
  # Unavailable geography
  expect_error(gini(geo = "zcta", state = "DC", year = 2020, quiet = TRUE))
  
  # Unavailable year
  expect_error(gini(state = "DC", year = 2005, quiet = TRUE))
  
  skip_if(Sys.getenv("CENSUS_API_KEY") == "")
  
  # Incorrect state
  expect_error(gini(state = "AB", year = 2020))
  
  # Unavailable geography for DC (only 1 'county' in DC so, alone, NDI cannot be computed)
  expect_error(gini(geo = "county", state = "DC", year = 2009, quiet = TRUE))
  
}
)   

test_that("gini works", {  
  
  skip_if(Sys.getenv("CENSUS_API_KEY") == "")
  
  expect_message(gini(state = "DC", year = 2020)) 
  
  expect_silent(gini(state = "DC", year = 2020, quiet = TRUE))
  
}
)  
